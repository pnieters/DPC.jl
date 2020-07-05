using ADSP, Distributions, DataFrames, VideoIO, ImageFiltering, ProgressMeter, Plots
import LinearAlgebra, ColorTypes

struct VideoBuffer{T,C}
    timestamps::Vector{T}
    frames::Array{C,3}
end
(get_time_indices(v::VideoBuffer{T,C}, (t1,t2)::Tuple{T,T}) where {T,C}) = searchsortedlast(v.timestamps, t1):searchsortedfirst(v.timestamps, t2)
get_time_indices(v,ts)=ts
Base.getindex(v::VideoBuffer, t_range, y_range, x_range) = v.frames[get_time_indices(v, t_range), y_range, x_range]
Base.checkbounds(::Type{Bool}, v::VideoBuffer, t_range, y_range, x_range) = checkbounds(Bool,v.frames, get_time_indices(v, t_range), y_range, x_range)

"""
    load_video(io::VideoIO.AVInput; c=ColorTypes.Gray)

Opens a video and converts all frames into a datafame with columns `timestamp` (end-time of the frames in seconds) and `frame` (2D-Array holding the frame data converted to color type `c`).
"""
function load_video(io::VideoIO.AVInput; c=ColorTypes.Gray, range=(:,:))
    vid = VideoIO.openvideo(io)

    timestamps = Float64[]
    frames = []

    while !eof(vid)
        # get time of the frame, read frame and convert to correct color type
        frame = c.(read(vid)[range...])
        t = gettime(vid)
        push!(frames, reshape(frame, (1,size(frame)...)))
        push!(timestamps, t)
    end

    buff = VideoBuffer(timestamps, cat(frames...; dims=1))
    buff.timestamps .+= mean(diff(buff.timestamps))
  
    close(vid)
    return buff
end


function ImageFiltering.imfilter(buff, filter)
    tnew = Base.return_types(*, (eltype(buff.frames),eltype(ganglion_rf)))[1]
    filtered_frames = Array{tnew, 3}(undef, size(buff.frames)...)
    
    new_buff = VideoBuffer(deepcopy(buff.timestamps), filtered_frames)
    for i ∈ 1:size(buff.frames, 1)
        new_buff.frames[i,:,:] = imfilter(buff.frames[i,:,:], filter)
    end

    return new_buff
end

"""
    frames_to_spikes(buff; transform=identity)

Converts each pixel of a video into a Poisson spike-train.
"""
function frames_to_spikes(buff; transform=identity)
    @assert size(buff.frames,1) > 0 "`df` must not be empty!"
    dims = size(buff.frames)[2:3]

    frame_durations = diff([0.0; buff.timestamps])
    
    # Utility function to draw a Poisson-spike-train within the current frame interval [a,b]
    # according to the neuron-specific firing rates
    draw_spikes(rate, (a,b)) = sort(rand(Uniform(a,b), rand(Poisson(max(0,rate*firing_rate_scale * (b-a))))))

    spikes = [Float64[] for i ∈ 1:prod(dims)]

    for (i,(t1,t2)) ∈ enumerate(zip([0.0; buff.timestamps[1:end-1]], buff.timestamps))
        # For each receptive field location...
        for (j,inp) ∈ enumerate(buff[i,:,:])
            new_spikes = draw_spikes(transform(inp), (t1,t2))
            # ... generate spikes for center- and surround-cells
            append!(spikes[j],new_spikes)
        end
    end
    
    return spikes
end

"""
    sta(buff, spikes::Vector{Tuple{Tuple{Int,Int}, Float64}}, slice)

Calculate the spike-triggered average for a `slice` (t_slice,y_slice,x_slice) in time and space
over all `spikes` (tuples ((y,x),t) of  positions `(y,x)` and times `t`).
"""
function calculate_sta(buff::VideoBuffer{T,C}, spikes, (t_slice,y_slice,x_slice)) where {T,C}
    sta = zeros(length(t_slice), length(y_slice), length(x_slice))

    c = 0
    for ((y,x),t) ∈ spikes
        t0 = searchsortedfirst(buff.timestamps, t)
        if checkbounds(Bool, buff, t_slice.+t0, y_slice.+y, x_slice.+x)
            c += 1
        else
            continue
        end

        sta .+= buff[t_slice.+t0, y_slice.+y, x_slice.+x]
    end

    sta ./= c

    return sta
end

struct GaborPatch
    scale::Float64
    grating::Vector{Float64}
    phase::Float64
    aperture::Float64
end

struct GaussianPatch
    scale::Float64
    precision::Matrix{Float64}
end

function (g::GaussianPatch)(relative_position)
    return g.scale*exp(-relative_position'*g.precision*relative_position)[]
end
(g::GaussianPatch)(dy,dx) = g([dy,dx])

function (g::GaborPatch)(relative_position)
    scale = g.scale*exp(-sum(relative_position.^2)/(2*g.aperture^2))
    offset = 1.0/g.aperture * g.grating' * relative_position + g.phase
    return scale*cos(π*offset)
end
(g::GaborPatch)(dy,dx) = g([dy,dx])

function configure_rfs!(net, configuration, input_positions)
    function configure_rf!(segment, location, g; connection_threshold=0.0)
        # filter and weight relevant sources for synaptic connections to segment
        for (source_id,position) ∈ input_positions
            synaptic_strength = g(position .- location)
            synapseID = SynapseID(segment.id,source_id,1)
            if synaptic_strength < -connection_threshold && startswith(String(source_id), "off")
                segment.synapses[source_id] = [Synapse(synapseID,-synaptic_strength,false)]
            elseif synaptic_strength > connection_threshold && startswith(String(source_id), "on")
                segment.synapses[source_id] = [Synapse(synapseID, synaptic_strength,false)]
            end
        end
    end

    for (nid,neuron_cfg) ∈ configuration
        nid=Symbol(nid)
        neuron = net.neurons[nid]
        for (bid,branch_cfg) ∈ neuron_cfg
            bid=Symbol(bid)
            location = branch_cfg["location"]
            segment = neuron.segments[bid]
            connection_threshold = get(branch_cfg, "connection_threshold", 0.0)
            generator = if "gabor" ∈ keys(branch_cfg)
                # create Gabor patch
                GaborPatch(
                    branch_cfg["gabor"]["scale"],
                    [sin(branch_cfg["gabor"]["orientation"]*2π), cos(branch_cfg["gabor"]["orientation"]*2π)].*branch_cfg["gabor"]["frequency"],
                    branch_cfg["gabor"]["phase"],
                    branch_cfg["gabor"]["aperture"]
                )
            elseif "gaussian" ∈ keys(branch_cfg)
                pc1 = [sin(branch_cfg["gaussian"]["orientation"]*2π), cos(branch_cfg["gaussian"]["orientation"]*2π)]
                pc2 = [pc1[2],-pc1[1]]
                cov = [pc1 pc2]' * (branch_cfg["gaussian"]["stds"].^2 .* [pc1 pc2])
                GaussianPatch(
                    branch_cfg["gaussian"]["scale"],
                    inv(cov)
                )
            else
                continue
            end

            configure_rf!(segment, location, generator; connection_threshold=connection_threshold)
        end
    end
    return net
end


function summarize_spikes(times, dims, events, pixel_positions; start_event=:spike_start, end_event=:spike_end)
    frames = zeros(Int, length(times), dims...)
    state = zeros(Int, dims...)
    j=1
    for (i,t) ∈ enumerate(times)
        while j <= length(events) && events[j].t <= t
            ev = events[j]
            if isa(ev, Event{T, start_event, V} where {T,V}) && ev.value ∈ keys(pixel_positions)
                state[pixel_positions[ev.value]...] += 1
            elseif isa(ev, Event{T,end_event, V} where {T,V}) && ev.value ∈ keys(pixel_positions)
                state[pixel_positions[ev.value]...] -= 1
            end
            j += 1
        end
        frames[i,:,:] = state
    end
    return frames
end

function render_video(filename, frames, fps=30; crange=(minimum(frames),maximum(frames)), colorizer=val->(c=(val-crange[1])/(crange[2]-crange[1]);ColorTypes.RGB{ColorTypes.N0f8}(c,c,c)))
    props = [:priv_data => ("crf"=>"22","preset"=>"medium")]
    encodevideo(filename,[colorizer.(frames[i,:,:]) for i ∈ 1:size(frames,1)],framerate=fps,AVCodecContextProperties=props)
end

function plot_segment_weights(plt, segment::Segment, source_locations; cmap=:viridis, randomize=x->(rand(2).-0.5).*0.2.+x, markerscale=2)
    points = Tuple{Float64,Float64,Any,Float64}[]
    for (source,synapses) ∈ pairs(segment.synapses)
        if source ∉ keys(source_locations)
            continue
        end
        location = source_locations[source]
        p = 1.0 - prod(synapse->1.0-synapse.p, synapses)
        
        push!(points, (randomize(location)..., p, length(synapses)))
    end

    if ~isempty(points)
        (x,y,c,s)=collect.(zip(points...))
        scatter!(plt, x, y, markercolor=cmap, marker_z=c, markersize=s.*markerscale)
    end
    return plt
end
