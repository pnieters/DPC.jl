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
function load_video(io::VideoIO.AVInput; c=ColorTypes.Gray)
    vid = VideoIO.openvideo(io)

    timestamps = Float64[]
    frames = []

    while !eof(vid)
        # get time of the frame, read frame and convert to correct color type
        frame = c.(read(vid))
        t = gettime(vid)
        push!(frames, reshape(frame, (1,size(frame)...)))
        push!(timestamps, t)
    end

    buff = VideoBuffer(timestamps, cat(frames...; dims=1))
    buff.timestamps .+= mean(diff(buff.timestamps))
  
    close(vid)
    return buff
end


"""
    frames_to_spikes(buff, ganglion_rf)

Applies a receptive field to each individual frame of a VideoBuffer and converts the result into a Poisson spike-train.
"""
function frames_to_spikes(buff, ganglion_rf)
    @assert size(buff.frames,1) > 0 "`df` must not be empty!"
    dims = size(buff.frames)[2:3]

    frame_durations = diff([0.0; buff.timestamps])
    
    # Utility function to draw a Poisson-spike-train within the current frame interval [a,b]
    # according to the neuron-specific firing rates
    draw_spikes(rate, (a,b)) = rand(Uniform(a,b), rand(Poisson(max(0,rate*firing_rate_scale * (b-a)))))

    center_spikes = [Float64[] for i ∈ 1:prod(dims)]
    surround_spikes = [Float64[] for i ∈ 1:prod(dims)]

    for (i,(t1,t2)) ∈ enumerate(zip([0.0; buff.timestamps[1:end-1]], buff.timestamps))
        # For each receptive field location...
        for (j,inp) ∈ enumerate(imfilter(buff[i,:,:], ganglion_rf))
            new_center_spikes = draw_spikes(inp, (t1,t2))
            new_surround_spikes = draw_spikes(-inp, (t1,t2))
            # ... generate spikes for center- and surround-cells
            append!(center_spikes[j],new_center_spikes)
            append!(surround_spikes[j],new_surround_spikes)
        end
    end
    
    return (center_spikes=center_spikes, surround_spikes=surround_spikes)
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

struct Gabor
    grating::Vector{Float64}
    phase::Float64
    aperture::Float64
end

function (g::Gabor)(relative_position)
    scale = exp(-sum(relative_position.^2)/(2*g.aperture^2))
    offset = 1.0/g.aperture * g.grating' * relative_position + g.phase
    return scale*cos(π*offset)
end
(g::Gabor)(dy,dx) = g([dy,dx])

function configure_rfs!(net, configuration)
    for (nid,neuron_cfg) ∈ configuration
        nid=Symbol(nid)
        neuron = net.neurons[nid]
        for (bid,branch_cfg) ∈ neuron_cfg
            bid=Symbol(bid)
            segment = neuron.segments[bid]
            
            # create Gabor patch
            g = Gabor(
                [sin(branch_cfg["orientation"]*2π), cos(branch_cfg["orientation"]*2π)],
                branch_cfg["phase"],
                branch_cfg["aperture"]
            )
            connection_threshold = 1e-2
            
            # filter and weight relevant sources for synaptic connections to segment
            for (source_id,position) ∈ retinotopic_positions
                synaptic_strength = g(position - branch_cfg["location"])
                synapseID = SynapseID(SegmentID(NeuronID(nid),bid),source_id,1)
                if synaptic_strength < -connection_threshold && startswith(String(source_id), "on")
                    segment.synapses[source_id] = [Synapse(synapseID,-synaptic_strength,false)]
                elseif synaptic_strength > connection_threshold && startswith(String(source_id), "o")
                    segment.synapses[source_id] = [Synapse(synapseID, synaptic_strength,false)]
                end
            end
        end
    end
    return net
end

firing_rate_scale = 200.0

#= Set up receptive fields =#
bipolar_weight, horizontal_weight = 1,1
bipolar_width, horizontal_width = 1,2
# narrow ON-center receptive field
bipolar_rf = Kernel.gaussian((bipolar_width,bipolar_width),(horizontal_width,horizontal_width).*4 .+1)
# wider OFF-center receptive field
horizontal_rf = Kernel.gaussian((horizontal_width,horizontal_width),(horizontal_width,horizontal_width).*4 .+1)
# ON-OFF ganglion cell receptive field
ganglion_rf = bipolar_weight*bipolar_rf - horizontal_weight*horizontal_rf
ganglion_rf ./= LinearAlgebra.norm(ganglion_rf)


# # parse video input
println("Parsing video...")
io = VideoIO.testvideo("annie_oakley")

buff = load_video(io);
center_spikes, surround_spikes = frames_to_spikes(buff, ganglion_rf)
# assign a location to each RGC receptive field
dims = size(buff.frames)[2:3]
retinotopic_positions = Dict(Symbol("$(onoff)$(id)") => collect(Tuple(c) .-(dims .+1)./2) for onoff ∈ (:on,:off) for (id,c) ∈ enumerate(CartesianIndices((1:dims[1],1:dims[2]))))
pixel_positions = Dict(Symbol("$(onoff)$(id)") => Tuple(c) for onoff ∈ (:on,:off) for (id,c) ∈ enumerate(CartesianIndices((1:dims[1],1:dims[2]))))

println("Collecting spikes...")
# # convert spikes to events for simulation
input_spike_events = generateEventsFromSpikes([center_spikes,surround_spikes], [:on, :off])

println("Setting up neurons...")
# load network template
net,metainfo = load_yaml(Network{Float64}, "examples/receptive_fields/cfg/template_simple.yaml")
configuration = metainfo[:configuration]
for (key,value) ∈ pairs(configuration)
    key = Symbol(key)
    value = value["soma"]["location"]
    retinotopic_positions[key] = value
    pixel_positions[key] = Int.(ceil.(Tuple(value) .+ (dims .+1)./2))
end


# fill network configuration with life given the metainformation
configure_rfs!(net, configuration)
# store instantiated network
save_yaml(net, "examples/receptive_fields/cfg/full_models/simple.yaml")

println("Running simulation...")
log_data=simulate!(net, input_spike_events; filter_events=ev->isa(ev, Event{T, :spike_start, V} where {T,V}) && ev.value==:simple, show_progress=true)

println("Estimating receptive fields...")

# collect all spikes and locations of RGCs
on_center_spikes_locations = [(pixel_positions[ev.value],ev.t) for ev ∈ filter(ev->isa(ev, Event{T, :spike_start, V} where {T,V}) && startswith(String(ev.value),"on"), input_spike_events)]
off_center_spikes_locations = [(pixel_positions[ev.value],ev.t) for ev ∈ filter(ev->isa(ev, Event{T, :spike_start, V} where {T,V}) && startswith(String(ev.value),"off"), input_spike_events)]
# collect all spikes and locations of simple cells
simple_cell_spikes_locations = [(pixel_positions[ev.value],ev.t) for ev ∈ log_data.event]

sta_range = (-4:0, -5:5,-5:5)
# Calculate STA of all on-center RGCs
on_center_sta = calculate_sta(buff, on_center_spikes_locations, sta_range)
off_center_sta = calculate_sta(buff, off_center_spikes_locations, sta_range)
# Calculate STA of simple cells
simple_cell_sta = calculate_sta(buff, simple_cell_spikes_locations, sta_range)

