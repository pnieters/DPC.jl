using ADSP, Distributions, DataFrames, VideoIO, ImageFiltering, ProgressMeter, Plots
import LinearAlgebra, ColorTypes


function load_video(io)
    vid = VideoIO.openvideo(io)
    
    df = DataFrame(timestamp=Float64[], frame=Matrix{<:Gray}[])
    while !eof(vid)
        # get time of the frame
        # read frame and convert to grayscale
        frame = ColorTypes.Gray.(read(vid))
        t = gettime(vid)
        push!(df, (t, frame))
    end
    df.timestamp .+= mean(diff(df.timestamp))
  
    close(vid)
    return df
end



function frames_to_spikes(df, ganglion_rf)
    @assert size(df,1) > 0 "`df` must not be empty!"
    dims = size(df[1,:frame])

    frame_durations = diff([0.0; df[:timestamp]])
    
    # Utility function to draw a Poisson-spike-train within the current frame interval [a,b]
    # according to the neuron-specific firing rates
    draw_spikes(rate, (a,b)) = rand(Uniform(a,b), rand(Poisson(max(0,rate*firing_rate_scale * (b-a)))))

    center_spikes = [Float64[] for i ∈ 1:prod(dims)]
    surround_spikes = [Float64[] for i ∈ 1:prod(dims)]

    for (t1,t2,f) ∈ zip([0.0; df.timestamp[1:end-1]], df.timestamp, df.frame)
        # For each receptive field location...
        for (i,inp) ∈ enumerate(imfilter(f, ganglion_rf))
            new_center_spikes = draw_spikes(inp, (t1,t2))
            new_surround_spikes = draw_spikes(-inp, (t1,t2))
            # ... generate spikes for center- and surround-cells
            append!(center_spikes[i],new_center_spikes)
            append!(surround_spikes[i],new_surround_spikes)
        end
    end
    
    return (center_spikes=center_spikes, surround_spikes=surround_spikes)
end


function sta(vid, duration, rf_region::Tuple{UnitRange{Int},UnitRange{Int}}, spikes::Vector{Tuple{Tuple{Int,Int}, Float64}})
    (y_rng,x_rng) = rf_region
    
    # parse the video
    fps = vid.framerate
    dims = vid.height,vid.width

    num_frames = Int(ceil(fps*duration))
    mean_frames = [zeros(Float64,length(y_rng),length(x_rng)) for i ∈ 1:num_frames]
    
    spike_counts = 0
    seekstart(vid)
    read(vid)
    for ((y,x),t) ∈ spikes
        if t < duration + 0.1
            continue
        end

        frames = []
        seek(vid, t-duration-0.1)
        t0=gettime(vid)
        while gettime(vid) <= t
            print("$(t) ")
            push!(frames, Gray.(read(vid))[y .+ y_rng, x .+ x_rng])
        end
        println("From $(t0) (desired: $(t-duration-0.1)) to $(gettime(vid)) (desired: $(t))")
        
        if length(frames) < num_frames
            println("Not enough frames ($(length(frames)) instead of $(num_frames))")
            continue
        else
            for (i,frame) ∈ enumerate(frames[end-num_frames+1:end])
                mean_frames[i] .+= frame
            end
            spike_counts += 1
        end
    end

    for i ∈ 1:num_frames
        mean_frames[i] ./= spike_counts
    end

    return mean_frames
end

struct Gabor
    grating::Vector{Float64}
    phase::Float64
    aperture::Float64
end

function (g::Gabor)(relative_position)
    scale = exp(-sum(relative_position.^2)/(2*g.aperture^2))
    offset = g.grating' * relative_position + g.phase
    return scale*cos(π*offset)
end

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
bipolar_width, horizontal_width = 2,4
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

df = load_video(io);
center_spikes, surround_spikes = frames_to_spikes(df, ganglion_rf)
#=
vid = VideoIO.openvideo(io)
center_spikes, surround_spikes, processed_frames, dims = video_to_spikes(vid, ganglion_rf)

# # write processed video
# props = [:priv_data => ("crf"=>"22","preset"=>"medium")]
# encodevideo("video.mp4",processed_frames,framerate=24,AVCodecContextProperties=props)

println("Collecting spikes...")
# # convert spikes to events for simulation
input_spike_events = generateEventsFromSpikes([center_spikes,surround_spikes], [:on, :off])

println("Setting up neurons...")
# assign a location to each RGC receptive field
retinotopic_positions = Dict(Symbol("$(onoff)$(id)") => collect(Tuple(c) .-(dims .+1)./2) for onoff ∈ (:on,:off) for (id,c) ∈ enumerate(CartesianIndices((1:dims[1],1:dims[2]))))

# load network template
net,metainfo = load_yaml(Network{Float64}, "examples/receptive_fields/cfg/template_simple.yaml")
configuration = metainfo[:configuration]

# fill network configuration with life given the metainformation
configure_rfs!(net, configuration)
# store instantiated network
save_yaml(net, "examples/receptive_fields/cfg/full_models/simple.yaml")

println("Running simulation...")
log_data=simulate!(net, input_spike_events; filter_events=ev->isa(ev, Event{T, :spike_start, V} where {T,V}) && ev.value==:simple, show_progress=true)

println("Estimating receptive fields...")
# close(vid)