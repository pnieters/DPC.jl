using ADSP, Distributions, VideoIO, ImageFiltering, ProgressMeter, Plots
import LinearAlgebra

function video_to_spikes(io, ganglion_rf)
    processed_frames = []
    
    # Utility function to draw a Poisson-spike-train within the current frame interval [a,b]
    # according to the neuron-specific firing rates
    draw_spikes(rate, (a,b)) = rand(Uniform(a,b), rand(Poisson(max(0,rate*firing_rate_scale * (b-a)))))
    
    meta_info = unsafe_load(io.apFormatContext[1])
    duration = VideoIO.get_duration(meta_info)
    
    # parse the video
    vid=VideoIO.openvideo(io)
    fps = vid.framerate
    dims = vid.height,vid.width
    center_spikes = [Float64[] for i ∈ 1:prod(dims)]
    surround_spikes = [Float64[] for i ∈ 1:prod(dims)]
    prog = Progress(Int(round(duration*fps)), 0.5, "Processing: ")
    let t=0.0, t′=0.0
        while !eof(vid)
            t = gettime(vid)
            next!(prog, showvalues=[("current time", t), ("video duration", duration)])
            # update time-interval of current frame
            t′=t+1.0/fps

            # read raw frame and convert it to grayscale
            frame = RGB.(read(vid))

            # filter the image to approximate retinal ganglion cell inputs
            ganglion_input = imfilter(Gray.(frame), ganglion_rf)

            # For each receptive field location...
            for (i,inp) ∈ enumerate(ganglion_input)
                # ... generate spikes for center- and surround-cells
                new_center_spikes = draw_spikes(inp, (t,t′))
                new_surround_spikes = draw_spikes(-inp, (t,t′))
                append!(center_spikes[i],new_center_spikes)
                append!(surround_spikes[i],new_surround_spikes)
                
                # ... draw spikes in the video frame
                α_g = 0.5^length(new_center_spikes)
                α_r = 0.5^length(new_surround_spikes)
                frame[i] = RGB(1.0-α_r + α_r*frame[i].r, 1.0-α_g+α_g*frame[i].g, frame[i].b)
            end

            # store the annotated video frame
            push!(processed_frames, RGB{ColorTypes.N0f8}.(frame))
        end
    end
    close(vid)
    
    return (center_spikes=center_spikes, surround_spikes=surround_spikes, processed_frames=processed_frames, dims=dims)
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

function configure_rfs!(net, configuration=net.metainfo[:configuration])
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

firing_rate_scale = 20.0

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
# io = VideoIO.testvideo("annie_oakley")
# center_spikes, surround_spikes, processed_frames, dims = video_to_spikes(io, ganglion_rf)

# # write processed video
# props = [:priv_data => ("crf"=>"22","preset"=>"medium")]
# encodevideo("video.mp4",processed_frames,framerate=24,AVCodecContextProperties=props)

# # convert spikes to events for simulation
# spike_events = generateEventsFromSpikes([center_spikes,surround_spikes], [:on, :off])
dims = 21,21

# assign a location to each RGC receptive field
retinotopic_positions = Dict(Symbol("$(onoff)$(id)") => collect(Tuple(c) .-(dims .+1)./2) for onoff ∈ (:on,:off) for (id,c) ∈ enumerate(CartesianIndices((1:dims[1],1:dims[2]))))

# load network template
net = load_yaml(Network{Float64, Any}, "examples/receptive_fields/cfg/template_simple.yaml")
# fill network configuration with life given the metainformation
configure_rfs!(net)
# store instantiated network
save_yaml(net, "examples/receptive_fields/cfg/full_models/simple.yaml")
