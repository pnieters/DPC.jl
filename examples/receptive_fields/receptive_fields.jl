using ADSP, Distributions, DataFrames, VideoIO, ImageFiltering, ProgressMeter, Plots
import LinearAlgebra, ColorTypes

include("./utils.jl")


firing_rate_scale = 400.0

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


##  parse video input to spikes ##
println("Parsing video...")
# load video
#io = VideoIO.testvideo("annie_oakley")
# load video into buffer
#buff = load_video(io; range=(1:100,1:100));
times=collect(0.01:0.2:100)
buff = VideoBuffer(times, rand(length(times), 20, 20))
# filter each frame with the RGC difference-of-gaussian kernel
filtered_buff = imfilter(buff, ganglion_rf)
# supersample each RGC (each pixel is represented by a `supersample`x`supersample` population)
upscale = 3
rgc_activations = VideoBuffer(filtered_buff.timestamps, repeat(filtered_buff.frames,inner=(1,upscale,upscale)))
# sample RGC spikes
center_spikes = frames_to_spikes(rgc_activations)
surround_spikes = frames_to_spikes(rgc_activations; transform = x -> -x)

# assign a location to each RGC's receptive field
dims = size(rgc_activations.frames)[2:3]
id = 1
on_center_ids = Symbol[]
off_center_ids = Symbol[]
pixel_positions = Dict{Symbol,Tuple{Int,Int}}()
retinotopic_positions = Dict{Symbol,Tuple{Float64,Float64}}()
for x ∈ 1:dims[1]
    for y ∈ 1:dims[2]
        for onoff ∈ (:on,:off)
            key = Symbol("$(onoff)$(id)")
            pixel_positions[key] = (div(y-1,upscale)+1,div(x-1,upscale)+1)
            retinotopic_positions[key] = ((y, x) .- (dims .+ upscale)./2)./upscale
            
            if onoff == :on
                push!(on_center_ids, key)
            else
                push!(off_center_ids, key)
            end
        end
        global id += 1
    end
end

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
    pixel_positions[key] = Tuple(Int.(round.(value .+ (div.(dims,upscale) .+ 1)./2))).+1
    retinotopic_positions[key] = Tuple(value)
end


# fill network configuration with life given the metainformation
configure_rfs!(net, configuration, retinotopic_positions)
# store instantiated network
save_yaml(net, "examples/receptive_fields/cfg/full_models/simple.yaml")

# # plot simple cell's receptive field
# p=plot()
# plot_segment_weights(p, net.neurons[:simple].segments[:soma], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
# plot_segment_weights(p, net.neurons[:simple].segments[:soma], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
# display(p)
# plot compartmentalized cell's receptive field
p=plot()
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:left_1], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:left_1], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:left_2], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:left_2], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:right_2], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:right_2], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:right_1], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:right_1], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:center_1], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:center_1], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:center_2], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:center_2], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:center_3], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
plot_segment_weights(p, net.neurons[:compartmentalized].segments[:center_3], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
display(p)
# plot_segment_weights(p, net.neurons[:compartmentalized].segments[:soma], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
# plot_segment_weights(p, net.neurons[:compartmentalized].segments[:soma], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
# plot_segment_weights(p, net.neurons[:compartmentalized].segments[:center_branch], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
# plot_segment_weights(p, net.neurons[:compartmentalized].segments[:center_branch], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
# plot_segment_weights(p, net.neurons[:compartmentalized].segments[:left_branch], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
# plot_segment_weights(p, net.neurons[:compartmentalized].segments[:left_branch], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)
# plot_segment_weights(p, net.neurons[:compartmentalized].segments[:right_branch], filter(((k,v),)->k∈on_center_ids, retinotopic_positions); cmap=:reds)
# plot_segment_weights(p, net.neurons[:compartmentalized].segments[:right_branch], filter(((k,v),)->k∈off_center_ids, retinotopic_positions); cmap=:blues)


println("Running simulation...")
log_data=simulate!(net, input_spike_events; show_progress=true) # filter_events=ev->isa(ev, Event{T, :spike_start, V} where {T,V}) && ev.value∈(:simple,:compartmentalized)

println("Estimating receptive fields...")

# collect all spikes and locations of RGCs
on_center_spikes_locations = [(pixel_positions[ev.value],ev.t) for ev ∈ filter(ev->isa(ev, Event{T, :spike_start, V} where {T,V}) && ev.value ∈ on_center_ids, input_spike_events)]
off_center_spikes_locations = [(pixel_positions[ev.value],ev.t) for ev ∈ filter(ev->isa(ev, Event{T, :spike_start, V} where {T,V}) && ev.value ∈ off_center_ids, input_spike_events)]
# collect all spikes and locations of simple cells
simple_cell_spikes_locations = [(pixel_positions[ev.value],ev.t) for ev ∈ filter(ev->isa(ev, Event{T, :spike_start, V} where {T,V}) && ev.value==:simple, log_data.event)]
compartmentalized_cell_spikes_locations = [(pixel_positions[ev.value],ev.t) for ev ∈ filter(ev->isa(ev, Event{T, :spike_start, V} where {T,V}) && ev.value==:compartmentalized, log_data.event)]

sta_range = (-4:0, -5:5,-5:5)
# Calculate STA of all on-center RGCs
on_center_sta = calculate_sta(buff, on_center_spikes_locations, sta_range)
off_center_sta = calculate_sta(buff, off_center_spikes_locations, sta_range)
# Calculate STA of simple cells
simple_cell_sta = calculate_sta(buff, simple_cell_spikes_locations, sta_range)
compartmentalized_cell_sta = calculate_sta(buff, compartmentalized_cell_spikes_locations, sta_range)
#=
# left branch
left_1_spikes_locations = [([11,11],ev.t) for ev ∈ filter(ev->isa(ev, Event{T,:plateau_start,Q} where {T,Q}) && ev.value.id==:left_1, log_data.event)]
left_2_spikes_locations = [([11,11],ev.t) for ev ∈ filter(ev->isa(ev, Event{T,:plateau_start,Q} where {T,Q}) && ev.value.id==:left_2, log_data.event)]
right_1_spikes_locations = [([11,11],ev.t) for ev ∈ filter(ev->isa(ev, Event{T,:plateau_start,Q} where {T,Q}) && ev.value.id==:right_1, log_data.event)]
right_2_spikes_locations = [([11,11],ev.t) for ev ∈ filter(ev->isa(ev, Event{T,:plateau_start,Q} where {T,Q}) && ev.value.id==:right_2, log_data.event)]
center_1_spikes_locations = [([11,11],ev.t) for ev ∈ filter(ev->isa(ev, Event{T,:plateau_start,Q} where {T,Q}) && ev.value.id==:center_1, log_data.event)]
center_2_spikes_locations = [([11,11],ev.t) for ev ∈ filter(ev->isa(ev, Event{T,:plateau_start,Q} where {T,Q}) && ev.value.id==:center_2, log_data.event)]
center_3_spikes_locations = [([11,11],ev.t) for ev ∈ filter(ev->isa(ev, Event{T,:plateau_start,Q} where {T,Q}) && ev.value.id==:center_3, log_data.event)]

left_1_sta = calculate_sta(buff, left_1_spikes_locations, sta_range)
left_2_sta = calculate_sta(buff, left_2_spikes_locations, sta_range)
right_1_sta = calculate_sta(buff, right_1_spikes_locations, sta_range)
right_2_sta = calculate_sta(buff, right_2_spikes_locations, sta_range)
center_1_sta = calculate_sta(buff, center_1_spikes_locations, sta_range)
center_2_sta = calculate_sta(buff, center_2_spikes_locations, sta_range)
center_3_sta = calculate_sta(buff, center_3_spikes_locations, sta_range)

# heatmap(left_1_sta[end,:,:])
# heatmap(left_2_sta[end,:,:])
# heatmap(right_1_sta[end,:,:])
# heatmap(right_2_sta[end,:,:])
# heatmap(center_1_sta[end,:,:])
# heatmap(center_2_sta[end,:,:])
# heatmap(center_3_sta[end,:,:])
=#
# duration = 25
# fps = 1000
heatmap(compartmentalized_cell_sta[end,:,:])

# frames = summarize_spikes(0:0.01:10.0, div.(dims, upscale), input_spike_events, pixel_positions)
# render_video("spikes.mp4", frames)


