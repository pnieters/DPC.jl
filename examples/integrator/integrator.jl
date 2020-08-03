using ADSP, DifferentialEquations, DataFrames, DataStructures, Plots, ProgressMeter
include("utils.jl")




## Set up neuron
const num_inputs = 1000
const num_segments = 1000
const synapses_per_segment = 20
const transmission_probability = 0.5
const θ_seg = 300
const plateau_duration = 100e-3
const spike_duration = 5e-3
const trange = (0.0,3.0)
const λ = 200.0
const background_rate = 10.0
const neuron_id = NeuronID(:integrator)
const input_names = Symbol.("in",1:num_inputs)
const seg_names = Symbol.("seg",1:num_segments)
const neuron_name = :integrator
const logged_segments = OrderedDict(segment=>PropertyID(SegmentID(neuron_id, segment), :active) for segment ∈ seg_names) 
const volleysize=10
const θ_syn = 5


segments = Dict{Symbol,Segment{Float64}}()
for (i,seg_name) ∈ enumerate(seg_names)
    seg_id = SegmentID(neuron_id, seg_name)
    # random_inputs = sample(input_names, synapses_per_segment, replace=false)
    # random_inputs = input_names[mod.(rand(1:num_inputs).+ (0:synapses_per_segment-1), Ref(1:num_inputs))]
    random_inputs = input_names[mod.(i.+ (0:synapses_per_segment-1), Ref(1:num_inputs))]
    synapses = Dict(source=>[Synapse(SynapseID(seg_id,source,i), transmission_probability, false)] for (i,source) ∈ enumerate(random_inputs))
    
    segments[seg_name] = Segment(seg_id, θ_syn, 0, false, plateau_duration, synapses, :soma, Symbol[])
end
segments[:soma] = Segment(SegmentID(neuron_id, :soma), 0, θ_seg, false, spike_duration, Dict{Symbol,Vector{Synapse}}(), nothing, seg_names)
n = Neuron(neuron_id, spike_duration, segments)
net = Network(Dict(neuron_name => n))
save_yaml(net, "examples/integrator/cfg/integrator.yaml")

## Set up stimulus

"""
Generate contiguous volleys at a high rate in 3 distinct periods.
"""
stimulus(t) = λ * if 0.25 ≤ t < 0.5
    1.0
elseif 0.75 ≤ t < 1.0
    2.0
elseif 1.25 ≤ t < 1.5
    3.0
elseif 0 ≤ t < 3.0
    0.1
else
    0.0
end

"""
Integral of the stimulus function (for plotting purposes).
"""
int_stimulus(t) = λ * if t < 0
    0
elseif 0.0 ≤ t < 0.25
    0.1(t-0.0)
elseif 0.25 ≤ t < 0.5
    1.0(t-0.25)+0.025
elseif 0.5 ≤ t < 0.75
    0.1(t-0.5)+0.25+0.025
elseif 0.75 ≤ t < 1.00
    2.0(t-0.75)+0.025+0.25+0.025
elseif 1.00 ≤ t < 1.25
    0.1(t-1.00)+0.5+0.025+0.25+0.025
elseif 1.25 ≤ t < 1.5
    3.0(t-1.25)+0.025+0.5+0.025+0.25+0.025
elseif 1.5 ≤ t < 3.5
    0.1(t-1.5)+0.75+0.025+0.5+0.025+0.25+0.025
else
    0.05+1.5+0.05+1.0+0.05+0.5+0.05
end


## Simulate!
# Generate the input spike-trains
spike_trains,v = generateContiguousVolleys(trange, num_inputs, stimulus, volleysize, background_rate; dt_spike = spike_duration)
spike_trains,v = [spike_trains], [v]
spike_events = sort!(generateEventsFromSpikes(spike_trains, ["in"]; spike_duration=spike_duration))

# simuluate the neuron
logger=simulate!(net, spike_events; show_progress=true, log=logged_segments, filter_events=ev->ev isa Union{Event{T,:plateau_start},Event{T,:plateau_end}} where {T})

# calculate plateau counts
plateau_counts = sum.(eachrow(logger[r"seg*"]))

# extract plateau strating times
plateau_starts = groupby(select!(filter(:event => (ev->isa(ev,Event{T,:plateau_start} where T)), logger), :t, :event => (ev -> map(x->x.value.id, ev)) => :segment), :segment)
plateau_starts_vec = [[ (seg_name,) ∈ keys(plateau_starts.keymap) ? plateau_starts[(seg_name,)][!,:t] : Float64[] for seg_name ∈ seg_names]]

# extract the neuron's spike times
neuron_spikes = plateau_starts[(:soma,)][:,:t]

## Plot!

# Plot the inputs, active segments, the ideal "rectangular integrator", spike threshold in regions of interest
p_input = plot([Shape([xoff, xoff+0.125, xoff+0.125, xoff],[0,0,750,750]) for xoff ∈ [0.25, 0.75, 1.25]], linecolor=nothing, fillcolor="#EEEEEE", title="dendritic integration of evidence", xlabel="time [s]", xlims=trange, ylims=(0,750), legend=:topright, label="", xticks=[0,0.25,0.5,0.75,1.0,1.25,1.5,2,3])
plot!(stimulus, trange..., linewidth=2, color=mypalette.orange, label="volleys/s")
hline!([θ_seg], linewidth=2, color=mypalette.red, label="spike threshold")
plot!([trange[1];logger.t;trange[2]], [0;plateau_counts;0], linewidth=2, color=mypalette.blue, seriestype=:steppost, xlims=trange, label="active segments")
plot!(t->(int_stimulus(t)-int_stimulus(t-100e-3))/100e-3, trange..., linewidth=2, color=:black, linestyle=:dash, label="rect. integrator")

# Plot two panels of spike-trains
v_cropped = [[(time=vv.time,idxs=filter(x->x∈1:50, vv.idxs)) for vv ∈ v[1]]]
spike_trains_cropped = [spike_trains[1][1:50]]
p_spikes1=plot_spike_raster(trange, spike_trains_cropped, spike_duration; colors=[:gray], highlight=v_cropped, highligh_color=mypalette.orange, title="input spike-trains", xlims=trange, xticks=false)
yticks!(collect(-9.5:-20.0:-49.5), string.(collect(10:20:50)))
plot!(xaxis=:off)
v_cropped = [[(time=vv.time,idxs=filter(x->x∈951:1000, vv.idxs).-950) for vv ∈ v[1]]]
spike_trains_cropped = [spike_trains[1][951:1000]]
p_spikes2=plot_spike_raster(trange, spike_trains_cropped, spike_duration; colors=[:gray], highlight=v_cropped, highligh_color=mypalette.orange, xlims=trange, xticks=false)
yticks!(collect(-9.5:-20.0:-49.5), string.(collect(960:20:1000)))

# Plot the "zoomed-in" regions of interest
subs = Plots.Plot[]
for xoff ∈ [0.25, 0.75, 1.25]
    sub_range = (xoff,xoff+0.125)
    pp=plot(#stimulus, sub_range..., color=mypalette.orange, label="volley rate", 
        background_color_inside="#EEEEEE",
        framestyle=:box,
        ylims=(0,750), xlims=sub_range, legend=false, xticks=[xoff, xoff+0.0625, xoff+0.125])
    hline!([θ_seg],linewidth=3,color=mypalette.red)
    vline!([neuron_spikes],linewidth=1,color=mypalette.purple)
    plot!([sub_range[1];logger.t;sub_range[2]], [0;plateau_counts;0], color=mypalette.blue, linewidth=3, seriestype=:steppost, label="num. of plateaus")
    plot!(t->(int_stimulus(t)-int_stimulus(t-100e-3))/100e-3, sub_range..., color=:black, linewidth=3, linestyle=:dash, label="rect. integrator")
    push!(subs,pp)
end

# Compose and save the figure
l = @layout [
    a1{0.15h}
    a2{0.15h}
    b{0.3h}
    [c{0.33w} d{0.33w} e{0.33w}]
]
plt=plot(p_spikes1, p_spikes2, p_input, subs..., layout=l, size=(500,600))
savefig("examples/integrator/figures/integrator.svg")
