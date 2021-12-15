using DPC, DataStructures, Plots, DifferentialEquations, Distributions, DataFrames
# using DifferentialEquations, DataFrames, DataStructures, Plots, ProgressMeter
include("./integrator/utils.jl")

const mypalette = (orange="#fcaf3e", 
                   green="#8ae234", 
                   purple="#ad7fa8", 
                   turquise="#729f7e", 
                   yellow="#fce94f", 
                   blue="#729fcf", 
                   red="#ef2929", 
                   gray="#cccccc")



## Set up neuron
const num_inputs = 1000
const num_segments = 1000
const synapses_per_segment = 20
const transmission_probability = 0.5
const θ_seg = 300
const plateau_duration = 100e-3
const spike_duration = 5e-3
const refractory_duration = 5.01e-3
const trange = (0.0,3.0)
const λ = 200.0
const background_rate = 10.0
const input_names = Symbol.("in",1:num_inputs)
const seg_names = Symbol.("seg",1:num_segments)
const neuron_id = :integrator
const volleysize=10
const θ_syn = 5

net = Network(
    id_type=Symbol,
    time_type=Float64,
    weight_type=BernoulliSynapseWeight{Float64},
    synaptic_input_type=Int,
    default_spike_duration=spike_duration, 
    default_refractory_duration=refractory_duration, 
    default_plateau_duration=plateau_duration
)

objects = Dict{Symbol,Any}()
n = Neuron(neuron_id, net; 
            θ_syn=0, 
            θ_seg=θ_seg 
           )
objects[neuron_id] = n

for input_id in input_names
    inp = Input(input_id, net)
    objects[input_id] = inp
end

for (i, seg_id) in enumerate(seg_names)
    seg = Segment(seg_id, n; θ_syn=5, θ_seg=0)

    syn_inputs = input_names[mod.(i .+ (0:synapses_per_segment), Ref(1:num_inputs))]
    for (j, syn_inp) in enumerate(syn_inputs)
        syn_id = Symbol("syn$(i)$(j)_" * String(syn_inp))
        source = objects[syn_inp]
        target = seg

        syn = Synapse(syn_id, source, target; 
                      weight=BernoulliSynapseWeight(1.0,transmission_probability))
        objects[syn_id] = syn
    end
end

# save_network(net,"integrator/network.yaml")

# Set up stimulus

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


# ## Simulate!
spike_trains,v = generateContiguousVolleys(trange, num_inputs, stimulus, volleysize, background_rate; dt_spike = spike_duration)

all_spikes = Event{Float64}[]
for (i, pop_spikes) in enumerate(spike_trains)
    for spike_t in pop_spikes
        push!(all_spikes, Event(:input_spikes, 0.0, spike_t, objects[input_names[i]]) )
    end
end

# don't ask
v = [v]
spike_trains = [spike_trains]

# # simuluate the neuron
function fltr(t, tp, id, x)
    if tp == :plateau_starts || tp == :plateau_ends || tp == :spikes
        return true
    end
    return false
end
logger=simulate!(net, all_spikes; logger! =Logger(net,filter=fltr))
logger = logger.data

# # calculate plateau counts
# plateau_counts = sum.(eachrow(logger[r"seg*"]))

# # extract plateau strating times
plateau_starts = groupby(select!(filter(:event => (ev -> (ev == :plateau_starts)), logger),
                                :t, 
                                :object), 
                         :object)
plateau_starts_vec = [[ (seg_name,) in keys(plateau_starts.keymap) ? 
                        plateau_starts[(seg_name,)][!,:t] : Float64[] for seg_name in seg_names]]
# # extract the neuron's spike times
neuron_spikes = filter(:event => (ev -> (ev == :spikes)), logger)[:,:t] 

# ## Plot!

# Plot the inputs, active segments, the ideal "rectangular integrator", spike threshold in regions of interest
p_input = plot( [Shape([xoff, xoff+0.125, xoff+0.125, xoff],[0,0,750,750]) for xoff ∈ [0.25, 0.75, 1.25]], 
                 linecolor=nothing, 
                 fillcolor="#EEEEEE", 
                 title="dendritic integration of evidence", 
                 xlabel="time [s]", 
                 xlims=trange, 
                 ylims=(0,750), 
                 legend=:topright, 
                 label="", 
                 xticks=[0,0.25,0.5,0.75,1.0,1.25,1.5,2,3])

plot!(stimulus, trange..., linewidth=2, color=mypalette.orange, label="volleys/s")

hline!([θ_seg], linewidth=2, color=mypalette.red, label="spike threshold")

# plot!([trange[1];logger.t;trange[2]], 
#       [0;plateau_counts;0], 
#       linewidth=2, 
#       color=mypalette.blue, 
#       seriestype=:steppost, 
#       xlims=trange, 
#       label="active segments")

plot!(t->(int_stimulus(t)-int_stimulus(t-100e-3))/100e-3, 
      trange..., 
      linewidth=2, 
      color=:black, 
      linestyle=:dash, 
      label="rect. integrator")

# Plot two panels of spike-trains
v_cropped = [[(time=vv.time,idxs=filter(x->x∈1:50, vv.idxs)) for vv ∈ v[1]]]
spike_trains_cropped = [spike_trains[1][1:50]]
p_spikes1 = plot_spike_raster(trange, 
                              spike_trains_cropped, 
                              spike_duration; 
                              colors=[:gray], 
                              highlight=v_cropped, 
                              highligh_color=mypalette.orange, 
                              title="input spike-trains", 
                              xlims=trange, 
                              xticks=false)
yticks!(collect(-9.5:-20.0:-49.5), string.(collect(10:20:50)))
plot!(xaxis=:off)

v_cropped = [[(time=vv.time,idxs=filter(x->x∈951:1000, vv.idxs).-950) for vv ∈ v[1]]]
spike_trains_cropped = [spike_trains[1][951:1000]]
p_spikes2 = plot_spike_raster(trange, 
                              spike_trains_cropped, 
                              spike_duration; 
                              colors=[:gray], 
                              highlight=v_cropped, 
                              highligh_color=mypalette.orange, 
                              xlims=trange, 
                              xticks=false)
yticks!(collect(-9.5:-20.0:-49.5), string.(collect(960:20:1000)))

# Plot the "zoomed-in" regions of interest
subs = Plots.Plot[]
for xoff ∈ [0.25, 0.75, 1.25]
    sub_range = (xoff,xoff+0.125)
    pp = plot(#stimulus, sub_range..., color=mypalette.orange, label="volley rate", 
              background_color_inside="#EEEEEE",
              framestyle=:box,
              ylims=(0,750), 
              xlims=sub_range, 
              legend=false, 
              xticks=[xoff, xoff+0.0625, xoff+0.125])
    hline!([θ_seg],linewidth=3,color=mypalette.red)
    vline!([neuron_spikes],linewidth=1,color=mypalette.purple)

    # plot!([sub_range[1];logger.t;sub_range[2]], 
    #       [0;plateau_counts;0], 
    #       color=mypalette.blue, 
    #       linewidth=3, 
    #       seriestype=:steppost, 
    #       label="num. of plateaus")
    plot!(t->(int_stimulus(t)-int_stimulus(t-100e-3))/100e-3, 
          sub_range..., 
          color=:black, 
          linewidth=3, 
          linestyle=:dash, 
          label="rect. integrator")
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