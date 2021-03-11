#ENV["JULIA_DEBUG"] = "ADSP"
using ADSP, CairoMakie
include("utils.jl")

fig = Figure(resolution = (1200, 900))
ax1 = fig[1, 1] = Axis(fig, title = "Multi-neuron latch storage cells")
ax2 = fig[2, 1] = Axis(fig, title = "Readout cells")

config = """
refractory_duration: 5.01
inputs:
- id: i_set
- id: i_reset
- id: i_read
outputs:
neurons:
- id: n1
  θ_syn: 2
- id: n2
  θ_syn: 2
- id: n3
  θ_syn: 2
- id: n4
  θ_syn: 2
- id: n5
  θ_syn: 2
- id: n_readout_1
  θ_syn: 2
  branches:
    - id: seg_readout_1
      θ_syn: 2
- id: n_readout_2
  θ_syn: 2
  branches:
  - id: seg_readout_2
    θ_syn: 2
- id: n_readout_3
  θ_syn: 2
  branches:
    - id: seg_readout_3
      θ_syn: 2
- id: n_readout_4
  θ_syn: 2
  branches:
    - id: seg_readout_4
      θ_syn: 2
- id: n_readout_5
  θ_syn: 2
  branches:
    - id: seg_readout_5
      θ_syn: 2
       
synapses:
- {id: syn_set_1, source: i_set, target: n1, weight: [5.0, 0.75]}
- {id: syn_set_2, source: i_set, target: n2, weight: [5.0, 0.75]}
- {id: syn_set_3, source: i_set, target: n3, weight: [5.0, 0.75]}
- {id: syn_set_4, source: i_set, target: n4, weight: [5.0, 0.75]}
- {id: syn_set_5, source: i_set, target: n5, weight: [5.0, 0.75]}
- {id: syn_reset_1, source: i_reset, target: n1, weight: -10.0, spike_duration: 10}
- {id: syn_reset_2, source: i_reset, target: n2, weight: -10.0, spike_duration: 10}
- {id: syn_reset_3, source: i_reset, target: n3, weight: -10.0, spike_duration: 10}
- {id: syn_reset_4, source: i_reset, target: n4, weight: -10.0, spike_duration: 10}
- {id: syn_reset_5, source: i_reset, target: n5, weight: -10.0, spike_duration: 10}
- {id: syn_read_1, source: i_read, target: n_readout_1, weight: [5.0, 0.75]}
- {id: syn_read_2, source: i_read, target: n_readout_2, weight: [5.0, 0.75]}
- {id: syn_read_3, source: i_read, target: n_readout_3, weight: [5.0, 0.75]}
- {id: syn_read_4, source: i_read, target: n_readout_4, weight: [5.0, 0.75]}
- {id: syn_read_5, source: i_read, target: n_readout_5, weight: [5.0, 0.75]}
- {id: syn_recurrent_1_2, source: n1, target: n2, weight: [1.0, 0.75]}
- {id: syn_recurrent_1_3, source: n1, target: n3, weight: [1.0, 0.75]}
- {id: syn_recurrent_1_4, source: n1, target: n4, weight: [1.0, 0.75]}
- {id: syn_recurrent_1_5, source: n1, target: n5, weight: [1.0, 0.75]}
- {id: syn_recurrent_2_1, source: n2, target: n1, weight: [1.0, 0.75]}
- {id: syn_recurrent_2_3, source: n2, target: n3, weight: [1.0, 0.75]}
- {id: syn_recurrent_2_4, source: n2, target: n4, weight: [1.0, 0.75]}
- {id: syn_recurrent_2_5, source: n2, target: n5, weight: [1.0, 0.75]}
- {id: syn_recurrent_3_1, source: n3, target: n1, weight: [1.0, 0.75]}
- {id: syn_recurrent_3_2, source: n3, target: n2, weight: [1.0, 0.75]}
- {id: syn_recurrent_3_4, source: n3, target: n4, weight: [1.0, 0.75]}
- {id: syn_recurrent_3_5, source: n3, target: n5, weight: [1.0, 0.75]}
- {id: syn_recurrent_4_1, source: n4, target: n1, weight: [1.0, 0.75]}
- {id: syn_recurrent_4_2, source: n4, target: n2, weight: [1.0, 0.75]}
- {id: syn_recurrent_4_3, source: n4, target: n3, weight: [1.0, 0.75]}
- {id: syn_recurrent_4_5, source: n4, target: n5, weight: [1.0, 0.75]}
- {id: syn_recurrent_5_1, source: n5, target: n1, weight: [1.0, 0.75]}
- {id: syn_recurrent_5_2, source: n5, target: n2, weight: [1.0, 0.75]}
- {id: syn_recurrent_5_3, source: n5, target: n3, weight: [1.0, 0.75]}
- {id: syn_recurrent_5_4, source: n5, target: n4, weight: [1.0, 0.75]}
- {id: syn_forward_1_2, source: n1, target: seg_readout_2, weight: [1.0, 0.75]}
- {id: syn_forward_1_3, source: n1, target: seg_readout_3, weight: [1.0, 0.75]}
- {id: syn_forward_1_4, source: n1, target: seg_readout_4, weight: [1.0, 0.75]}
- {id: syn_forward_1_5, source: n1, target: seg_readout_5, weight: [1.0, 0.75]}
- {id: syn_forward_2_1, source: n2, target: seg_readout_1, weight: [1.0, 0.75]}
- {id: syn_forward_2_3, source: n2, target: seg_readout_3, weight: [1.0, 0.75]}
- {id: syn_forward_2_4, source: n2, target: seg_readout_4, weight: [1.0, 0.75]}
- {id: syn_forward_2_5, source: n2, target: seg_readout_5, weight: [1.0, 0.75]}
- {id: syn_forward_3_1, source: n3, target: seg_readout_1, weight: [1.0, 0.75]}
- {id: syn_forward_3_2, source: n3, target: seg_readout_2, weight: [1.0, 0.75]}
- {id: syn_forward_3_4, source: n3, target: seg_readout_4, weight: [1.0, 0.75]}
- {id: syn_forward_3_5, source: n3, target: seg_readout_5, weight: [1.0, 0.75]}
- {id: syn_forward_4_1, source: n4, target: seg_readout_1, weight: [1.0, 0.75]}
- {id: syn_forward_4_2, source: n4, target: seg_readout_2, weight: [1.0, 0.75]}
- {id: syn_forward_4_3, source: n4, target: seg_readout_3, weight: [1.0, 0.75]}
- {id: syn_forward_4_5, source: n4, target: seg_readout_5, weight: [1.0, 0.75]}
- {id: syn_forward_5_1, source: n5, target: seg_readout_1, weight: [1.0, 0.75]}
- {id: syn_forward_5_2, source: n5, target: seg_readout_2, weight: [1.0, 0.75]}
- {id: syn_forward_5_3, source: n5, target: seg_readout_3, weight: [1.0, 0.75]}
- {id: syn_forward_5_4, source: n5, target: seg_readout_4, weight: [1.0, 0.75]}
"""

(net,objects) = load_network(YAML_source=config, weight_type=BernoulliSynapseWeight{Float64})

input=sort!([
    [Event(:input_spikes, 0.0,  t, objects[:i_set]) for t ∈ [105.0, 605.0]];
    [Event(:input_spikes, 0.0,  t, objects[:i_reset]) for t ∈ [400.0, 800.0]];
    [Event(:input_spikes, 0.0,  t, objects[:i_read]) for t ∈ [150.0, 300.0, 450.0, 600.0, 750.0, 900.0, 1050.0, 1200.0, 1350.0, 1500.0]];
])

logger=simulate!(net, input, 1000)

n1_spikes = filter(x->(x.object==:n1 && x.event == :spikes), logger.data)
n2_spikes = filter(x->(x.object==:n2 && x.event == :spikes), logger.data)
n3_spikes = filter(x->(x.object==:n3 && x.event == :spikes), logger.data)
n4_spikes = filter(x->(x.object==:n4 && x.event == :spikes), logger.data)
n5_spikes = filter(x->(x.object==:n5 && x.event == :spikes), logger.data)

syn_set_1 = get_trace(:syn_set_1, logger.data)
syn_set_2 = get_trace(:syn_set_2, logger.data)
syn_set_3 = get_trace(:syn_set_3, logger.data)
syn_set_4 = get_trace(:syn_set_4, logger.data)
syn_set_5 = get_trace(:syn_set_5, logger.data)

syn_reset_1 = get_trace(:syn_reset_1, logger.data)
syn_reset_2 = get_trace(:syn_reset_2, logger.data)
syn_reset_3 = get_trace(:syn_reset_3, logger.data)
syn_reset_4 = get_trace(:syn_reset_4, logger.data)
syn_reset_5 = get_trace(:syn_reset_5, logger.data)

syn_read_1 = get_trace(:syn_read_1, logger.data)
syn_read_2 = get_trace(:syn_read_2, logger.data)
syn_read_3 = get_trace(:syn_read_3, logger.data)
syn_read_4 = get_trace(:syn_read_4, logger.data)
syn_read_5 = get_trace(:syn_read_5, logger.data)

seg_read_1 = get_trace(:seg_readout_1, logger.data)
seg_read_2 = get_trace(:seg_readout_2, logger.data)
seg_read_3 = get_trace(:seg_readout_3, logger.data)
seg_read_4 = get_trace(:seg_readout_4, logger.data)
seg_read_5 = get_trace(:seg_readout_5, logger.data)

n_readout_1_spikes = filter(x->(x.object==:n_readout_1 && x.event == :spikes), logger.data)
n_readout_2_spikes = filter(x->(x.object==:n_readout_2 && x.event == :spikes), logger.data)
n_readout_3_spikes = filter(x->(x.object==:n_readout_3 && x.event == :spikes), logger.data)
n_readout_4_spikes = filter(x->(x.object==:n_readout_4 && x.event == :spikes), logger.data)
n_readout_5_spikes = filter(x->(x.object==:n_readout_5 && x.event == :spikes), logger.data)

steps!(ax1, [0;syn_set_1.t;1000],    0 .+ 0.9 .* [0;Int.(syn_set_1.state);0],   fill=color_1_50)
steps!(ax1, [0;syn_reset_1.t;1000],  0 .+ 0.9 .* [0;Int.(syn_reset_1.state);0], fill=color_2_50)
linesegments!(ax1, repeat(n1_spikes.t,inner=2), repeat(     [0,0.9], outer=length(n1_spikes.t)))

steps!(ax1, [0;syn_set_2.t;1000],    1 .+ 0.9 .* [0;Int.(syn_set_2.state);0],   fill=color_1_50)
steps!(ax1, [0;syn_reset_2.t;1000],  1 .+ 0.9 .* [0;Int.(syn_reset_2.state);0], fill=color_2_50)
linesegments!(ax1, repeat(n2_spikes.t,inner=2), repeat(1 .+ [0,0.9], outer=length(n2_spikes.t)))

steps!(ax1, [0;syn_set_3.t;1000],    2 .+ 0.9 .* [0;Int.(syn_set_3.state);0],   fill=color_1_50)
steps!(ax1, [0;syn_reset_3.t;1000],  2 .+ 0.9 .* [0;Int.(syn_reset_3.state);0], fill=color_2_50)
linesegments!(ax1, repeat(n3_spikes.t,inner=2), repeat(2 .+ [0,0.9], outer=length(n3_spikes.t)))

steps!(ax1, [0;syn_set_4.t;1000],    3 .+ 0.9 .* [0;Int.(syn_set_4.state);0],   fill=color_1_50)
steps!(ax1, [0;syn_reset_4.t;1000],  3 .+ 0.9 .* [0;Int.(syn_reset_4.state);0], fill=color_2_50)
linesegments!(ax1, repeat(n4_spikes.t,inner=2), repeat(3 .+ [0,0.9], outer=length(n4_spikes.t)))

steps!(ax1, [0;syn_set_5.t;1000],    4 .+ 0.9 .* [0;Int.(syn_set_5.state);0],   fill=color_1_50)
steps!(ax1, [0;syn_reset_5.t;1000],  4 .+ 0.9 .* [0;Int.(syn_reset_5.state);0], fill=color_2_50)
linesegments!(ax1, repeat(n5_spikes.t,inner=2), repeat(4 .+ [0,0.9], outer=length(n5_spikes.t)))


steps!(ax2, [0;seg_read_1.t;1000],    0 .+ 0.45 .* [0;Int.(seg_read_1.state);0],   fill=color_3_50)
steps!(ax2, [0;syn_read_1.t;1000],    0 .+ 0.9 .* [0;Int.(syn_read_1.state);0],   fill=color_1_50)
linesegments!(ax2, repeat(n_readout_1_spikes.t,inner=2), repeat(     [0,0.9], outer=length(n_readout_1_spikes.t)), linewidth=3, color=:yellow)

steps!(ax2, [0;seg_read_2.t;1000],    1 .+ 0.45 .* [0;Int.(seg_read_2.state);0],   fill=color_3_50)
steps!(ax2, [0;syn_read_2.t;1000],    1 .+ 0.9 .* [0;Int.(syn_read_2.state);0],   fill=color_1_50)
linesegments!(ax2, repeat(n_readout_2_spikes.t,inner=2), repeat(1 .+ [0,0.9], outer=length(n_readout_2_spikes.t)), linewidth=3, color=:yellow)

steps!(ax2, [0;seg_read_3.t;1000],    2 .+ 0.45 .* [0;Int.(seg_read_3.state);0],   fill=color_3_50)
steps!(ax2, [0;syn_read_3.t;1000],    2 .+ 0.9 .* [0;Int.(syn_read_3.state);0],   fill=color_1_50)
linesegments!(ax2, repeat(n_readout_3_spikes.t,inner=2), repeat(2 .+ [0,0.9], outer=length(n_readout_3_spikes.t)), linewidth=3, color=:yellow)

steps!(ax2, [0;seg_read_4.t;1000],    3 .+ 0.45 .* [0;Int.(seg_read_4.state);0],   fill=color_3_50)
steps!(ax2, [0;syn_read_4.t;1000],    3 .+ 0.9 .* [0;Int.(syn_read_4.state);0],   fill=color_1_50)
linesegments!(ax2, repeat(n_readout_4_spikes.t,inner=2), repeat(3 .+ [0,0.9], outer=length(n_readout_4_spikes.t)), linewidth=3, color=:yellow)

steps!(ax2, [0;seg_read_5.t;1000],    4 .+ 0.45 .* [0;Int.(seg_read_5.state);0],   fill=color_3_50)
steps!(ax2, [0;syn_read_5.t;1000],    4 .+ 0.9 .* [0;Int.(syn_read_5.state);0],   fill=color_1_50)
linesegments!(ax2, repeat(n_readout_5_spikes.t,inner=2), repeat(4 .+ [0,0.9], outer=length(n_readout_5_spikes.t)), linewidth=3, color=:yellow)

save("stochastic_latch.svg", fig)