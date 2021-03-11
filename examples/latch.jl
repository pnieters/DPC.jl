#ENV["JULIA_DEBUG"] = "ADSP"
using ADSP, CairoMakie
include("utils.jl")

fig = Figure(resolution = (1200, 900))
ax1 = fig[1, 1] = Axis(fig, title = "Single neuron latch")
ax2 = fig[2, 1] = Axis(fig, title = "Two neuron latch")

################################################################################
config = """
refractory_duration: 5.01
inputs:
- id: i_set
- id: i_reset
outputs:
neurons:
- id: n_storage
synapses:
- {id: syn1, source: i_set, target: n_storage}
- {id: syn2, source: i_reset, target: n_storage, weight: -10, spike_duration: 10}
- {id: syn3, source: n_storage, target: n_storage}
"""

(net,objects) = load_network(YAML_source=config)

input=sort!([
    [Event(:input_spikes, 0.0,  t, objects[:i_set]) for t ∈ [105.0, 605.0]];
    [Event(:input_spikes, 0.0,  t, objects[:i_reset]) for t ∈ [400.0, 800.0]];
])

logger=simulate!(net, input, 1000)
n    = get_trace(:n_storage, logger.data)
n_storage_spikes = filter(x->(x.object==:n_storage && x.event == :spikes), logger.data)
syn1   = get_trace(:syn1, logger.data)
syn2   = get_trace(:syn2, logger.data)

steps!(ax1, [0;syn1.t;1000], 1 .+ 0.45 .* [0;Int.(syn1.state);0] , fill=color_1_50)
steps!(ax1, [0;syn2.t;1000], 2 .+ 0.45 .* [0;Int.(syn2.state);0] , fill=color_2_50)
# steps!(ax1, [0;n.t;1000], 0.45 .* [0;Int.(n.state);0] , fill=color_3_50)
linesegments!(ax1, repeat(n_storage_spikes.t,inner=2), repeat([0,0.9], outer=length(n_storage_spikes.t)))

################################################################################

config = """
refractory_duration: 5.01
inputs:
- id: i_set
- id: i_reset
- id: i_read
outputs:
neurons:
- id: n_storage
- id: n_readout
  branches:
    - id: seg_readout
synapses:
- {id: syn1, source: i_set, target: n_storage}
- {id: syn2, source: i_reset, target: n_storage, weight: -10, spike_duration: 10}
- {id: syn3, source: n_storage, target: n_storage}
- {id: syn4, source: n_storage, target: seg_readout}
- {id: syn5, source: i_read, target: n_readout}
"""

(net,objects) = load_network(YAML_source=config)

input=sort!([
    [Event(:input_spikes, 0.0,  t, objects[:i_set]) for t ∈ [105.0, 1105.0]];
    [Event(:input_spikes, 0.0,  t, objects[:i_reset]) for t ∈ [400.0, 1400.0]];
    [Event(:input_spikes, 0.0,  t, objects[:i_read]) for t ∈ [150.0, 300.0, 450.0, 600.0, 750.0, 900.0, 1050.0, 1200.0, 1350.0, 1500.0]];
])

logger=simulate!(net, input)

n_storage   = get_trace(:n_storage, logger.data)
n_storage_spikes = filter(x->(x.object==:n_storage && x.event == :spikes), logger.data)
n_readout   = get_trace(:n_readout, logger.data)
syn1   = get_trace(:syn1, logger.data)
syn2   = get_trace(:syn2, logger.data)
syn5   = get_trace(:syn5, logger.data)


steps!(ax2, [0;syn1.t;1000], 3 .+ 0.45 .* [0;Int.(syn1.state);0] , fill=color_1_50)
steps!(ax2, [0;syn2.t;1000], 4 .+ 0.45 .* [0;Int.(syn2.state);0] , fill=color_2_50)
steps!(ax2, [0;syn5.t;1000], 1 .+ 0.45 .* [0;Int.(syn5.state);0] , fill=RGBAf0(1,0,1,0.5))
linesegments!(ax2, repeat(n_storage_spikes.t,inner=2), repeat([2,2.9], outer=length(n_storage_spikes.t)))
#steps!(ax2, [0;n_storage.t;1000], 2 .+ 0.45 .* [0;Int.(n_storage.state);0] )#, fill=color_3_50)
steps!(ax2, [0;n_readout.t;1000], 0.45 .* [0;Int.(n_readout.state);0] , fill=RGBAf0(0,1,1,0.5))

################################################################################
save("latch.svg", fig)
