using ADSP, CairoMakie
#include(joinpath(@__DIR__, "utils.jl"))
include("utils.jl")


fig = Figure(resolution = (1200, 900))
ax1 = fig[1, 1] = Axis(fig, title = "Coincidence at soma")
ax2 = fig[2, 1] = Axis(fig, title = "Single segment")
ax3 = fig[3, 1] = Axis(fig, title = "Two segments")
ax4 = fig[4, 1] = Axis(fig, title = "Two segments, bursts")

################################################################################

config = """
refractory_duration: 5.01
inputs:
- id: i1
- id: i2
outputs:
- id: o
neurons:
- id: n
  θ_syn: 2
synapses:
- {id: syn1, source: i1, target: n}
- {id: syn2, source: i2, target: n}
- {id: syn_out, source: n, target: o, spike_duration: 0, delay: 0}
"""

(net,objects) = load_network(YAML_source=config)

input=[
    Event(:input_spikes, 0.0,  10.0, objects[:i1]),
    Event(:input_spikes, 0.0,   0.0, objects[:i2]),
    Event(:input_spikes, 0.0, 150.0, objects[:i1]),
    Event(:input_spikes, 0.0, 152.0, objects[:i2]),
    Event(:input_spikes, 0.0, 300.0, objects[:i1]),
    Event(:input_spikes, 0.0, 310.0, objects[:i2]),
]

logger=simulate!(net, input)

syn1 = get_trace(:syn1, logger.data)
syn2 = get_trace(:syn2, logger.data)
n    = get_trace(:n, logger.data)

steps!(ax1, [0;syn1.t;450], 1 .+ 0.9 .* [0;Int.(syn1.state);0], fill=color_1)
steps!(ax1, [0;syn2.t;450], 2 .+ 0.9 .* [0;Int.(syn2.state);0], fill=color_2)
steps!(ax1, [0;n.t;450], 0.45 .* [0;Int.(n.state);0] , fill=color_3_50)

################################################################################
config = """
refractory_duration: 5.01
inputs:
- id: i1
- id: i2
outputs:
- id: o
neurons:
- id: n
  branches: 
    - id: seg
synapses:
- {id: syn1, source: i1, target: seg}
- {id: syn2, source: i2, target: n}
- {id: syn_out, source: n, target: o, spike_duration: 0, delay: 0}
"""

(net,objects) = load_network(YAML_source=config)

input=[
    Event(:input_spikes, 0.0,  10.0, objects[:i1]),
    Event(:input_spikes, 0.0,   0.0, objects[:i2]),
    Event(:input_spikes, 0.0, 150.0, objects[:i1]),
    Event(:input_spikes, 0.0, 152.0, objects[:i2]),
    Event(:input_spikes, 0.0, 300.0, objects[:i1]),
    Event(:input_spikes, 0.0, 310.0, objects[:i2]),
]

logger=simulate!(net, input)

syn1 = get_trace(:syn1, logger.data)
syn2 = get_trace(:syn2, logger.data)
seg = get_trace(:seg, logger.data)
n    = get_trace(:n, logger.data)

steps!(ax2, [0;syn1.t;450], 1 .+ 0.9 .* [0;Int.(syn1.state);0], fill=color_1)
steps!(ax2, [0;seg.t;450], 1 .+ 0.45 .* [0;Int.(seg.state);0], fill=color_1_50)
steps!(ax2, [0;syn2.t;450], 2 .+ 0.9 .* [0;Int.(syn2.state);0], fill=color_2)
steps!(ax2, [0;n.t;450], 0.45 .* [0;Int.(n.state);0] , fill=color_3_50)

################################################################################

config = """
refractory_duration: 5.01
inputs:
- id: i1
- id: i2
outputs:
- id: o
neurons:
- id: n
  θ_syn: 1
  θ_seg: 2
  branches: 
    - id: seg1
    - id: seg2
synapses:
- {id: syn1, source: i1, target: seg1}
- {id: syn2, source: i2, target: seg2}
- {id: syn3, source: i1, target: n}
- {id: syn4, source: i2, target: n}
- {id: syn_out, source: n, target: o, spike_duration: 0, delay: 0}
"""

(net,objects) = load_network(YAML_source=config)

input=[
    Event(:input_spikes, 0.0,  10.0, objects[:i1]),
    Event(:input_spikes, 0.0,   0.0, objects[:i2]),
    Event(:input_spikes, 0.0, 150.0, objects[:i1]),
    Event(:input_spikes, 0.0, 152.0, objects[:i2]),
    Event(:input_spikes, 0.0, 300.0, objects[:i1]),
    Event(:input_spikes, 0.0, 310.0, objects[:i2]),
]

logger=simulate!(net, input)

syn1 = get_trace(:syn1, logger.data)
syn2 = get_trace(:syn2, logger.data)
seg1 = get_trace(:seg1, logger.data)
seg2 = get_trace(:seg2, logger.data)
n    = get_trace(:n, logger.data)

steps!(ax3, [0;syn1.t;450], 1 .+ 0.9 .* [0;Int.(syn1.state);0], fill=color_1)
steps!(ax3, [0;seg1.t;450], 1 .+ 0.45 .* [0;Int.(seg1.state);0], fill=color_1_50)
steps!(ax3, [0;syn2.t;450], 2 .+ 0.9 .* [0;Int.(syn2.state);0], fill=color_2)
steps!(ax3, [0;seg2.t;450], 2 .+ 0.45 .* [0;Int.(seg2.state);0], fill=color_2_50)
steps!(ax3, [0;n.t;450], 0.45 .* [0;Int.(n.state);0] , fill=color_3_50)

################################################################################

config = """
refractory_duration: 1.01
inputs:
- id: i1
- id: i2
outputs:
- id: o
neurons:
- id: n
  θ_syn: 1
  θ_seg: 2
  branches: 
    - id: seg1
    - id: seg2
synapses:
- {id: syn1, source: i1, target: seg1}
- {id: syn2, source: i2, target: seg2}
- {id: syn3, source: i1, target: n}
- {id: syn4, source: i2, target: n}
- {id: syn_out, source: n, target: o, spike_duration: 0, delay: 0}
"""

(net,objects) = load_network(YAML_source=config)

input=[
    Event(:input_spikes, 0.0,  10.0, objects[:i1]),
    Event(:input_spikes, 0.0,   0.0, objects[:i2]),
    Event(:input_spikes, 0.0, 150.0, objects[:i1]),
    Event(:input_spikes, 0.0, 152.0, objects[:i2]),
    Event(:input_spikes, 0.0, 300.0, objects[:i1]),
    Event(:input_spikes, 0.0, 310.0, objects[:i2]),
]

logger=simulate!(net, input)

syn1 = get_trace(:syn1, logger.data)
syn2 = get_trace(:syn2, logger.data)
seg1 = get_trace(:seg1, logger.data)
seg2 = get_trace(:seg2, logger.data)
n    = get_trace(:n, logger.data)

steps!(ax4, [0;syn1.t;450], 1 .+ 0.9 .* [0;Int.(syn1.state);0], fill=color_1)
steps!(ax4, [0;seg1.t;450], 1 .+ 0.45 .* [0;Int.(seg1.state);0], fill=color_1_50)
steps!(ax4, [0;syn2.t;450], 2 .+ 0.9 .* [0;Int.(syn2.state);0], fill=color_2)
steps!(ax4, [0;seg2.t;450], 2 .+ 0.45 .* [0;Int.(seg2.state);0], fill=color_2_50)
steps!(ax4, [0;n.t;450], 0.45 .* [0;Int.(n.state);0] , fill=color_3_50)

save("coincidence.svg", fig)
