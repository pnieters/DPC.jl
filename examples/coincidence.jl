using ADSP, Makie
#include(joinpath(@__DIR__, "utils.jl"))
include("utils.jl")

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

syn1 = filter(x->(x.object==:syn1), logger.data)
syn2 = filter(x->(x.object==:syn2), logger.data)
n    = filter(x->(x.object==:n), logger.data)

steps([0;syn1.t;320], 1 .+ 0.9 .* [0;Int.(syn1.state);0], fill=:red)
steps!([0;syn2.t;320], 2 .+ 0.9 .* [0;Int.(syn2.state);0], fill=:blue)
steps!([0;n.t;320], 0.9 .* [0;Int.(n.state);0], fill=:green)

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

syn1 = filter(x->(x.object==:syn1), logger.data)
syn2 = filter(x->(x.object==:syn2), logger.data)
n    = filter(x->(x.object==:n), logger.data)

steps([0;syn1.t;320], 1 .+ 0.9 .* [0;Int.(syn1.state);0], fill=:red)
steps!([0;syn2.t;320], 2 .+ 0.9 .* [0;Int.(syn2.state);0], fill=:blue)
steps!([0;n.t;320], 0.9 .* [0;Int.(n.state);0], fill=:green)

################################################################################

config = """
refractory_duration: 100.01
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

syn1 = filter(x->(x.object==:syn1), logger.data)
syn2 = filter(x->(x.object==:syn2), logger.data)
n    = filter(x->(x.object==:n), logger.data)

steps([0;syn1.t;320], 1 .+ 0.9 .* [0;Int.(syn1.state);0], fill=:red)
steps!([0;syn2.t;320], 2 .+ 0.9 .* [0;Int.(syn2.state);0], fill=:blue)
steps!([0;n.t;320], 0.9 .* [0;Int.(n.state);0], fill=:green)
