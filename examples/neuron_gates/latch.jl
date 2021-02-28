#ENV["JULIA_DEBUG"] = "ADSP"
using ADSP, Makie

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
- {id: syn2, source: i_reset, target: n_storage, weight: -10}
- {id: syn3, source: n_storage, target: n_storage}
"""

(net,objects) = load_network(YAML_source=config)

input=sort!([
    [Event(:input_spikes, 0.0,  t, objects[:i_set]) for t ∈ [105.0, 605.0]];
    [Event(:input_spikes, 0.0,  t, objects[:i_reset]) for t ∈ [400.0, 800.0]];
])

logger=simulate!(net, input)

n    = filter(x->(x.object==:n_storage && x.event==:spikes), logger.data)[!, :t]
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
- {id: syn2, source: i_reset, target: n_storage, weight: -10}
- {id: syn3, source: n_storage, target: n_storage}
- {id: syn4, source: n_storage, target: seg_readout}
- {id: syn5, source: i_read, target: n_readout}
"""

(net,objects) = load_network(YAML_source=config)

input=sort!([
    [Event(:input_spikes, 0.0,  t, objects[:i_set]) for t ∈ [105.0, 1105.0]];
    [Event(:input_spikes, 0.0,  t, objects[:i_reset]) for t ∈ [400.0, 1400.0]];
    [Event(:input_spikes, 0.0,  t, objects[:i_read]) for t ∈ [150.0, 300.0, 450.0, 600.0, 750.0, 900.0, 1050.0, 1300.0, 1450.0, 1600.0]];
])

logger=simulate!(net, input)

n_storage   = filter(x->(x.object==:n_storage && x.event==:spikes), logger.data)[!, :t]
n_readout   = filter(x->(x.object==:n_readout && x.event==:spikes), logger.data)[!, :t]
################################################################################
