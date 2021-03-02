#ENV["JULIA_DEBUG"] = "ADSP"
using ADSP, CairoMakie
include("utils.jl")

fig = Figure(resolution = (1200, 900))
# ax1 = fig[1, 1] = Axis(fig, title = "Single neuron latch")
# ax2 = fig[2, 1] = Axis(fig, title = "Two neuron latch")

config = """
refractory_duration: 5.01
inputs:
- id: i_set
- id: i_reset
outputs:
neurons:
- id: n1
- id: n2
- id: n3
- id: n4
- id: n5
synapses:
- {id: syn_set_1, source: i_set, target: n1, weight: [1.0, 0.5]}
- {id: syn_set_2, source: i_set, target: n2, weight: [1.0, 0.5]}
- {id: syn_set_3, source: i_set, target: n3, weight: [1.0, 0.5]}
- {id: syn_set_4, source: i_set, target: n4, weight: [1.0, 0.5]}
- {id: syn_set_5, source: i_set, target: n5, weight: [1.0, 0.5]}
- {id: syn_reset_1, source: i_reset, target: n1, weight: -10, spike_duration: 10}
- {id: syn_reset_2, source: i_reset, target: n2, weight: -10, spike_duration: 10}
- {id: syn_reset_3, source: i_reset, target: n3, weight: -10, spike_duration: 10}
- {id: syn_reset_4, source: i_reset, target: n4, weight: -10, spike_duration: 10}
- {id: syn_reset_5, source: i_reset, target: n5, weight: -10, spike_duration: 10}
- {id: syn_recurrent_1_2, source: n1, target: n2, weight: [1.0, 0.5]}
- {id: syn_recurrent_1_3, source: n1, target: n3, weight: [1.0, 0.5]}
- {id: syn_recurrent_1_4, source: n1, target: n4, weight: [1.0, 0.5]}
- {id: syn_recurrent_1_5, source: n1, target: n5, weight: [1.0, 0.5]}
- {id: syn_recurrent_2_1, source: n2, target: n1, weight: [1.0, 0.5]}
- {id: syn_recurrent_2_3, source: n2, target: n3, weight: [1.0, 0.5]}
- {id: syn_recurrent_2_4, source: n2, target: n4, weight: [1.0, 0.5]}
- {id: syn_recurrent_2_5, source: n2, target: n5, weight: [1.0, 0.5]}
- {id: syn_recurrent_3_1, source: n3, target: n1, weight: [1.0, 0.5]}
- {id: syn_recurrent_3_2, source: n3, target: n2, weight: [1.0, 0.5]}
- {id: syn_recurrent_3_4, source: n3, target: n4, weight: [1.0, 0.5]}
- {id: syn_recurrent_3_5, source: n3, target: n5, weight: [1.0, 0.5]}
- {id: syn_recurrent_4_1, source: n4, target: n1, weight: [1.0, 0.5]}
- {id: syn_recurrent_4_2, source: n4, target: n2, weight: [1.0, 0.5]}
- {id: syn_recurrent_4_3, source: n4, target: n3, weight: [1.0, 0.5]}
- {id: syn_recurrent_4_5, source: n4, target: n5, weight: [1.0, 0.5]}
- {id: syn_recurrent_5_1, source: n5, target: n1, weight: [1.0, 0.5]}
- {id: syn_recurrent_5_2, source: n5, target: n2, weight: [1.0, 0.5]}
- {id: syn_recurrent_5_3, source: n5, target: n3, weight: [1.0, 0.5]}
- {id: syn_recurrent_5_4, source: n5, target: n4, weight: [1.0, 0.5]}
"""

(net,objects) = load_network(YAML_source=config, weight_type=BernoulliSynapseWeight{Float64})

