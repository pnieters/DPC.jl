using ADSP, Makie

config = """
refractory_duration: 5.01
inputs:
- id: i1
- id: i2
- id: sel1
- id: sel2
outputs:
neurons:
- id: n
  branches: 
    - id: seg1b
      branches: 
        - id: seg1a
    - id: seg2b
      branches: 
        - id: seg2a
synapses:
- {id: syn1, source: i1, target: seg1b}
- {id: syn2, source: i2, target: seg2b}
- {id: syn3, source: i1, target: n}
- {id: syn4, source: i2, target: n}
- {id: syn5, source: sel1, target: seg1a}
- {id: syn6, source: sel2, target: seg2a}
"""

(net,objects) = load_network(YAML_source=config)

input=sort!([
    [Event(:input_spikes, 0.0,  t, objects[:i1]) for t ∈ 0.0:50:1000.0];
    [Event(:input_spikes, 0.0,  t, objects[:i2]) for t ∈ 10.0:71:1000.0];
    [Event(:input_spikes, 0.0,  t, objects[:sel1]) for t ∈ [105.0, 805.0]];
    [Event(:input_spikes, 0.0,  t, objects[:sel2]) for t ∈ [400.0, 600.0]];
])

logger=simulate!(net, input)

syn1 = filter(x->(x.object==:syn1 && x.event==:epsp_starts), logger.data)[!, :t]
syn2 = filter(x->(x.object==:syn2 && x.event==:epsp_starts), logger.data)[!, :t]
n    = filter(x->(x.object==:n && x.event==:spikes), logger.data)[!, :t]

################################################################################

config = """
refractory_duration: 5.01
inputs:
- id: i1
- id: i2
- id: sel1
- id: sel2
outputs:
neurons:
- id: n1
  branches: 
    - id: seg1
- id: n2
  branches: 
    - id: seg2
synapses:
- {id: syn1, source: i1, target: n1}
- {id: syn2, source: i2, target: n2}
- {id: syn3, source: sel1, target: seg1}
- {id: syn4, source: sel2, target: seg2}
"""

(net,objects) = load_network(YAML_source=config)

logger=simulate!(net, input)

syn1 = filter(x->(x.object==:syn1 && x.event==:epsp_starts), logger.data)[!, :t]
syn2 = filter(x->(x.object==:syn2 && x.event==:epsp_starts), logger.data)[!, :t]
n1    = filter(x->(x.object==:n1 && x.event==:spikes), logger.data)[!, :t]
n2    = filter(x->(x.object==:n2 && x.event==:spikes), logger.data)[!, :t]

################################################################################


