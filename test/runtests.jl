using ADSP, YAML

# ###############################
# net = Network()

# i1 = Input(:i1, net)
# i2 = Input(:i2, net)
# i3 = Input(:i3, net)

# n = Neuron(:test, net)


# s1 = Segment(:s1, n)
# s2 = Segment(:s2, n)
# s21 = Segment(:s21, s2)
# s22 = Segment(:s22, s2)

# syn1 = Synapse(:syn1, i1, s21)
# syn2 = Synapse(:syn2, i2, s22)
# syn3 = Synapse(:syn3, i3, s2)
# syn4 = Synapse(:syn4, n, s1; delay=10.0)




net="""
inputs:
    - id: i1
    - id: i2
    - id: i3

neurons:
    - id: n
      branches:
        - id: s1
        - id: s2
          branches:
            - id: s21
            - id: s22

synapses:
    - {id: syn1, source: i1, target: s21}
    - {id: syn2, source: i2, target: s22}
    - {id: syn3, source: i3, target:  s2}
    - {id: syn4, source:  n, target:  s1, delay: 10.0}
"""

(net,obj) = load_network(YAML_source=net)
inp=Event{Float64,Symbol,<:Any}[Event(1.0,:spike_start,obj[:i1])]
simulate!(net, inp; callback! = println)
