@testset "Loading from YAML tests" begin
    test_net="""
spike_duration: 0.1
T_refrac: 0.2
delay: 0.3
plateau_duration: 50.0

# Comments should work, tool
inputs:
    - id: i1 # even inline comments
    - id: i2
    - id: i3

random_unsupported_network_stuff1: nothing_of_relevance
random_unsupported_network_stuff2: nothing_of_relevance

neurons:
    - id: n1
      random_unsupported_neuron_stuff: nothing_of_relevance
    - id: n2
      T_refrac: 0.4
      θ_syn: 2
      θ_seg: 2
      branches:
        - id: s1
          θ_syn: 10
        - id: s2
          branches:
            - id: s21
              plateau_duration: 20.0
              random_unsupported_segment_stuff: nothing_of_relevance
              θ_seg: 0
            - id: s22

synapses:
    - {id: syn1, source: i1, target: n1}
    - {id: syn2, source: i2, target: s22}
    - {id: syn3, source: i3, target:  s2, spike_duration: 3.0, random_unsupported_synapse_stuff: nothing_of_relevance}
    - {id: syn4, source: n1, target:  s1, delay: 10.0, weight: 2}
    """

    (net,obj) = (
        @test_logs (
                :warn, "The following network parameters could not be parsed: random_unsupported_network_stuff1, random_unsupported_network_stuff2"
            ) (
                :warn, "The following neuron parameters could not be parsed: random_unsupported_neuron_stuff"
            ) (
                :warn, "The following segment parameters could not be parsed: random_unsupported_segment_stuff"
            ) (
                :warn, "The following synapse parameters could not be parsed: random_unsupported_synapse_stuff"
            ) load_network(YAML_source=test_net)
        )
    
    @testset "Network structure tests" begin
        @test net.neurons[1] == obj[:n1]
        @test net.neurons[2] == obj[:n2]
        @test net.neurons[2].next_upstream[1] == obj[:s1]
        @test net.neurons[2].next_upstream[2].next_upstream[1] == obj[:s21]
        @test net.inputs[1] == obj[:i1]
        @test net.inputs[1].synapses[1] == obj[:syn1]
        @test net.neurons[1].synapses[1] == obj[:syn4]
    end

    @testset "Neuron parameters tests" begin
        @test obj[:n1].id == :n1
        @test obj[:n1].T_refrac == 0.2
        @test obj[:n1].θ_syn == 1
        @test obj[:n1].θ_seg == 1
        @test obj[:n2].T_refrac == 0.4
        @test obj[:n2].θ_syn == 2
        @test obj[:n2].θ_seg == 2
    end

    @testset "Segment parameters tests" begin
        @test obj[:s1].id == :s1
        @test obj[:s1].plateau_duration == 50.0
        @test obj[:s1].θ_syn == 10
        @test obj[:s1].θ_seg == 1
        @test obj[:s21].plateau_duration == 20.0
        @test obj[:s21].θ_syn == 1
        @test obj[:s21].θ_seg == 0
    end

    @testset "Synapse parameters tests" begin
        @test obj[:syn3].id == :syn3
        @test obj[:syn3].target == obj[:s2]
        @test obj[:syn3].delay == 0.3
        @test obj[:syn3].spike_duration == 3.0
        @test obj[:syn3].weight == 1
        @test obj[:syn4].target == obj[:s1]
        @test obj[:syn4].delay == 10.0
        @test obj[:syn4].weight == 2
    end

    @testset "Saving tests" begin
        tmp_str = save_network(net)
        (net2,obj2) = @test_logs  load_network(YAML_source=tmp_str)
        @test all(pairs(obj2)) do (key,value)
            key ∈ keys(obj)
        end
    end
end

