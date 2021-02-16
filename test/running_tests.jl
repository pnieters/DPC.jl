using ADSP, Test
#ENV["JULIA_DEBUG"] = ADSP

@testset "Test delay & refractoriness" begin
    config = """
    spike_duration: 4.99
    inputs:
    - id: i1
    neurons:
    - id: n1
    - id: n2
    synapses:
    - {id: syn1, source: i1, target: n1}
    - {id: syn2, source: n1, target: n2}
    """
    
    (net,objects) = load_network(YAML_source=config)
    
    inp=Event{Float64,Symbol,<:Any}[Event(1.0,:activate,objects[:i1])]

    event_times = []
    simulate!(net, inp, 10.0; show_progress=false, callback! = ev -> (if (ev.value.id==:syn2 && ev.id==:activate) push!(event_times,ev.t) end))
    # should spike once every 1ms for the whole spike duration of 5ms with 2ms delay 
    # (due to the 2 synapses) after the input spike at 1ms, i.e. after 3ms
    @test event_times ≈ [3.0,4.0,5.0,6.0,7.0]

    # if the refractoriness = 5.0ms > the duration of a spike, the neuron should only spike once 
    config = """
    spike_duration: 4.99
    inputs:
    - id: i1
    neurons:
    - id: n1
      T_refrac: 5.0
    - id: n2
    synapses:
    - {id: syn1, source: i1, target: n1}
    - {id: syn2, source: n1, target: n2}
    """
    
    (net,objects) = load_network(YAML_source=config)
    inp=Event{Float64,Symbol,<:Any}[Event(1.0,:activate,objects[:i1])]
    event_times = []
    simulate!(net, inp, 10.0; show_progress=false, callback! = ev -> if ev.value.id==:syn2 && ev.id==:activate push!(event_times,ev.t) end)
    @test event_times ≈ [3.0]
end
