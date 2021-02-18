using ADSP, Test
#ENV["JULIA_DEBUG"] = ADSP

@testset "Test delay & refractoriness" begin
    config = """
    spike_duration: 4.99
    inputs:
    - id: i1
    outputs:
    - id: o1
    neurons:
    - id: n1
    synapses:
    - {id: syn1, source: i1, target: n1}
    - {id: syn2, source: n1, target: o1, spike_duration: 0}
    """
    
    (net,objects) = load_network(YAML_source=config)
    
    inp=Event{Float64,Symbol,<:Any}[Event(1.0,:activate,objects[:i1])]

    logger=simulate!(net, inp, 10.0; show_progress=false)
    # should spike once every 1ms for the whole spike duration of 5ms with 2ms delay 
    # (due to the 2 synapses) after the input spike at 1ms, i.e. after 3ms
    @test length(logger.t)==5 && logger.t ≈ [3.0,4.0,5.0,6.0,7.0]

    # if the refractoriness = 5.0ms > the duration of a spike, the neuron should only spike once 
    config = """
    spike_duration: 4.99
    inputs:
    - id: i1
    outputs:
    - id: o1
    neurons:
    - id: n1
      refractory_duration: 5.0 # <-----------------
    synapses:
    - {id: syn1, source: i1, target: n1}
    - {id: syn2, source: n1, target: o1, spike_duration: 0}
    """
    
    (net,objects) = load_network(YAML_source=config)
    inp=Event{Float64,Symbol,<:Any}[Event(1.0,:activate,objects[:i1])]
    logger=simulate!(net, inp, 10.0; show_progress=false)
    @test length(logger.t)==1 && logger.t ≈ [3.0]
end

@testset "Reset test" begin
    config = """
    refractory_duration: 10
    plateau_duration: 10
    spike_duration: 10
    inputs:
    - id: i1
    outputs:
    - id: o1
    neurons:
    - id: n1
      branches:
        - id: s1 
          branches:
            - id: s2
              branches: 
                - id: s3
    synapses:
    - {id: syn1, source: i1, target: n1}
    - {id: syn2, source: i1, target: s1}
    - {id: syn3, source: i1, target: s2}
    - {id: syn4, source: i1, target: s3}
    - {id: syn5, source: n1, target: o1, spike_duration: 0}
    """
    
    (net,objects) = load_network(YAML_source=config)
    inp=Event{Float64,Symbol,<:Any}[Event(1.0,:activate,objects[:i1])]
    logger=simulate!(net, inp, 10.0; show_progress=false)
    @test length(logger.t)==1 && logger.t ≈ [3.0]
    logger=simulate!(net, inp, 10.0; show_progress=false, reset=false)
    @test isempty(logger.t)
    logger=simulate!(net, inp, 10.0; show_progress=false)
    @test length(logger.t)==1 && logger.t ≈ [3.0]
end

@testset "Test cascades" begin
  config = """
  spike_duration: 0.99
  plateau_duration: 5.0
  inputs:
  - id: i1
  - id: i2
  - id: i3
  - id: i4
  outputs:
  - id: o1
  neurons:
  - id: n1
    branches:
      - id: s1 
        branches:
          - id: s2
            branches: 
              - id: s3
  synapses:
  - {id: syn1, source: i1, target: n1}
  - {id: syn2, source: i2, target: s1}
  - {id: syn3, source: i3, target: s2}
  - {id: syn4, source: i4, target: s3}
  - {id: syn5, source: n1, target: o1, spike_duration: 0}
  """

  (net,objects) = load_network(YAML_source=config)

  cascade_trials=[
    (
      summary="Slow sequence in wrong direction",
      input=Event{Float64,Symbol,<:Any}[
        Event(1.0,:activate,objects[:i1]), # should do nothing
        Event(2.0,:activate,objects[:i2]), # should do nothing
        Event(3.0,:activate,objects[:i3]), # should do nothing
        Event(4.0,:activate,objects[:i4]), # should trigger a plateau on s3, but no spike!
      ],
      output=[]
    ),
    (
      summary="Re-triggering the neuron",
      input=Event{Float64,Symbol,<:Any}[
        Event(0.0,:activate,objects[:i1]), # should do nothing
        Event(0.2,:activate,objects[:i1]), # should do nothing
        Event(0.4,:activate,objects[:i1]), # should do nothing
      ],
      output=[]
    ),
    (
      summary="Quick sequence in wrong direction",
      input=Event{Float64,Symbol,<:Any}[
        Event(0.0,:activate,objects[:i1]), # should do nothing
        Event(0.2,:activate,objects[:i2]), # should do nothing
        Event(0.4,:activate,objects[:i3]), # should do nothing
        Event(0.6,:activate,objects[:i4]), # triggers a spike at 1.6 -> spike on syn5 at 2.6
      ],
      output=[2.6]
    ),
    (
      summary="Slow sequence in correct direction",
      input=Event{Float64,Symbol,<:Any}[
        Event(0.0,:activate,objects[:i4]), # should do nothing
        Event(2.0,:activate,objects[:i3]), # should do nothing
        Event(4.0,:activate,objects[:i2]), # should do nothing
        Event(6.0,:activate,objects[:i1]), # trigger a spike at 7.0 -> spike on syn5 at 8.0
        ],
      output=[8.0]
    ),
    (
      summary="Too slow sequence in correct direction",
      input=Event{Float64,Symbol,<:Any}[
        Event(0.0,:activate,objects[:i4]), # should do nothing
        Event(2.0,:activate,objects[:i3]), # should do nothing
        Event(4.0,:activate,objects[:i2]), # should do nothing
        Event(10.0,:activate,objects[:i1]), # should do nothing
        ],
      output=[]
    )      
  ]

  @testset "Test cascade: $(trial.summary)" for trial in cascade_trials
    logger=simulate!(net, trial.input; show_progress=false)
    @test length(logger.t) == length(trial.output) && logger.t ≈ Vector{Float64}(trial.output)
  end
end



@testset "Test multi-branch" begin
    config = """
    refractory_duration: 10
    inputs:
    - id: i1
    - id: i2
    - id: i3
    outputs:
    - id: o1
    neurons:
    - id: n1
      branches:
        - id: s1 
        - id: s2
    synapses:
    - {id: syn1, source: i1, target: n1}
    - {id: syn2, source: i2, target: s1}
    - {id: syn3, source: i3, target: s2}
    - {id: syn4, source: n1, target: o1, spike_duration: 0}
    """
    
    (net,objects) = load_network(YAML_source=config)

    trials=[
      (
        summary="Branch 1",
        input=Event{Float64,Symbol,<:Any}[
          Event(0.0,:activate,objects[:i2]), # should do nothing
          Event(1.0,:activate,objects[:i1]), # should trigger a spike!
        ],
        output=[3.0]
      ),
      (
        summary="Branch 2",
        input=Event{Float64,Symbol,<:Any}[
          Event(0.0,:activate,objects[:i3]), # should do nothing
          Event(1.0,:activate,objects[:i1]), # should trigger a spike!
        ],
        output=[3.0]
      ),
      (
        summary="No branch",
        input=Event{Float64,Symbol,<:Any}[
          Event(1.0,:activate,objects[:i3]), # should do nothing
          Event(2.0,:activate,objects[:i2]), # should not trigger a spike!
        ],
        output=[]
      ),
    ]
    
    @testset "Test cascade: $(trial.summary)" for trial in trials
        logger=simulate!(net, trial.input; show_progress=false)
        @test length(logger.t) == length(trial.output) && logger.t ≈ Vector{Float64}(trial.output)
    end
end