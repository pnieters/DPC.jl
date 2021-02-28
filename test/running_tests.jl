#ENV["JULIA_DEBUG"] = "ADSP"
using ADSP, Test

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
    
    inp=[Event(:input_spikes, 0.0, 1.0, objects[:i1])]

    logger=simulate!(net, inp, 10.0; show_progress=false, logger! = Logger(net; filter=(t,tp,id,x)->id==:o1))
    # should spike once every 1ms for the whole spike duration of 5ms with 2ms delay 
    # (due to the 2 synapses) after the input spike at 1ms, i.e. after 3ms
    @test length(logger.data.t)==5 && logger.data.t ≈ [3.0,4.0,5.0,6.0,7.0]

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
    inp=[Event(:input_spikes, 0.0, 1.0, objects[:i1])]
    logger=simulate!(net, inp, 10.0; show_progress=false, logger! = Logger(net; filter=(t,tp,id,x)->id==:o1))
    @test length(logger.data.t)==1 && logger.data.t ≈ [3.0]
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
    inp=[Event(:input_spikes, 0.0, 1.0, objects[:i1])]
    logger=simulate!(net, inp, 10.0; show_progress=false, logger! = Logger(net; filter=(t,tp,id,x)->id==:o1))
    @test length(logger.data.t)==1 && logger.data.t ≈ [3.0]
    logger=simulate!(net, inp, 10.0; show_progress=false, reset=false, logger! = Logger(net; filter=(t,tp,id,x)->id==:o1))
    @test isempty(logger.data.t)
    logger=simulate!(net, inp, 10.0; show_progress=false, logger! = Logger(net; filter=(t,tp,id,x)->id==:o1))
    @test length(logger.data.t)==1 && logger.data.t ≈ [3.0]
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
      input=[
        Event(:input_spikes, 0.0, 1.0, objects[:i1]), # should do nothing
        Event(:input_spikes, 0.0, 2.0, objects[:i2]), # should do nothing
        Event(:input_spikes, 0.0, 3.0, objects[:i3]), # should do nothing
        Event(:input_spikes, 0.0, 4.0, objects[:i4]), # should trigger a plateau on s3, but no spike!
      ],
      output=[]
    ),
    (
      summary="Re-triggering the neuron",
      input=[
        Event(:input_spikes, 0.0, 0.0, objects[:i1]), # should do nothing
        Event(:input_spikes, 0.0, 0.2, objects[:i1]), # should do nothing
        Event(:input_spikes, 0.0, 0.4, objects[:i1]), # should do nothing
      ],
      output=[]
    ),
    (
      summary="Quick sequence in wrong direction",
      input=[
        Event(:input_spikes, 0.0, 0.0, objects[:i1]), # should do nothing
        Event(:input_spikes, 0.0, 0.2, objects[:i2]), # should do nothing
        Event(:input_spikes, 0.0, 0.4, objects[:i3]), # should do nothing
        Event(:input_spikes, 0.0, 0.6, objects[:i4]), # triggers a spike at 1.6 -> spike on syn5 at 2.6
      ],
      output=[2.6]
    ),
    (
      summary="Slow sequence in correct direction",
      input=[
        Event(:input_spikes, 0.0, 0.0, objects[:i4]), # should do nothing
        Event(:input_spikes, 0.0, 2.0, objects[:i3]), # should do nothing
        Event(:input_spikes, 0.0, 4.0, objects[:i2]), # should do nothing
        Event(:input_spikes, 0.0, 6.0, objects[:i1]), # trigger a spike at 7.0 -> spike on syn5 at 8.0
        ],
      output=[8.0]
    ),
    (
      summary="Too slow sequence in correct direction",
      input=[
        Event(:input_spikes, 0.0, 0.0, objects[:i4]), # should do nothing
        Event(:input_spikes, 0.0, 2.0, objects[:i3]), # should do nothing
        Event(:input_spikes, 0.0, 4.0, objects[:i2]), # should do nothing
        Event(:input_spikes, 0.0, 10.0, objects[:i1]), # should do nothing
        ],
      output=[]
    )      
  ]

  @testset "Test cascade: $(trial.summary)" for trial in cascade_trials
    logger=simulate!(net, trial.input; show_progress=false, logger! = Logger(net; filter=(t,tp,id,x)->id==:o1))
    @test length(logger.data.t) == length(trial.output) && logger.data.t ≈ Vector{Float64}(trial.output)
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
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i2]), # should do nothing
          Event(:input_spikes, 0.0, 1.0, objects[:i1]), # should trigger a spike!
        ],
        output=[3.0]
      ),
      (
        summary="Branch 2",
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i3]), # should do nothing
          Event(:input_spikes, 0.0, 1.0, objects[:i1]), # should trigger a spike!
        ],
        output=[3.0]
      ),
      (
        summary="No branch",
        input=[
          Event(:input_spikes, 0.0, 1.0, objects[:i3]), # should do nothing
          Event(:input_spikes, 0.0, 2.0, objects[:i2]), # should not trigger a spike!
        ],
        output=[]
      ),
    ]
    
    @testset "Test cascade: $(trial.summary)" for trial in trials
        logger=simulate!(net, trial.input; show_progress=false, logger! = Logger(net; filter=(t,tp,id,x)->id==:o1))
        @test length(logger.data.t) == length(trial.output) && logger.data.t ≈ Vector{Float64}(trial.output)
    end
end


@testset "Test inhibition at soma" begin
    config = """
    spike_duration: 4.99
    inputs:
    - id: i1
    - id: i2
    - id: i3
    outputs:
    - id: o1
    neurons:
    - id: n1
      refractory_duration: 5.0 # <-----------------
    synapses:
    - {id: syn1, source: i1, target: n1}
    - {id: syn2, source: i2, target: n1, weight: -1}
    - {id: syn3, source: i3, target: n1, weight: -1, spike_duration: 10}
    - {id: syn4, source: n1, target: o1, spike_duration: 0}
    """
    
    (net,objects) = load_network(YAML_source=config)
    
    trials=[
      (
        summary="Inhibit after spike",
        input=[
          Event(:input_spikes, 0.0, 1.0, objects[:i1]),
          Event(:input_spikes, 0.0, 1.5, objects[:i2])
        ],
        output=[3.0]
      ),
      (
        summary="Inhibit before spike",
        input=[
          Event(:input_spikes, 0.0, 1.0, objects[:i2]),
          Event(:input_spikes, 0.0, 1.5, objects[:i1])
        ],
        output=[7.99]
      ),
      (
        summary="Inhibit with long EPSP",
        input=[
          Event(:input_spikes, 0.0, 1.0, objects[:i3]),
          Event(:input_spikes, 0.0, 1.5, objects[:i1])
        ],
        output=[]
      ),
      (
        summary="Inhibition can prevent rebound spike",
        input=[
          Event(:input_spikes, 0.0, 1.0, objects[:i3]),
          Event(:input_spikes, 0.0, 1.0, objects[:i2]),
          Event(:input_spikes, 0.0, 1.5, objects[:i1])
        ],
        output=[]
      ),
    ]
    @testset "$(trial.summary)" for trial in trials
      logger=simulate!(net, trial.input; show_progress=false, logger! = Logger(net; filter=(t,tp,id,x)->id==:o1))
      @test length(logger.data.t) == length(trial.output) && logger.data.t ≈ Vector{Float64}(trial.output)
    end
end

#####

@testset "Test inhibition on multiple branches" begin
    config = """
    refractory_duration: 10
    inputs:
    - id: i1
    - id: i2
    - id: i3
    - id: i4
    - id: i5
    - id: i6
    - id: i7
    - id: i8
    outputs:
    - id: o1
    neurons:
    - id: n1
      branches:
        - id: s0
          branches:
            - id: s1 
            - id: s2
    synapses:
    - {id: syn1, source: i1, target: n1}
    - {id: syn2, source: i2, target: n1, weight: -1, spike_duration: 10}
    - {id: syn3, source: i3, target: s1}
    - {id: syn4, source: i4, target: s1, weight: -1}
    - {id: syn5, source: i5, target: s2}
    - {id: syn6, source: i6, target: s2, weight: -1}
    - {id: syn7, source: i7, target: s0}
    - {id: syn8, source: i8, target: s0, weight: -1}
    - {id: syn9, source: n1, target: o1, spike_duration: 0}
    """

    
    (net,objects) = load_network(YAML_source=config)

    trials=[
      (
        summary="Cascade",
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i3]), # trigger plateau on s1 
          Event(:input_spikes, 0.0, 7.5, objects[:i7]), # trigger plateau on s0
          Event(:input_spikes, 0.0, 10.0, objects[:i1]), # should trigger a spike!
        ],
        output=[12.0]
      ),
      (
        summary="Cascade interrupted 1",
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i3]), # trigger plateau on s1
          Event(:input_spikes, 0.0, 5.0, objects[:i4]), # end plateau on s1
          Event(:input_spikes, 0.0, 7.5, objects[:i7]), # fails to trigger plateau on s0
          Event(:input_spikes, 0.0, 10.0, objects[:i1]), # should not trigger a spike!
        ],
        output=[]
      ),
      (
        summary="Cascade not interrupted 2",
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i3]), # trigger plateau in s1 
          Event(:input_spikes, 0.0, 5.0, objects[:i7]), # trigger plateau in s0 
          Event(:input_spikes, 0.0, 7.5, objects[:i4]), # end plateau in top s1 --> doesn't turns off s0
          Event(:input_spikes, 0.0, 10.0, objects[:i1]), # should trigger a spike!
        ],
        output=[12.0]
      ),
      (
        summary="Cascade interrupted 3",
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i3]), # trigger plateau
          Event(:input_spikes, 0.0, 5.0, objects[:i7]), # trigger plateau
          Event(:input_spikes, 0.0, 7.5, objects[:i8]), # end plateau in lower segment
          Event(:input_spikes, 0.0, 10.0, objects[:i1]), # should not trigger a spike!
        ],
        output=[]
      ),
      (
        summary="Cascade interrupted 4",
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i3]), # trigger plateau
          Event(:input_spikes, 0.0, 5.0, objects[:i7]), # trigger plateau
          Event(:input_spikes, 0.0, 7.5, objects[:i2]), # inhibit soma
          Event(:input_spikes, 0.0, 10.0, objects[:i1]), # should not trigger a spike!
        ],
        output=[]
      ),
      (
        summary="Cross-branch uninterrupted",
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i3]), # trigger plateau on s1
          Event(:input_spikes, 0.0, 1.0, objects[:i5]), # trigger plateau on s2
          Event(:input_spikes, 0.0, 2.0, objects[:i7]), # trigger plateau on s0
          Event(:input_spikes, 0.0, 10.0, objects[:i1]), # should trigger a spike!
        ],
        output=[12.0]
      ),
      (
        summary="Cross-branch interrupted",
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i3]), # trigger plateau on s1
          Event(:input_spikes, 0.0, 1.0, objects[:i5]), # trigger plateau on s2
          Event(:input_spikes, 0.0, 2.0, objects[:i7]), # trigger plateau on s0
          Event(:input_spikes, 0.0, 3.0, objects[:i8]), # end plateau on s0
          Event(:input_spikes, 0.0, 10.0, objects[:i1]), # should not trigger a spike!
        ],
        output=[]
      ),
      (
        summary="Cross-branch not interrupted 1",
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i3]), # trigger plateau on s1
          Event(:input_spikes, 0.0, 1.0, objects[:i5]), # trigger plateau on s2
          Event(:input_spikes, 0.0, 2.0, objects[:i7]), # trigger plateau on s0
          Event(:input_spikes, 0.0, 3.0, objects[:i4]), # end plateau on s1
          Event(:input_spikes, 0.0, 10.0, objects[:i1]), # should still trigger a spike!
        ],
        output=[12.0]
      ),
      (
        summary="Cross-branch not interrupted 2",
        input=[
          Event(:input_spikes, 0.0, 0.0, objects[:i3]), # trigger plateau on s1
          Event(:input_spikes, 0.0, 1.0, objects[:i5]), # trigger plateau on s2
          Event(:input_spikes, 0.0, 2.0, objects[:i7]), # trigger plateau on s0
          Event(:input_spikes, 0.0, 3.0, objects[:i4]), # end plateau on s1
          Event(:input_spikes, 0.0, 5.0, objects[:i6]), # end plateau on s2
          Event(:input_spikes, 0.0, 10.0, objects[:i1]), # should not trigger a spike!
        ],
        output=[12.0]
      )
    ]

    @testset "Test cascade: $(trial.summary)" for trial in trials
        logger=simulate!(net, trial.input; show_progress=false, logger! = Logger(net; filter=(t,tp,id,x)->id==:o1))
        @test length(logger.data.t) == length(trial.output) && logger.data.t ≈ Vector{Float64}(trial.output)
    end
end

@testset "Don't forget reset" begin
    #This test for a bug, where a plateau-reset would not happen if plateau conditions are still met at the time where 
    #the plateau should end. No new end would be scheduled, therefore the segment would stay on indefinitely.

    config = """
    plateau_duration: 100
    refractory_duration: 5.01
    inputs: 
    - id: i1
    - id: i2
    neurons:
    - id: n
      branches:
        - id: seg
    synapses:
    - {id: syn1, source: i1, target: seg}
    - {id: syn2, source: i2, target: n}
    """

    (net,objects) = load_network(YAML_source=config)
    input=[
      Event(:input_spikes, 0.0, 0.0, objects[:i1]), # trigger plateau on seg
      Event(:input_spikes, 0.0, 99.0, objects[:i1]), # satisfy plateau conditions at time of plateau-end
      Event(:input_spikes, 0.0, 1000.0, objects[:i2]), # this would trigger a spike much later
    ]
    logger=simulate!(net, input; show_progress=false, logger! = Logger(net; filter=(t,tp,id,x)->tp==:spikes))
    @test isempty(logger.data)
end

@testset "Prevent multiple triggers" begin
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
        Event(:input_spikes, 0.0,  0.0, objects[:i1]),
        Event(:input_spikes, 0.0,  2.0, objects[:i2]),
    ]

    logger=simulate!(net, input, 10.0; show_progress=false, logger! = Logger(net; filter=(t,tp,id,x)->id==:n))
    n = filter(x->(x.object==:n && x.event==:spikes), logger.data)[!, :t]
    @test length(logger.data.t)==1 && n ≈ [3.0]
end
