export Synapse, Segment, Neuron, Network, Input, Logger
using DataFrames: DataFrame


################################################################################
## model component definitions                                                ##
################################################################################

abstract type AbstractComponent{ID,T,WT,IT} end
abstract type NeuronOrSegmentOrOutput{ID,T,WT,IT} <: AbstractComponent{ID,T,WT,IT} end
abstract type NeuronOrSegment{ID,T,WT,IT} <: NeuronOrSegmentOrOutput{ID,T,WT,IT} end
abstract type AbstractNetwork{ID,T,WT,IT} <:AbstractComponent{ID,T,WT,IT} end

struct Segment{ID,T,WT,IT} <: NeuronOrSegment{ID,T,WT,IT}
    id::ID
    θ_syn::IT
    θ_seg::Int
    plateau_duration::T

    state_syn::Ref{IT}
    state_seg::Ref{Int}
    state::Ref{Symbol}
    
    next_downstream::NeuronOrSegment{ID,T,WT,IT}
    next_upstream::Vector{Segment{ID,T,WT,IT}}
    net::AbstractNetwork{ID,T,WT,IT}

    function Segment(id::ID, root::NeuronOrSegment{ID,T,WT,IT}; θ_syn=1, θ_seg=1, plateau_duration=root.net.default_plateau_duration) where {ID,T,WT,IT}
        this = new{ID,T,WT,IT}(id, θ_syn, θ_seg, plateau_duration, Ref(zero(IT)), Ref(0), Ref(:off), root, Segment{ID,T,WT,IT}[], root.net)
        push!(root.next_upstream, this)
        return this
    end
end

struct Synapse{ID,T,WT,IT} <: AbstractComponent{ID,T,WT,IT}
    id::ID
    target::NeuronOrSegmentOrOutput{ID,T,WT,IT}
    delay::T
    spike_duration::T
    weight::WT
    state::Ref{IT}
    on::Ref{Bool}

    net::AbstractNetwork{ID,T,WT,IT}

    function Synapse(
                id::ID, source, target::NeuronOrSegmentOrOutput{ID,T,WT,IT}; 
                delay::T = source.net.default_delay, spike_duration::T = source.net.default_spike_duration, weight::WT=one(IT)
            ) where {ID,T,WT,IT}
        this = new{ID,T,WT,IT}(id, target, delay, spike_duration, weight, Ref(zero(IT)), Ref(false), source.net)
        push!(source.synapses, this)
        return this
    end
end

struct Neuron{ID,T,WT,IT} <: NeuronOrSegment{ID,T,WT,IT}
    id::ID
    θ_syn::IT
    θ_seg::Int
    refractory_duration::T

    state_syn::Ref{IT}
    state_seg::Ref{Int}
    state::Ref{Symbol}

    next_upstream::Vector{Segment{ID,T,WT,IT}}
    synapses::Vector{Synapse{ID,T,WT,IT}}
    net::AbstractNetwork{ID,T,WT,IT}
    
    function Neuron(id::ID, net::AbstractComponent{ID,T,WT,IT}; θ_syn=one(IT), θ_seg=1, refractory_duration::T=net.default_refractory_duration) where {ID,T,WT,IT}
        this = new{ID,T,WT,IT}(id,θ_syn,θ_seg,refractory_duration,Ref(zero(IT)),Ref(0),Ref(:off),Segment{ID,T,WT,IT}[],Synapse{ID,T,WT,IT}[],net)
        push!(net.neurons, this)
        return this
    end
end

struct Input{ID,T,WT,IT} <: AbstractComponent{ID,T,WT,IT}
    id::ID
    synapses::Vector{Synapse{ID,T,WT,IT}}
    net::AbstractNetwork{ID,T,WT,IT}

    function Input(id::ID, net::AbstractComponent{ID,T,WT,IT}) where {ID,T,WT,IT}
        this = new{ID,T,WT,IT}(id::ID,Synapse{ID,T,WT,IT}[],net)
        push!(net.inputs, this)
        return this
    end
end

struct Output{ID,T,WT,IT} <: NeuronOrSegmentOrOutput{ID,T,WT,IT}
    id::ID
    net::AbstractNetwork{ID,T,WT,IT}
    state_syn::Ref{IT}

    function Output(id::ID, net::AbstractComponent{ID,T,WT,IT}) where {ID,T,WT,IT}
        this = new{ID,T,WT,IT}(id::ID, net, Ref(zero(IT)))
        push!(net.outputs, this)
        return this
    end
end

struct Network{ID,T,WT,IT} <: AbstractNetwork{ID,T,WT,IT}
    inputs::Vector{Input{ID,T,WT,IT}}
    outputs::Vector{Output{ID,T,WT,IT}}
    neurons::Vector{Neuron{ID,T,WT,IT}}

    default_refractory_duration::T
    default_delay::T
    default_spike_duration::T
    default_plateau_duration::T

    Network(;
        id_type::Type=Symbol, time_type::Type=Float64, weight_type::Type=Int, synaptic_input_type::Type=Int, 
        default_refractory_duration = 1.0, default_delay = 1.0, default_spike_duration = 5.0, default_plateau_duration = 100.0
    ) = new{id_type, time_type, weight_type, synaptic_input_type}(
            Input{id_type, time_type, weight_type, synaptic_input_type}[], 
            Output{id_type, time_type, weight_type, synaptic_input_type}[], 
            Neuron{id_type, time_type, weight_type, synaptic_input_type}[],
            default_refractory_duration, default_delay, default_spike_duration, default_plateau_duration
        )
end

################################################################################
## model behavior                                                             ##
################################################################################
function reset!(net::Network{ID,T,WT,IT}) where {ID,T,WT,IT}
    foreach(reset!, net.inputs)
    foreach(reset!, net.outputs)
    foreach(reset!, net.neurons)
end

function reset!(x::Neuron{ID,T,WT,IT}) where {ID,T,WT,IT}
    x.state_syn[] = zero(IT)
    x.state_seg[] = 0
    x.state[] = :off

    foreach(reset!, x.synapses)
    foreach(reset!, x.next_upstream)
end

function reset!(x::Segment{ID,T,WT,IT}) where {ID,T,WT,IT}
    x.state_syn[] = zero(IT)
    x.state_seg[] = 0
    x.state[] = :off

    foreach(reset!, x.next_upstream)
end

function reset!(x::Input{ID,T,WT,IT}) where {ID,T,WT,IT}
    foreach(reset!, x.synapses)
end

function reset!(x::Output{ID,T,WT,IT}) where {ID,T,WT,IT}
    x.state_syn[] = zero(IT)
end

function reset!(x::Synapse{ID,T,WT,IT}) where {ID,T,WT,IT}
    x.on[]=false
    x.state[]=zero(IT)
end

"""Turn all upstream branches on."""
function backprop_on!(s::Segment)
    # only backprop if the segment is currently off; if it's on or clamped, we don't need to do anything
    if s.state[] == :off
        for s_up in s.next_upstream
            backprop_on!(s_up)
        end
    end

    # this segment is now clamped to the downstream and all upstream branches are on
    s.state_seg[] = length(s.next_upstream)
    s.state[] = :clamped
    return nothing
end

"""Turn all upstream branches off"""
function backprop_off!(s::Segment)
    for s_up in s.next_upstream
        backprop_off!(s_up)
    end

    # this segment is now off all so are all upstream branches
    s.state_seg[] = 0
    s.state[] = :off
    return nothing
end

"""Check if this segment was just turned on; if so, backpropagate!"""
function trigger!( t,  ev::Val{:activate}, s::Segment,  queue!, logger!)
    @debug "$(t): Triggered event $(ev) on object $(s)"
    # if the segment was already on or clamped, do nothing
    # if the segment is currently off, but the plateau conditions are satisfied, turn on
    if s.state[] == :off && s.state_syn[] >= s.θ_syn && (isempty(s.next_upstream) || s.state_seg[] >= s.θ_seg)
        # propagate to the upstream segments
        for s_up in s.next_upstream
            backprop_on!(s_up)
        end

        # the segment is now on
        s.state[] = :on
        s.next_downstream.state_seg[] += 1

        # turning on this segment might have also turned on the next_downstream segment
        trigger!( t,  ev, s.next_downstream,  queue!, logger!)

        # if we are still on (and not clamped by a downstream segment), don't forget to turn off this segment!
        if s.state[] == :on
            queue!(s.plateau_duration, :deactivate, s)
        end
    end
    return nothing
end

"""Check if this segment should switch off, since its plateau has timed out"""
function trigger!( t,  ev::Val{:deactivate}, s::Segment,  queue!, logger!)
    @debug "$(t): Triggered event $(ev) on object $(s)"
    # if the segment is clamped or off, do nothing
    # if the segment is currently on by itself, turn off
    if s.state[]==:on
        # propagate this info upwards
        backprop_off!(s)

        # update next_downstream segment
        s.next_downstream.state_seg[] -= 1
    end
    return nothing
end


"""Check if this Input was triggered to spike"""
function trigger!( t,  ev::Val{:activate}, n::Input,  queue!, logger!)
    @debug "$(t): Triggered event $(ev) on object $(n)"
    for s in n.synapses
        queue!(s.delay, :activate, s)
    end
    return nothing
end

"""Check if this neuron was triggered to spike"""
function trigger!( t,  ev::Val{:activate}, n::Neuron,  queue!, logger!)
    @debug "$(t): Triggered event $(ev) on object $(n)"
    # if the neuron was off but enough segments are on, trigger a spike
    if n.state[] == :off  && n.state_syn[] >= n.θ_syn && (isempty(n.next_upstream) || n.state_seg[] >= n.θ_seg)
        n.state[] == :refractory

        # trigger all synapses
        for s in n.synapses
            queue!(s.delay, :activate, s)
        end

        # don't forget to turn of :refractoriness
        queue!(n.refractory_duration, :deactivate, n)
    end
    return nothing
end

"""Check if this neuron was re-triggered to spike after refractoriness"""
function trigger!( t,  ev::Val{:deactivate}, n::Neuron,  queue!, logger!)
    @debug "$(t): Triggered event $(ev) on object $(n)"
    n.state[] = :off

    # check if the neuron should spike again, right away
    trigger!( t,  Val(:activate), n,  queue!, logger!)
end

"""Check if an incoming spike should trigger an EPSP"""
function trigger!( t,  ev::Val{:activate}, s::Synapse,  queue!, logger!)
    @debug "$(t): Triggered event $(ev) on object $(s)"
    if ~s.on[]
        s.on[] = true
        # set the synapse's state for this spike
        s.state[] = s.weight
        # inform the target segment about this new EPSP
        s.target.state_syn[] += s.state[]

        # don't forget to turn off EPSP
        queue!(s.spike_duration, :deactivate, s)

        # this EPSP might have triggered a plateau in the target
        trigger!( t,  Val(:activate), s.target,  queue!, logger!)
    end
end

"""Turn off the EPSP"""
function trigger!( t,  ev::Val{:deactivate}, s::Synapse,  queue!, logger!)
    @debug "$(t): Triggered event $(ev) on object $(s)"
    # only update if the synapse isn't already off (shouldn't happen!)
    if s.on[]
        # inform the target, that the EPSP is over
        s.target.state_syn[] -= s.state[]
        # reset the synapse's state
        s.state[] = zero(s.weight)
        s.on[] = false
    end
    return nothing
end

"""Trigger output - logs event by default"""
function trigger!( t,  ev::Val{X}, o::Output,  queue!, logger!) where X
    @debug "$(t): Triggered event $(ev) on object $(s)"
    logger!(t, X, o)
    return nothing
end


################################################################################
## Logging                                                                    ##
################################################################################

struct Logger
    data::DataFrame

    function Logger(net::Network{ID,T,WT,IT}) where {ID,T,WT,IT}
        data = DataFrame(:t=>T[], :output=>ID[], :event=>Symbol[], :state=>IT[])
        new(data)
    end
end

"""Log state of 'Output' object at each event"""
function (l::Logger)(t, ev, obj::Output)
    push!(getfield(l,:data), (t=t, output=obj.id, event=ev, state=obj.state_syn[]))
end

Base.getproperty(l::Logger, sym::Symbol) = getproperty(getfield(l,:data),sym)
Base.propertynames(l::Logger) = propertynames(getfield(l,:data))

################################################################################
## pretty printing                                                            ##
################################################################################

Base.show(io::IO, x::Network)  = print(io, "Network with $(length(x.inputs)) inputs and $(length(x.neurons)) neurons")
Base.show(io::IO, x::Neuron)   = print(io, "Neuron '$(x.id)' with $(length(x.next_upstream)) child-segments and $(length(x.synapses)) outgoing synapses")
Base.show(io::IO, x::Input)    = print(io, "Input '$(x.id)' with $(length(x.synapses)) outgoing synapses")
Base.show(io::IO, x::Output)   = print(io, "Output '$(x.id)'")
Base.show(io::IO, x::Segment)  = print(io, "Segment '$(x.id)' with $(length(x.next_upstream)) child-segments")
Base.show(io::IO, x::Synapse)  = print(io, "Synapse '$(x.id)' (connects to $(x.target.id))")
