export Synapse, Segment, Neuron, Network, Input

####################################
#  Definitions of network objects  #
const default_refrac_duration = 1.0
const default_spike_duration = 5.0
const default_plateau_duration = 100.0
const default_delay = 1.0

################################################################################
## model component definitions                                                ##
################################################################################

abstract type NeuronOrSegment{ID,IT} end
abstract type AbstractComponent{ID,T,WT,IT} end

struct Segment{ID,IT} <: NeuronOrSegment{ID,IT}
    id::ID
    θ_syn::IT
    θ_seg::Int
    state_syn::Ref{IT}
    state_seg::Ref{Int}
    state::Ref{Symbol}
    
    next_downstream::NeuronOrSegment{ID,IT}
    next_upstream::Vector{Segment{ID,IT}}

    function Segment(id::ID, root::NeuronOrSegment{ID,IT}; θ_syn=1, θ_seg=1) where {ID,IT}
        this = new{ID,IT}(id, θ_syn, θ_seg, Ref(zero(IT)), Ref(0), Ref(:off), root, Segment{ID,IT}[])
        push!(root.next_upstream, this)
        return this
    end
end

struct Synapse{ID,T,WT,IT} <: AbstractComponent{ID,T,WT,IT}
    id::ID
    target::Segment{ID,IT}
    delay::T
    w::WT
    state::Ref{IT}
    on::Ref{Bool}
    function Synapse(id::ID, source, target::Segment{ID,IT}; delay::T = default_delay, w::WT=one(IT)) where {ID,T,WT,IT}
        this = new{ID,T,WT,IT}(id, target, delay, w, Ref(zero(IT)), Ref(false))
        push!(source.synapses, this)
        return this
    end
end

struct Neuron{ID,T,WT,IT} <: NeuronOrSegment{ID,IT}
    id::ID
    T_refrac::T
    θ_seg::Int
    state_seg::Ref{Int}
    state::Ref{Symbol}
    next_upstream::Vector{Segment{ID,IT}}
    synapses::Vector{Synapse{ID,T,WT,IT}}
    
    function Neuron(id::ID, net::AbstractComponent{ID,T,WT,IT}; T_refrac::T=default_refrac_duration) where {ID,T,WT,IT}
        this = new{ID,T,WT,IT}(id,T_refrac,1,Ref(0),Ref(:off),Segment{ID,IT}[],Synapse{ID,T,WT,IT}[])
        push!(net.neurons, this)
        return this
    end
end

struct Input{ID,T,WT,IT} <: AbstractComponent{ID,T,WT,IT}
    id::ID
    synapses::Vector{Synapse{ID,T,WT,IT}}

    function Input(id::ID, net::AbstractComponent{ID,T,WT,IT}) where {ID,T,WT,IT}
        this = new{ID,T,WT,IT}(id::ID,Synapse{ID,T,WT,IT}[])
        push!(net.inputs, this)
        return this
    end
end

struct Network{ID,T,WT,IT} <: AbstractComponent{ID,T,WT,IT}
    inputs::Vector{Input{ID,T,WT,IT}}
    neurons::Vector{Neuron{ID,T,WT,IT}}

    Network(id_type::Type=Symbol, time_type::Type=Float64, weight_type::Type=Int, synaptic_input_type::Type=Int) =
        new{id_type, time_type, weight_type, synaptic_input_type}(Input{id_type, time_type, weight_type, synaptic_input_type}[], Neuron{id_type, time_type, weight_type, synaptic_input_type}[])
end

################################################################################
## model behavior                                                             ##
################################################################################
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
function maybe_on!(s::Segment, queue!)::Union{Nothing,Segment}
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
        maybe_on!(s.next_downstream, queue!)

        # if we are still on (and not clamped by a downstream segment), don't forget to turn off this segment!
        if s.state[] == :on
            queue!(:plateau_end, default_plateau_duration, s)
        end
    end
    return nothing
end

"""Check if this neuron was triggered to spike"""
function maybe_on!(n::Input, queue!)
    # if the neuron was off but enough segments are on, trigger a spike
    # trigger all synapses
    for s in n.synapses
        queue!(:spike_start, s.delay, s)
    end
    return nothing
end

"""Check if this neuron was re-triggered to spike after refractoriness"""
function maybe_off!(n::Input, queue!)
    nothing
end

"""Check if this neuron was triggered to spike"""
function maybe_on!(n::Neuron, queue!)
    # if the neuron was off but enough segments are on, trigger a spike
    if n.state[] == :off  && n.state_seg[] >= n.θ_seg
        n.state[] == :refractory

        # trigger all synapses
        for s in n.synapses
            queue!(:spike_start, s.delay, s)
        end

        # don't forget to turn of :refractoriness
        queue!(:refractory_end, n.T_refrac, n)
    end
    return nothing
end

"""Check if this neuron was re-triggered to spike after refractoriness"""
function maybe_off!(n::Neuron, queue!)
    n.state[] = :off

    # check if the neuron should spike again, right away
    maybe_on!(n, queue!)
end

"""Check if this segment should switch off, since its plateau has timed out"""
function maybe_off!(s::Segment, queue!)
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

"""Check if an incoming spike should trigger an EPSP"""
function maybe_on!(s::Synapse, queue!)
    if ~s.on[]
        s.on[] = true
        # set the synapse's state for this spike
        s.state[] = s.w
        # inform the target segment about this new EPSP
        s.target.state_syn[] += s.state[]

        # don't forget to turn off EPSP
        queue!(:spike_end, default_spike_duration, s)

        # this EPSP might have triggered a plateau in the target
        maybe_on!(s.target, queue!)
    end
end

"""Turn off the EPSP"""
function maybe_off!(s::Synapse, queue!)
    # only update if the synapse isn't already off (shouldn't happen!)
    if s.on[]
        # inform the target, that the EPSP is over
        s.target.state_syn[] -= s.state[]
        # reset the synapse's state
        s.state[] = zero(s.w)
        s.on[] = false
    end
    return nothing
end
