using DataFrames, ProgressMeter
using DataStructures: BinaryMinHeap

export simulate!, Event, is_superthreshold

struct Event{T, K, V}
    t::T
    value::V
end
Base.isless(ev1::Event, ev2::Event) = isless(ev1.t, ev2.t)

##############################
# Handle events:
#   0. fallback: do nothing
(ev::Event)(context...) = nothing

#   1. spikes
"""
Handle the rising edge of a spike by stochastically setting the corresponding signal of all receiving synapses to `true` according to probability `p`.
"""
function (ev::Event{T,:spike_start,Symbol})(net, queue) where T
    for (n_id,neuron) ∈ net.neurons
        needs_update = false
        for (seg_id,segment) ∈ neuron.segments
            if ev.value ∈ keys(segment.synapses)
                for synapse ∈ segment.synapses[ev.value]
                    if rand() < synapse.p
                        needs_update = true
                        synapse.active = true
                    end
                end
            end
        end
        
        if needs_update
            update!(neuron, ev.t, queue)
        end
    end
    nothing
end

"""
Handle the falling edge of a spike by setting the corresponding signal of all receiving synapses to `false`.
"""
function (ev::Event{T,:spike_end,Symbol})(net, queue) where T
    for (n_id,neuron) ∈ net.neurons
        needs_update = false
        for (seg_id,segment) ∈ neuron.segments
            if ev.value ∈ keys(segment.synapses)
                needs_update = true
                for synapse ∈ segment.synapses[ev.value]
                    synapse.active = false
                end
            end
        end
        
        if needs_update
            update!(neuron, ev.t, queue)
        end
    end
    nothing
end

#   2. plateaus
"""
Handle the rising edge of a plateau by setting the corresponding signal to `true`.
"""
function (ev::Event{T,:plateau_start,SegmentID})(net, queue) where T
    segment = net[ev.value]
    segment.active = true

    neuron = net[ev.value.neuron]
    update!(neuron, ev.t, queue)
    nothing
end

"""
Handle the falling edge of a plateau by setting the corresponding signal to `false`.
"""
function (ev::Event{T,:plateau_end,SegmentID})(net, queue) where T
    segment = net[ev.value]
    segment.active = false

    neuron = net[ev.value.neuron]
    update!(neuron, ev.t, queue)
    nothing
end
##############################


"""
Check if the neuron's segment is superthreshold and could thus emit a plateau potential now.
"""
function is_superthreshold(neuron::Neuron, s::Segment)
    Σ_synapses = 0
    
    # add synaptic contributions
    for (sym, synapses) ∈ s.synapses
        for syn ∈ synapses
            Σ_synapses += syn.active
        end
    end
    
    Σ_segments = isa(s.next_downstream, Nothing) ? 0 : neuron.segments[s.next_downstream].active
    # add contributions from upstream branches
    for seg_id ∈ s.next_upstream
        Σ_segments += neuron.segments[seg_id].active
    end
    
    Σ_synapses ≥ s.θ_syn && Σ_segments ≥ s.θ_seg
end

"""
Update the neuron's state, potentially triggering a plateau potential.
"""
function update!(neuron::Neuron{T}, t, queue) where {T}
    for (seg_id, segment) ∈ neuron.segments
        if ~segment.active && is_superthreshold(neuron, segment)
            ev = Event{T, :plateau_start, SegmentID}(t, segment.id)
            if !(ev in queue.valtree) # guard against repeated events
                push!(queue, ev)
                push!(queue, Event{T,:plateau_end,SegmentID}(t+neuron.plateau_duration, segment.id))
            end
            
            if seg_id == :soma
                ev = Event{T, :spike_start, Symbol}(t, neuron.id.id)
                if !(ev in queue.valtree) # guard against repeated events
                    push!(queue, Event{T,:spike_start,Symbol}(t, neuron.id.id))
                    push!(queue, Event{T,:spike_end,Symbol}(t+neuron.spike_duration, neuron.id.id))
                end
            end
        end
    end
end



function simulate!(sim_state, inputs; log=Dict{Symbol,Any}(), filter_events=ev->true, show_progress=true, T=Float64) where V
    event_queue = BinaryMinHeap(convert(Vector{Event{T}},inputs))

    # attach logging callback to signals tagged for logging
    logged_signal_names = keys(log)
    get_logged_signals(sim_state) = [sim_state[prop] for prop ∈ values(log)]
    
    log_data = DataFrame(:t=>Float64[], :event=>Event[], Base.:(=>).(logged_signal_names, [typeof(v)[] for v ∈ get_logged_signals(sim_state)])...)
    
    p=nothing
    if show_progress
        p = ProgressThresh(0)
    end
    
    while ~isempty(event_queue)
        if show_progress
            ProgressMeter.update!(p, length(event_queue))
        end
        next_event! = pop!(event_queue)
        next_event!(sim_state, event_queue)
        
        if filter_events(next_event!)
            push!(log_data, [next_event!.t, next_event!, get_logged_signals(sim_state)...])
        end
    end

    return log_data
end