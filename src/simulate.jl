using DataFrames, ProgressMeter
using DataStructures: BinaryMinHeap

export simulate!, Event, is_superthreshold


struct Event{T, ID, V}
    t::T
    id::ID
    value::V
end
Base.isless(ev1::Event, ev2::Event) = isless(ev1.t, ev2.t)

##############################
# Handle events:


#TODO simplify
function handle!(ev::Event{T,Symbol,V}, net, queue!) where {T,V}
    if ev.id == :spike_start
        #Handle the rising edge of a spike 
        maybe_on!(ev.value, queue!)
    elseif ev.id == :spike_end
        #Handle the falling edge of a spike
        maybe_off!(ev.value, queue!)
    elseif ev.id == :plateau_end 
        #Handle the falling edge of a plateau
        maybe_off!(ev.value, queue!)
    elseif ev.id == :refractory_end
        #Handle the end of a refractory period
        maybe_off!(ev.value, queue!)
    end
    nothing
end


##############################
"""
Dict{Symbol,Any}(), filter_events=ev->true, 

# attach logging callback to signals tagged for logging
logged_signal_names = keys(log)
get_logged_signals(sim_state) = [sim_state[prop] for prop ∈ values(log)]

log_data = DataFrame(:t=>Float64[], :event=>Event[], Base.:(=>).(logged_signal_names, [typeof(v)[] for v ∈ get_logged_signals(sim_state)])...)
"""

function simulate!(sim_state, inputs::Vector{Event{T,Symbol,<:Any}}, tstop=Inf; callback! = x->nothing, show_progress=true) where T
    event_queue = BinaryMinHeap(inputs)
    sim_clock = Ref(zero(T))

    queue!(id, Δt, obj) = push!(event_queue, Event(sim_clock[] + Δt, id, obj))

    p=nothing
    if show_progress
        p = ProgressThresh(0)
    end
    
    while ~isempty(event_queue)
        if show_progress
            ProgressMeter.update!(p, length(event_queue))
        end

        # Handle next event
        next_event = pop!(event_queue)
        sim_clock[] = next_event.t
        if sim_clock[] > tstop
            return nothing
        end

        handle!(next_event, sim_state, queue!)

        callback!(next_event)
    end

    return nothing
end