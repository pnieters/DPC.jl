using DataStructures: BinaryMinHeap
import ProgressMeter
export simulate!, Event, is_superthreshold


struct Event{T, ID, V}
    t::T
    id::ID
    value::V
end
Base.isless(ev1::Event, ev2::Event) = isless(ev1.t, ev2.t)

##############################
trigger!(obj::O,T,x,t) where {O} = @warn "No handler implemented for event-type $(T) on object-type $(O)!"

function simulate!(sim_state, inputs::Vector{Event{T,Symbol,<:Any}}, tstop=Inf; callback! = x->nothing, show_progress=true) where T
    event_queue = BinaryMinHeap(inputs)
    sim_clock = Ref(zero(T))

    queue!(id, Δt, obj) = push!(event_queue, Event(sim_clock[] + Δt, id, obj))

    p=nothing
    if show_progress
        p = ProgressMeter.ProgressThresh(0)
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

        trigger!(next_event.value, Val(next_event.id), queue!, sim_clock[])

        callback!(next_event)
    end

    return nothing
end