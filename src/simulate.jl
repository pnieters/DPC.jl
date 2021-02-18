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
trigger!(t,T,obj::O,x,y) where {O} = @warn "No handler implemented for event-type $(T) on object-type $(O)!"
function Logger(::T) where T
    @warn "No output_logger defined for type $(T)!"
    function output_logger(args...)
    end
end

function simulate!(net, inputs::Vector{Event{T,Symbol,<:Any}}, tstop=Inf; logger! = Logger(net), show_progress=true, reset=true) where T
    event_queue = BinaryMinHeap(inputs)
    sim_clock = Ref(zero(T))

    queue!(Δt, id, obj) = push!(event_queue, Event(sim_clock[] + Δt, id, obj))

    if reset
        reset!(net)
    end

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
            break
        end

        trigger!(sim_clock[], Val(next_event.id), next_event.value, queue!, logger!)
    end

    return logger!
end