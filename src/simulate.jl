using DataStructures: BinaryMinHeap
import ProgressMeter
export simulate!, Event, is_superthreshold


##############################

function simulate!(net, inputs::Vector{<:Event{T}}, tstop=Inf; logger! = Logger(net), show_progress=true, reset=true, handle! = handle!) where {T}
    inputs = convert(Vector{Event{T}}, inputs)
    event_queue = BinaryMinHeap(inputs)
    sim_clock = Ref(zero(T))

    queue!(ev) = push!(event_queue, ev)

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
        
        handle!(next_event, queue!, logger!)
    end

    return logger!
end