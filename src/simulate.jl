using DataStructures: BinaryMinHeap
import ProgressMeter
export simulate!, Event, is_superthreshold


##############################
"""simulate!(net, inputs, tstop=Inf; logger! = Logger(net), show_progress=true, reset=true, handle! = handle!)

Simulates a `Network` object `net` for all the input events in `inputs` up to maximum time `tstop` (defaults to `Inf`).

All events are handled through the `handle!` function (should be callable as `handle!(event, event_queue!, logger!)`).
Here, `event_queue` is a `BinaryMinHeap` with eltype is given by `inputs`.

An additional `logger!` (should be callable as `logger!(current_time, event_name, target_object_name, target_state)`) 
can be specified to  log all state transitions.

With `show_progress=false`, no progress-bar will be displayed during simultion.

With `reset=false`, the network state will not be reset before simulation (useful for iterative simulations).

"""
function simulate!(net, inputs::Vector{<:Event{T}}, tstop=Inf; logger! = Logger(net), show_progress=true, reset=true, handle! = handle!) where {T}
    inputs = convert(Vector{Event{T}}, inputs)
    event_queue = BinaryMinHeap(inputs)
    sim_clock = zero(T)

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
        sim_clock = next_event.t
        if sim_clock > tstop
            break
        end
        
        handle!(next_event, queue!, logger!)
    end

    return logger!
end