import YAML
using DataStructures: OrderedDict

export P_superthreshold, P_fire, load_network, save_network


"""P_superthreshold(threshold, event_probabilities, event_weights)

Compute the probability with which a sum of events with given `event_weights` and `event_probabilities` 
would exceed a given `threshold`.
"""
function P_superthreshold(threshold, event_probabilities, event_weights)
    @assert length(event_probabilities) == length(event_weights)
    if isempty(event_probabilities)
        return threshold ≤ 0 ? 1.0 : 0.0
    end
    
    d = zeros(Bool, length(event_weights))
    P = zero(Float64)
    for i in 1:2^length(event_weights)
        digits!(d, i, base=2)
        if event_weights'*d ≥ threshold
            p = prod(ifelse.(d, event_probabilities, one(Float64) .- event_probabilities))
            P += p 
        end
    end
    P
end

"""P_fire(obj, volleys)

Compute the probability with which the `volleys` (provided they occur in the correct temporal sequence) 
would trigger a spike/plateau in `obj`
"""
function P_fire(obj, volleys)
    (probabilities, weights) = get(volleys, obj.id, ([],[]))
    P_syn = P_superthreshold(obj.θ_syn, probabilities, weights)
    return if isempty(obj.next_upstream) 
        P_syn
    else
        P_seg = P_superthreshold(obj.θ_seg, P_fire.(obj.next_upstream, Ref(volleys)), ones(length(obj.next_upstream)))
        P_syn*P_seg
    end
end

"""
    load_network(YAML_source)

Loads a network from a configuration file or string `YAML_source`.
"""
function load_network(YAML_file=nothing; YAML_source=nothing, id_type=Symbol, time_type=Float64, weight_type=Int, synaptic_input_type=Int) where {T}
    @assert (YAML_file === nothing) ⊻ (YAML_source === nothing) "Must provide exactly one of `YAML_file` or `YAML_source` "
    obj = if YAML_file !== nothing
        YAML.load_file(YAML_file)
    else
        YAML.load(YAML_source)
    end
    
    """ Helper function to parse optional keyword arguments, should they be specified"""
    parse_kwarg!(kwargs, obj, key, parser, key_str=String(key)) = if key_str in keys(obj) kwargs[key]=parser(pop!(obj,key_str)) end
    
    # get basic blocks
    inputs = pop!(obj, "inputs", [])
    outputs = pop!(obj, "outputs", [])
    neurons = pop!(obj, "neurons", [])
    synapses = pop!(obj, "synapses", [])

    if inputs === nothing
        inputs = []
    end
    if outputs === nothing
        outputs = []
    end
    if neurons === nothing
        neuron = []
    end
    if synapses === nothing
        synapses = []
    end

    # get default parameters
    kwargs = Dict{id_type,Any}()
    parse_kwarg!(kwargs, obj, :default_refractory_duration, time_type, "refractory_duration")
    parse_kwarg!(kwargs, obj, :default_delay, time_type, "delay")
    parse_kwarg!(kwargs, obj, :default_spike_duration, time_type, "spike_duration")
    parse_kwarg!(kwargs, obj, :default_plateau_duration, time_type, "plateau_duration")
    
    if ~isempty(obj)
        @warn "The following network parameters could not be parsed: $(join(keys(obj), ", "))"
    end

    net = Network(;
        id_type=id_type, time_type=time_type, weight_type=weight_type, synaptic_input_type=synaptic_input_type, 
        kwargs...
    )
    
    
    # internal flat store for all onject (needed for snypases)
    all_segments = Dict{id_type,Any}()

    function recurse_branches(parent, branch)
        id = id_type(pop!(branch, "id"))
        
        kwargs = Dict{Symbol,Any}()
        parse_kwarg!(kwargs, branch, :plateau_duration, time_type)
        parse_kwarg!(kwargs, branch, :θ_syn, synaptic_input_type)
        parse_kwarg!(kwargs, branch, :θ_seg, Int)

        seg = Segment(id, parent; kwargs...)
        all_segments[id] = seg

        for subbranch in pop!(branch, "branches", [])
            recurse_branches(seg, subbranch)
        end

        if ~isempty(branch)
            @warn "The following segment parameters could not be parsed: $(join(keys(branch), ", "))"
        end
    end

    for neuron in neurons
        id = id_type(pop!(neuron, "id"))
        
        kwargs = Dict{Symbol,Any}()
        parse_kwarg!(kwargs, neuron, :refractory_duration, time_type)
        parse_kwarg!(kwargs, neuron, :θ_syn, synaptic_input_type)
        parse_kwarg!(kwargs, neuron, :θ_seg, Int)

        n = Neuron(id, net; kwargs...)
        all_segments[id] = n


        for branch in pop!(neuron, "branches", [])
            recurse_branches(n, branch)
        end

        if ~isempty(neuron)
            @warn "The following neuron parameters could not be parsed: $(join(keys(neuron), ", "))"
        end
    end

    for input in inputs
        id = id_type(input["id"])
        inp = Input(id, net)
        all_segments[id] = inp
    end

    for output in outputs
        id = id_type(output["id"])
        out = Output(id, net)
        all_segments[id] = out
    end

    for synapse in synapses
        id = id_type(pop!(synapse, "id"))
        source = all_segments[id_type(pop!(synapse, "source"))]
        target = all_segments[id_type(pop!(synapse, "target"))]

        kwargs = Dict{Symbol,Any}()
        parse_kwarg!(kwargs, synapse, :delay, time_type)
        parse_kwarg!(kwargs, synapse, :weight, weight_type)
        parse_kwarg!(kwargs, synapse, :spike_duration, time_type)

        syn = Synapse(id, source, target; kwargs...)

        if ~isempty(synapse)
            @warn "The following synapse parameters could not be parsed: $(join(keys(synapse), ", "))"
        end
        all_segments[id] = syn
    end

    return (network=net, objects=all_segments)
end

"""
    save_network!(net, filename=nothing)

Saves a network to a YAML string that is returned and saved to file if `filename` is not `nothing`
"""
function save_network(net::Network{id_type, time_type, weight_type, synaptic_input_type}, filename=nothing) where {id_type, time_type, weight_type, synaptic_input_type}
    data = OrderedDict{String,Any}()

    parse_kwarg!(kwargs, obj, key, parser, key_str=String(key)) = if key_str in keys(obj) kwargs[key]=parser(pop!(obj,key_str)) end
    
    # create basic blocks
    data["inputs"] = []
    data["outputs"] = []
    data["neurons"] = []
    data["synapses"] = []
    
    function parse_synapse(obj, synapse::Synapse)
        syn = OrderedDict{String,Any}()
        syn["id"] = synapse.id
        syn["source"] = obj.id
        syn["target"] = synapse.target.id
        syn["delay"] = synapse.delay
        syn["spike_duration"] = synapse.spike_duration
        syn["weight"] = serialize(synapse.weight)
        return syn
    end

    function recurse_branches!(obj, branch)
        branch["id"] = obj.id
        branch["plateau_duration"] = obj.plateau_duration
        branch["θ_syn"] = obj.θ_syn
        branch["θ_seg"] = obj.θ_seg

        if ~isempty(obj.next_upstream)
            branch["branches"] = []
            for subobj in obj.next_upstream
                tmp = OrderedDict{String,Any}()
                recurse_branches!(subobj, tmp)
                push!(branch["branches"], tmp)
            end
        end
    end

    for neuron in net.neurons
        branch = OrderedDict{String,Any}()
        branch["id"] = neuron.id
        branch["refractory_duration"] = neuron.refractory_duration
        branch["θ_syn"] = neuron.θ_syn
        branch["θ_seg"] = neuron.θ_seg
        
        for synapse in neuron.synapses
            push!(data["synapses"], parse_synapse(neuron, synapse))
        end

        if ~isempty(neuron.next_upstream)
            branch["branches"] = []
            for subobj in neuron.next_upstream
                tmp = Dict{String,Any}()
                push!(branch["branches"], tmp)
                recurse_branches!(subobj, tmp)
            end
        end
        push!(data["neurons"],branch)
    end

    for input in net.inputs
        branch = OrderedDict{String,Any}()
        branch["id"] = input.id
        
        for synapse in input.synapses
            push!(data["synapses"], parse_synapse(input, synapse))
        end
        push!(data["inputs"],branch)
    end

    for output in net.outputs
        branch = OrderedDict{String,Any}()
        branch["id"] = output.id
        push!(data["outputs"],branch)
    end

    return if filename === nothing
        YAML.write(data)
    else
        YAML.write_file(filename, data)
    end
end
