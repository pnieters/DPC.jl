import YAML
using DataStructures: OrderedDict

export generateRFSpikes, generateSpikeVolleys, plot_spike_raster, sample_inhomogeneous_poisson, generateEventsFromSpikes, load_network, save_network

"""
    generateRFSpikes(trange, input, rfs, λ; group_size=ones(Int,length(rfs)), dt=1e-3, background_rate=0.0)

Generates outputs resembling the activity of multiple homogeneous neuron populations.
Emits volleys of spikes at a constant rate `λ` and encodes the receptive field responses
rfs[i](input(t)) to the input at those times in the magnitudes of the volleys.
"""
function generateRFSpikes(trange, input, rfs, λ; group_size=ones(Int,length(rfs)), dt=0.0, background_rate=0.0)
    T = (trange[end]-trange[1])

    spikes = Vector{Vector{Vector{Float64}}}(undef, length(rfs))
    for (i,rf) ∈ enumerate(rfs)
        N = rand(Poisson(λ*T))
        
        master_spikes = rand(N) .* T .+ trange[1]
        
        p_accept = rf.(input.(master_spikes))
        
        spikes[i] = Vector{Vector{Float64}}(undef, group_size[i])
        for j ∈ 1:group_size[i]
            background_spikes = rand(rand(Poisson(background_rate*T))) .* T .+ trange[1]
            spikes[i][j] = sort!([master_spikes[rand(length(master_spikes)) .≤ p_accept] ; background_spikes])
            dels = Int[]
            for (k,t) ∈ enumerate(spikes[i][j][2:end])
                if t ≤ spikes[i][j][k] + dt
                    push!(dels, k+1)
                end
            end
            deleteat!(spikes[i][j], dels)
        end
    end

    return spikes
end

"""
function sample_inhomogeneous_poisson(rate, trange; kwargs...)
    spiketimes = Float64[]
    rate_wrapped(u,p,t) = rate(t)
    effect!(integrator) = push!(integrator.p, integrator.t)

    prob = ODEProblem((du,u,p,t)->du.=0, Float64[],trange, spiketimes)
    jump = VariableRateJump(rate_wrapped, effect!)
    jump_prob=JumpProblem(prob, Direct(),jump)
    solve(jump_prob,Tsit5(); kwargs...)

    return spiketimes
end
"""

function generateSpikeVolleys(trange, population_sizes, λs, volley_sizes, background_rate=0; dt_spike=0.0, dt_volley=0, kwargs...)
    T = trange[2]-trange[1]
    spikes = [[rand(rand(Poisson(background_rate*T))).*T.+trange[1] for i ∈ 1:num_neurons] for num_neurons ∈ population_sizes]
    volleys = [NamedTuple{(:time,:idxs),Tuple{Float64,Vector{Int}}}[] for num_neurons ∈ population_sizes]
    for (i,(λ,v,N)) ∈ enumerate(zip(λs, volley_sizes,population_sizes))

        smp=if isa(λ,Real) 
            rand(rand(Poisson(λ*T))).*T.+trange[1]
        else
            sample_inhomogeneous_poisson(λ, trange; kwargs...)
        end

        sort!(smp)
        dels = Int[]
        for (k,t) ∈ enumerate(smp[2:end])
            if t ≤ smp[k] + dt_volley
                push!(dels, k+1)
            end
        end
        deleteat!(smp, dels)

        for st ∈ smp
            volley_size = if isa(v,Integer)
                v
            elseif isa(v,Real)
                rand(Binomial(N, v))
            else
                vv = v(st)
                if isa(vv,Integer)
                    vv
                else
                    rand(Binomial(N, vv))
                end
            end
            
            idxs = sample(1:N, volley_size; replace=false)
            for idx ∈ idxs
                push!(spikes[i][idx], st)
            end

            push!(volleys[i], (time=st, idxs=idxs))
        end
        
        foreach(sort!, spikes[i])
        for j ∈ 1:N
            sort!(spikes[i][j])
            dels = Int[]
            for (k,t) ∈ enumerate(spikes[i][j][2:end])
                if t ≤ spikes[i][j][k] + dt_spike
                    push!(dels, k+1)
                end
            end
            deleteat!(spikes[i][j], dels)
        end
    end
    return spikes, volleys
end


"""
    plot_spike_raster(trange, spikes, spike_width; p=plot(), colors=1:sum(length.(spikes)), kwargs...)

Plots spiketrains belonging to several different groups given by `spikes::Vector{Vector{Vector{Float64}}}` with the given `spike_width` over a time interval `trange`.
Each group is colorcoded by the corresponding element of `colors`.
"""
function plot_spike_raster(trange, spikes, spike_width; p=plot(), colors=1:length(spikes), highlight=[NamedTuple{(:time,:idxs),Tuple{Float64,Vector{Int}}}[] for i ∈ eachindex(spikes)], highlight_color=:red, ε=1e-3, kwargs...)
    plot!(p)
    k=0
    for (i,group) ∈ enumerate(spikes)
        
        # Set spike colors for each spike in each train of this group
        spike_colors = [convert(Vector{Any},fill(colors[i], length(spiketrain))) for spiketrain ∈ group]
        for (k,volley) ∈ enumerate(highlight[i])
            c = if isa(highlight_color, AbstractArray)
                highlight_color[i]
            else
                highlight_color
            end
            
            for j ∈ volley.idxs
                cc = if isa(c, AbstractArray)
                    c[j]
                else
                    c
                end
                
                # overwrite spike color if highlighted
                for s ∈ searchsortedfirst(group[j], volley.time-ε/2):searchsortedlast(group[j], volley.time+ε/2)
                    spike_colors[j][s] = c
                end
            end
        end

        # plot the actual spikes
        for (j,synapse) ∈ enumerate(group)
            plot!([Shape([t,t+spike_width,t+spike_width,t],[k-0.45, k-0.45, k+0.45, k+0.45]) for t ∈ synapse], color=spike_colors[j], linecolor=nothing)
            k-=1
        end
    end
    plot!(; legend=false, grid=false, yticks=false, xlim=trange, ylim=(k+1-0.5,0.5), kwargs...)
    return p
end

"""
    generateEventsFromSpikes(s, population_names; spike_duration=0.001)

Generates `Event` objects from groups of spiketrains given by `s::Vector{Vector{Vector{Float64}}}`.
Each event is named by `Symbol("\$(name)\$(sid)")`, where the `name`s are given by `population_names` and `sid` enumerates the spiketrains belonging to one group.
"""
function generateEventsFromSpikes(s, population_names; spike_duration=0.001)
    events = [
        Event{Float64,ev_type,Symbol}(ev_time, Symbol("$(groupname)$(sid)"))
        for (groupname,group) ∈ zip(population_names,s)
        for (sid,synapse) ∈ enumerate(group)
        for spike ∈ synapse
        for (ev_type,ev_time) ∈ ((:spike_start,spike),(:spike_end,spike+spike_duration))
    ]

    sort!(events)
    return events
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
