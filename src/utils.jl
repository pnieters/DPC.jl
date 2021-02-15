using Distributions, Plots, YAML, DifferentialEquations

export generateRFSpikes, generateSpikeVolleys, plot_spike_raster, sample_inhomogeneous_poisson, generateEventsFromSpikes, load_network

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


function _recursive_parse_branches!(ancestor_id, branch_id, branch_obj)
    sub_branches = haskey(branch_obj, "branches") && ~isa(branch_obj["branches"],Nothing) && ~isempty(branch_obj["branches"]) ?
        (_recursive_parse_branches!(branch_id, Symbol(subbranch_id), subbranch_obj) for (subbranch_id, subbranch_obj) ∈ branch_obj["branches"]) : []

    delete!(branch_obj, "branches")

    branch_obj["ancestor"] = ancestor_id
    branches = merge!(Dict(branch_id=>branch_obj), sub_branches...)
    return branches
end

function _load_synapse!(syn_id, syn, synapses; default_release_probability=0.5)
    synapse = if isa(syn, Dict)
        transmitter = Symbol(get(syn, "transmitter", "excitatory"))#get(Dict("excitatory"=>:excitatory,"inhibitory"=>:inhibitory,"reward"=>:reward,"punish"=>:punish), lowercase(get(syn, "transmitter", "excitatory")), :excitatory)
        release_probability = (get(syn, "release_probability", default_release_probability))
        Synapse(syn_id, Float64(release_probability), false)
    else
        Synapse(syn_id, Float64(syn), false)
    end
    push!(synapses, synapse)
end


"""
    load_network(YAML_source)

Loads a network from a configuration file or string `YAML_source`.
"""
function load_network(YAML_file=nothing; YAML_source=nothing) where {T}
    @assert (YAML_file === nothing) ⊻ (YAML_source === nothing) "Must provide exactly one of `YAML_file` or `YAML_source` "
    obj = if YAML_file !== nothing
        YAML.load_file(YAML_file)
    else
        YAML.load(YAML_source)
    end

    
    # get basic blocks
    inputs = pop!(obj, "inputs", [])
    neurons = pop!(obj, "neurons", [])
    synapses = pop!(obj, "synapses", [])

    # get network type params
    id_type = getproperty(Base, Symbol(pop!(obj, "id_type", "Symbol")))
    time_type = getproperty(Base, Symbol(pop!(obj, "time_type", "Float64")))
    weight_type = getproperty(Base, Symbol(pop!(obj, "weight_type", "Int")))
    synaptic_input_type = getproperty(Base, Symbol(pop!(obj, "synaptic_input_type", "Int")))
    
    # internal flat store for all onject (needed for snypases)
    all_segments = Dict{id_type,Any}()

    net = Network(id_type, time_type, weight_type, synaptic_input_type)

    parse_kwarg!(kwargs, obj, key, parser) = if key in keys(obj) kwargs[Symbol(key)]=parser(obj[key]) end

    function recurse_branches(parent, branch)
        id = id_type(branch["id"])
        
        kwargs = Dict{Symbol,Any}()
        parse_kwarg!(kwargs, branch, "θ_syn", synaptic_input_type)
        parse_kwarg!(kwargs, branch, "θ_seg", Int)

        seg = Segment(id, parent; kwargs...)
        all_segments[id] = seg

        for subbranch in get(branch, "branches", [])
            recurse_branches(seg, subbranch)
        end
    end

    for neuron in neurons
        id = id_type(neuron["id"])
        
        kwargs = Dict{Symbol,Any}()
        parse_kwarg!(kwargs, neuron, "T_refrac", time_type)

        n = Neuron(id, net; kwargs...)
        all_segments[id] = n

        for branch in get(neuron, "branches", [])
            recurse_branches(n, branch)
        end
    end

    for input in inputs
        id = id_type(input["id"])
        inp = Input(id, net)
        all_segments[id] = inp
    end

    for synapse in synapses
        id = id_type(synapse["id"])
        source = all_segments[id_type(synapse["source"])]
        target = all_segments[id_type(synapse["target"])]

        kwargs = Dict{Symbol,Any}()
        parse_kwarg!(kwargs, synapse, "delay", time_type)
        parse_kwarg!(kwargs, synapse, "w", weight_type)

        syn = Synapse(id, source, target; kwargs...)
        all_segments[id] = syn
    end

    return (network=net, objects=all_segments)
end

"""
    save_yaml(::Type{Network}, filename=nothing)

Saves a network to a configuration file or string, if no filename is given.
"""
function save_yaml(net::Network, filename=nothing, metainfo=Dict{Symbol,Any}())

    neurons = Dict{Symbol,Dict}()
    for (nid,neuron) ∈ net.neurons
        
        # parse all segments
        segments = Dict{Symbol,Dict}()
        for (segid,segment) ∈ neuron.segments

            synapse_groups = Dict{Symbol,Vector}()
            for (synid,synapse_group) ∈ segment.synapses
                synapses = Vector{Dict}()
                for synapse ∈ synapse_group
                    push!(synapses, Dict(
                        :release_probability => synapse.p
                    ))
                end
                synapse_groups[synid] = synapses
            end

            segments[segid] = Dict(
                :min_synapses => segment.θ_syn,
                :min_segments => segment.θ_seg,
                :plateau_duration => segment.plateau_duration,
                :synapses => synapse_groups,
                :branches => Dict{Symbol,Dict}()
            )
        end

        function _nest_recursive(segid, segments)
            seg = segments[segid]
            seg[:branches] = Dict(bid => _nest_recursive(bid, segments) for bid ∈ neuron.segments[segid].next_upstream)
            return seg
        end

        neurons[nid] = Dict(
            :spike_duration => neuron.spike_duration, 
            :soma => _nest_recursive(:soma, segments)
        )
    end

    data = deepcopy(metainfo)
    data[:neurons] = neurons

    return if filename === nothing
        YAML.write(data)
    else
        YAML.write_file(filename, data)
    end
end