using Distributions, Plots, YAML

export generateRFSpikes, plot_spike_raster, generateEventsFromSpikes, load_yaml, save_yaml

"""
    generateRFSpikes(trange, input, rfs, λ; group_size=ones(Int,length(rfs)), dt=1e-3, background_rate=0.0)

Generates outputs resembling the activity of multiple homogeneous neuron populations.
Emits volleys of spikes at a constant rate `λ` and encodes the receptive field responses
rfs[i](input(t)) to the input at those times in the magnitudes of the volleys.
"""
function generateRFSpikes(trange, input, rfs, λ; group_size=ones(Int,length(rfs)), dt=1e-3, background_rate=0.0)
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
        end
    end

    return spikes
end


"""
    plot_spike_raster(trange, spikes, spike_width; p=plot(), colors=1:sum(length.(spikes)), kwargs...)

Plots spiketrains belonging to several different groups given by `spikes::Vector{Vector{Vector{Float64}}}` with the given `spike_width` over a time interval `trange`.
Each group is colorcoded by the corresponding element of `colors`.
"""
function plot_spike_raster(trange, spikes, spike_width; p=plot(), colors=1:sum(length.(spikes)), kwargs...)
    plot!(p)
    k=0
    for (i,group) ∈ enumerate(spikes)
        for (j,synapse) ∈ enumerate(group)
            plot!([Shape([t,t+spike_width,t+spike_width,t],[k-0.45, k-0.45, k+0.45, k+0.45]) for t ∈ synapse], color=colors[i])
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
    load_yaml(::Type{Network}, YAML_source)

Loads a network from a configuration file or string `YAML_source`.
"""
function load_yaml(::Type{Network{T}}, YAML_file=nothing; YAML_source=nothing) where {T}
    @assert (YAML_file === nothing) ⊻ (YAML_source === nothing) "Must provide exactly one of `YAML_file` or `YAML_source` "
    obj = if YAML_file != nothing
        YAML.load_file(YAML_file)
    else
        YAML.load(YAML_source)
    end

    default_spike_duration = (get(obj,"spike_duration",1.0))
    default_plateau_duration = (get(obj,"plateau_duration",100.0))
    default_release_probability = (get(obj, "release_probability", 0.5))

    # collect all neurons belonging to the network
    neurons = Dict{Symbol,Neuron{T}}()
    for (nid,neuron) ∈ obj["neurons"]
        neuron_id = NeuronID(Symbol(nid))
        spike_duration = (get(neuron, "spike_duration", default_spike_duration))
        release_probability = (get(neuron, "release_probability", default_release_probability))

        branches = _recursive_parse_branches!(nothing, :soma, neuron["soma"])
        descendents = Dict{Symbol,Vector{Symbol}}(bid=>Symbol[] for bid ∈ keys(branches))

        # determine each branches' descendents
        for (bid,branch) ∈ branches
            if ~isa(branch["ancestor"], Nothing)
                # "I'm a direct descendent of my closest ancestor"
                push!(descendents[branch["ancestor"]], bid)
            end
        end

        # collect all segments belonging to the neuron
        segments = Dict{Symbol, Segment}()
        for (bid,branch) ∈ branches
            seg_id = SegmentID(neuron_id, Symbol(bid))
            min_synapses = branch["min_synapses"]
            min_segments = branch["min_segments"]
            plateau_duration = (get(branch, "plateau_duration", default_plateau_duration))

            next_upstream = descendents[bid]
            next_downstream = branch["ancestor"]

            synapses = Dict{Symbol,Vector{Synapse}}()
            if isa(branch["synapses"], Vector)
                for syn_id ∈ branch["synapses"]
                    synapse_id = SynapseID(seg_id, Symbol(syn_id), 1)
                    synapses[Symbol(syn_id)] = [Synapse(synapse_id, release_probability, false)]
                end
            else
                for (syn_id, syns) ∈ branch["synapses"]
                    synapses[Symbol(syn_id)] = Synapse[]

                    if isa(syns, Vector)
                        for (i,syn) ∈ enumerate(syns)
                            synapse_id = SynapseID(seg_id, Symbol(syn_id), i)
                            _load_synapse!(synapse_id, syn, synapses[Symbol(syn_id)]; default_release_probability=release_probability)
                        end
                    else
                        synapse_id = SynapseID(seg_id, Symbol(syn_id), 1)
                        _load_synapse!(synapse_id, syns, synapses[Symbol(syn_id)]; default_release_probability=release_probability)
                    end
                end
            end

            segments[bid] = Segment(seg_id, min_synapses, min_segments, false, plateau_duration, synapses, next_downstream, next_upstream)
        end

        neurons[Symbol(nid)] = Neuron(neuron_id, Float64(spike_duration), segments)
    end

    delete!(obj, "neurons")

    return (network=Network(neurons), metainfo=Dict(Symbol(key)=>value for (key,value) ∈ pairs(obj)))
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