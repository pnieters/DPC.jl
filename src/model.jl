export Synapse, Segment, Neuron, Network, PropertyID, NeuronID, SegmentID, SynapseID

####################################
# Keys to identify network objects #

# Key to retrieve a property of a network object
struct PropertyID{O}
    obj_id::O
    prop::Symbol
end

# Key to retrieve a neuron in a network
struct NeuronID
    id::Symbol
end

# Key to retrieve a neuron's dendrite segment in a network
struct SegmentID
    neuron::NeuronID
    id::Symbol
end

# Key to retrieve a dendrite segment's synapse in a network
struct SynapseID
    segment::SegmentID
    sym::Symbol
    id::Int
end
####################################


####################################
#  Definitions of network objects  #

# an individual synapse has an ID, a release probability and a state (active/inactive)
mutable struct Synapse
    id::SynapseID
    p::Float64
    active::Bool
end

# an individual segment has an ID, a synapse threshold, a segment threshold, a state (active/inactive) and incoming synaptic connections as well as information about the neighbouring segments.
mutable struct Segment
    id::SegmentID
    θ_syn::Int
    θ_seg::Int
    active::Bool
    
    synapses::Dict{Symbol,Vector{Synapse}}
    next_downstream::Union{Nothing,Symbol}
    next_upstream::Vector{Symbol}
end

# an individual neuron has an ID, stores spike and plateau durations and dendritic segments
struct Neuron{T}
    id::NeuronID
    spike_duration::T
    plateau_duration::T
    segments::Dict{Symbol, Segment}
end

# a network is a collection of neurons and allows for convenient indexing of objects by ID
struct Network{T}
    neurons::Dict{Symbol,Neuron{T}}
end
####################################


#####################################
# Methods to retrieve objects by id #
Base.getindex(n::Network, id::NeuronID) = n.neurons[id.id]
Base.getindex(n::Network, id::SegmentID) = n[id.neuron].segments[id.id]
Base.getindex(n::Network, id::SynapseID) = n[id.segment].synapses[id.sym][id.id]
Base.getindex(n::Network, id::PropertyID) = getproperty(n[id.obj_id], id.prop)
#######################################