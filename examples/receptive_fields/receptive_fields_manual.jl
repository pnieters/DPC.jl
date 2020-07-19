using OffsetArrays
include("utils.jl")


## Construct network
#          _________N_________           neuron's soma
#         /         |         \
#       B           B           B        dendritic branches
#     / | \       / | \       / | \ 
#  [S] [S] [S] [S] [S] [S] [S] [S] [S]   RGC ensembles (Poisson spiketrains)
#   |   |   |   |   |   |   |   |   |    
#   Y   Y   Y   Y   Y   Y   Y   Y   Y    convolution with D.O.G rfs, strided
#  /|\ /|\ /|\ /|\ /|\ /|\ /|\ /|\ /|\
# X X X X X X X X X X X X X X X X X X X  raw input pixels

rate = 200
duration = 0.01
population_size = 10
num_populations = 7
min_fire_possibility = 0.9

# 1.) linear single cell model: all input cells directly connect to the neuron, all cells must be active
# p_linear_single = pdf(volley_size_distribution(rate, duration, num_populations),5)
# -> value too low, but can't increase firing rate further -> can only increase population size

# 2.) linear population model: each cell is replaced by population, all populations connect to the neuron, enough cells must be active
# BUT: if a single neuron already fires with probability `p_single`, the number of neurons `n` that need to fire 
# out of the population of size `N` should be small enough such that the probability `p_population` > `p_single`
p_single = pdf(volley_size_distribution(rate, duration, 1),1)
max_threshold = quantile(volley_size_distribution(rate, duration,population_size), 1-p_single)
p1=plot(x->1-cdf(volley_size_distribution(rate, duration,population_size),x-1),0,population_size,label="p_population", xlabel="threshold", legend=:bottomleft)
hline!([p_single], label="p_single")
vline!([max_threshold], label="maximum threshold")
# with `num_populations` populations that means a proportionally larger population size, with each at least `max_threshold` neurons coactive for each
p2=plot(x->1-cdf(volley_size_distribution(rate, duration, num_populations*population_size),x-1),0,num_populations*population_size)
hline!([p_single], label="p_single", legend=:bottomleft)
vline!([num_populations*max_threshold], label="maximum threshold")
# Problem: one valid solution is that 2 populations don't fire at all!

# 3.) segregated dendrites model: each population projects to different segment, each segment must turn on
segment_threshold = quantile(volley_size_distribution(rate, duration,population_size), 1-min_fire_possibility)
p3 = plot(x->(1-cdf(volley_size_distribution(rate, duration,population_size),x-1))^num_populations, 0, population_size, legend=:bottomleft)
hline!([min_fire_possibility], label="minimum probability")
vline!([segment_threshold], label="maximum threshold")

stride = 2
# load input
X = VideoBuffer(times, rand(length(times), 20, 20))
# filter each frame with the RGC difference-of-gaussian kernel
Y = VideoBuffer(X.times,imfilter(X, ganglion_rf)[:, 1:stride:end,1:stride:end])

ensemble_size = 10

rgc_spiketrains = Dict{Symbol,Vector{Float64}}()
for x ∈ 1:size(Y.frames,2)
    for y ∈ 1:size(Y.frames,1)
        inputs = Y.frames[:,y,x]
        for centersurround ∈ (:center,:surround)
            for i ∈ 1:ensemble_size
                input_id = Symbol("rgc_$(centersurround)_$(y)_$(x)_$(i)")
                rgc_spiketrains[input_id] = 
            end
        end
    end
end