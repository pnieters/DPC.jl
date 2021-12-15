# spike generation function for integrator experiment
using Distributions, DifferentialEquations

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

function generateContiguousVolleys(trange, population_size, λ, blob_size, background_rate=0; dt_spike=0.0, dt_volley=0)
    T = trange[2]-trange[1]
    spikes = [rand(rand(Poisson(background_rate*T))).*T.+trange[1] for i ∈ 1:population_size]
    volleys = NamedTuple{(:time,:idxs),Tuple{Float64,Vector{Int}}}[]

    smp=if isa(λ,Real) 
        rand(rand(Poisson(λ*T))).*T.+trange[1]
    else
        sample_inhomogeneous_poisson(λ, trange)
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
        idxs = mod.(rand(1:population_size).+ (0:blob_size-1), Ref(1:population_size))
        for idx ∈ idxs
            push!(spikes[idx], st)
        end

        push!(volleys, (time=st, idxs=idxs))
    end
    
    foreach(sort!, spikes)
    for j ∈ 1:population_size
        sort!(spikes[j])
        dels = Int[]
        for (k,t) ∈ enumerate(spikes[j][2:end])
            if t ≤ spikes[j][k] + dt_spike
                push!(dels, k+1)
            end
        end
        deleteat!(spikes[j], dels)
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