using ADSP, Distributions, Plots

gr(titlefontsize=12, legendfontsize=10, labelfontsize=10)
const mypalette = (orange="#fcaf3e", green="#8ae234", purple="#ad7fa8", turquise="#729f7e", yellow="#fce94f", blue="#729fcf", red="#ef2929", gray="#cccccc")

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
