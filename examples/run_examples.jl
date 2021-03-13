using ADSP, CairoMakie
include("utils.jl")

files = ["coincidence.jl", "latch.jl", "morphologies.jl", "probabilistic.jl", "spatio_temporal_rf.jl", "stochastic_latch.jl"]
for file in files
    let
        println("Running $(file) ...")
        include(file)
    end
end

println("DONE.")
