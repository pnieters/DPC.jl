using ADSP, CairoMakie
include("utils.jl")

files = ["bio_schema.jl", "inhibition_rf.jl", "invariant_rf.jl", "multineuron_and.jl", "multineuron_or.jl", "place_cells.jl", "probabilistic.jl", "spatio_temporal_rf.jl"]
for file in files
    let
        println("Running $(file) ...")
        include(file)
    end
end

println("DONE.")
