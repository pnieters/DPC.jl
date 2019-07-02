using ADSP, DifferentialEquations, DataStructures, Plots, ProgressMeter
include("utils.jl")

## Configure the neuron and match segments to colors matching the grid ##
name = :unidirectional
cfg = ("examples/grid_cells/cfg/unidirectional.yaml")
segment_colors = Dict(
    :soma=>gridcell_colors[3],
    :d1=>gridcell_colors[1],
    :d2=>gridcell_colors[6]
)

## Run the experiment ##
println("-- Running experiment '$(name)'--")
println("Generating path (this may take a while)...")
logged_segments = OrderedDict(segment=>PropertyID(SegmentID(NeuronID(:readout), segment), :active) for segment ∈ keys(segment_colors)) 

while true
    # generate a single path and the response
    p,s,log = run_path(cfg, logged_segments; group_size=20, background_rate=10)

    print(".")
    # if this path led to somatic plateau, generate the plot
    if any(log[:soma])
        # plot the path itself
        plt1=plot(t->p(t)[1], t->p(t)[2], trange..., legend=false, grid=false, xlims=(domain[1][1],domain[2][1]), ylims=(domain[1][2],domain[2][2]), xticks=false, yticks=false, framestyle=:box, aspect_ratio=1, linewidth=3, color=:black, arrow=:arrow, title="trajectory")

        # plot each receptive field population response
        plt2=plot(legend=false,xlabel="time [s]", ylabel="rate [Hz]", title="receptive field responses")
        for (j,rf) ∈ enumerate(rfs)
            plot!(rf∘p, trange..., linewidth=2, color=gridcell_colors[j])
        end
        
        # plot the resulting spike trains
        plt3=plot_spike_raster(trange, s, 5e-3; colors=gridcell_colors, title="population spike trains")

        # plot the triggered plateaus on the dendritic segments
        plt4=plot(xlabel="time [s]", title="dendritic plateaus")
        for (i,segment) ∈ enumerate(names(log))
            if segment == :t
                continue
            end
            plot!(log[:t], log[segment].+i,   seriestype=:steppre, legend=false, yticks=false, color=segment_colors[segment], linewidth=2, xlims=trange, fillrange=i)
        end

        # compose the subplots
        display( plot(plt1,plt3,plt2,plt4, layout=grid(2,2, heights=[0.8, 0.2])))
        break
    end
end
savefig("examples/grid_cells/figures/$(name)_sample_path.svg")

println("Generating paths...")
plot()
let i=0
    prog = Progress(num_paths, 0.0)
    while i < num_paths
        # generate a single path and the response
        p,s,log = run_path(cfg, logged_segments; group_size=20, background_rate=10)

        # if this path ended in a somatic plateau -> generate the plot
        if any(log[:soma])
            i+=1
            next!(prog)
            plot!(t->p(t)[1], t->p(t)[2], trange..., legend=false, grid=false, xlims=(domain[1][1],domain[2][1]), ylims=(domain[1][2],domain[2][2]), xticks=false, yticks=false, framestyle=:box, aspect_ratio=1, linewidth=1, color="#cccccc80", arrow=:arrow)
        end
    end
end
savefig("examples/grid_cells/figures/$(name)_sample_paths.svg")
