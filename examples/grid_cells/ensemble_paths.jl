using ADSP, DataStructures, Plots, ProgressMeter
include("utils.jl")

cfg = ("examples/grid_cells/cfg/ensemble_unidirectional.yaml")

## Run the experiment ##
println("Generating paths...")

function run()
    logged_segments = OrderedDict(neuron=>PropertyID(SegmentID(NeuronID(neuron), :soma), :active) for neuron ∈ [:readout1, :readout2, :readout3]) 
    plots = [plot() for i ∈ 1:3]
    counts = zeros(Int, 3)
    prog = Progress(num_paths, 0.0)

    while any(counts .< num_paths)
        # generate a single path and the response
        p,s,log = run_path(cfg, logged_segments; group_size=20, background_rate=10)

        for (i,neuron) ∈ enumerate([:readout1, :readout2, :readout3])
            # if this path ended in a somatic plateau -> generate the plot
            if counts[i] < num_paths && any(log[neuron])
                counts[i] += 1
                ProgressMeter.update!(prog, minimum(counts))
                plot!(plots[i], t->p(t)[1], t->p(t)[2], trange..., 
                    legend=false, grid=false, xticks=false, yticks=false, 
                    xlims=(domain[1][1],domain[2][1]), 
                    ylims=(domain[1][2],domain[2][2]), 
                    framestyle=:box, 
                    aspect_ratio=1, 
                    linewidth=1, 
                    color="#cccccc80", 
                    arrow=:arrow
                )
            end
        end
    end
    plots
end

plots = run()

plot(plots..., layout=grid(1,3))
savefig("examples/grid_cells/figures/ensemble_sample_paths.svg")
