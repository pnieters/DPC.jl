using ADSP, Distributions, DataFrames, Plots, JLD2

# Create Poisson spike-trains
draw_spikes(population_firing_rates) = Dict(pop => [sort!(rand(rand(Poisson(firing_rate*(tspan[2]-tspan[1])))).*(tspan[2]-tspan[1]) .+ tspan[1]) for i ∈ 1:popsize] for (pop,(popsize,firing_rate)) in pairs(population_firing_rates))

# simulate and count responses
function count_spikes(net, population_firing_rates)
    # generate spikes
    spikes = draw_spikes(population_firing_rates)
    input_spike_events = generateEventsFromSpikes(values(spikes), keys(spikes); spike_duration=net.neurons[:main].spike_duration)
    
    # simulate
    println("Started simulation ($(population_firing_rates) on thread $(Threads.threadid()))")
    log_data=simulate!(net, input_spike_events; filter_events=ev-> ev isa Event{Float64, :spike_start} && ev.value==:main, show_progress=false)
    return size(log_data,1)
end

# general setup
population_size = 25
tspan = (0.0,250.0)
input_rates = LinRange(0,100,21)
input_rates2d = LinRange(0,50,11)


# experiment 1: single soma-segment
println("\n\n ========== EXPERIMENT 1 ========== \n\n")
# load network
net1,meta = load_yaml(Network{Float64}, "examples/rate_coding/cfg/example_1.yaml")

output_rates = zeros(length(input_rates))
Threads.@threads for (i,input_rate) ∈ collect(enumerate(input_rates))
    output_rates[i] = count_spikes(deepcopy(net1), (A=(population_size, input_rate),)) / (tspan[2]-tspan[1])
    println("$(input_rates[i]): $(output_rates[i])")
end

# experiment 2: single plateau-segment
println("\n\n ========== EXPERIMENT 2 ========== \n\n")
# load network
net2,meta = load_yaml(Network{Float64}, "examples/rate_coding/cfg/example_2.yaml")

output_rates2 = zeros(length(input_rates))
Threads.@threads for (i,input_rate) ∈ collect(enumerate(input_rates))
    output_rates2[i] = count_spikes(deepcopy(net2), (A=(population_size, input_rate),)) / (tspan[2]-tspan[1])
    println("$(input_rates[i]): $(output_rates2[i])")
end
p12 = plot(input_rates, output_rates)
plot!(twinx(), input_rates, output_rates2)
display(p12)
savefig(p12, "examples/rate_coding/experiment_1_2.svg")


# experiment 3: three segments in a chain
println("\n\n ========== EXPERIMENT 3 ========== \n\n")
# load network
net3,meta = load_yaml(Network{Float64}, "examples/rate_coding/cfg/example_3.yaml")
output_rates3 = zeros(length(input_rates2d),length(input_rates2d))
Threads.@threads for (i,input_rate_B) ∈ collect(enumerate(input_rates2d))
    for (j,input_rate_C) ∈ enumerate(input_rates2d)
        output_rates3[i,j] = count_spikes(deepcopy(net3), (A=(population_size, 100.0), B=(population_size, input_rate_B), C=(population_size, input_rate_C))) / (tspan[2]-tspan[1])
        println("($(input_rates2d[i]),$(input_rates2d[j])): $(output_rates3[i,j])")
    end
end
p3 = contourf(input_rates2d,input_rates2d,output_rates3, levels = 20)
display(p3)
savefig(p3, "examples/rate_coding/experiment_3.svg")


# experiment 4: three segments in 2 branches, `or` connected
println("\n\n ========== EXPERIMENT 4 ========== \n\n")
# load network
net4,meta = load_yaml(Network{Float64}, "examples/rate_coding/cfg/example_4.yaml")
output_rates4 = zeros(length(input_rates2d),length(input_rates2d))
Threads.@threads for (i,input_rate_B) ∈ collect(enumerate(input_rates2d))
    for (j,input_rate_C) ∈ enumerate(input_rates2d)
        output_rates4[i,j] = count_spikes(deepcopy(net4), (A=(population_size, 100.0), B=(population_size, input_rate_B), C=(population_size, input_rate_C))) / (tspan[2]-tspan[1])
        println("($(input_rates2d[i]),$(input_rates2d[j])): $(output_rates4[i,j])")
    end
end
p4 = contourf(input_rates2d,input_rates2d,output_rates4, levels = 20)
display(p4)
savefig(p4, "examples/rate_coding/experiment_4.svg")

# experiment 5: three segments in 2 branches, `and` connected
println("\n\n ========== EXPERIMENT 5 ========== \n\n")
# load network
net5,meta = load_yaml(Network{Float64}, "examples/rate_coding/cfg/example_5.yaml")
output_rates5 = zeros(length(input_rates2d),length(input_rates2d))
Threads.@threads for (i,input_rate_B) ∈ collect(enumerate(input_rates2d))
    for (j,input_rate_C) ∈ enumerate(input_rates2d)
        output_rates5[i,j] = count_spikes(deepcopy(net5), (A=(population_size, 100.0), B=(population_size, input_rate_B), C=(population_size, input_rate_C))) / (tspan[2]-tspan[1])
        println("($(input_rates2d[i]),$(input_rates2d[j])): $(output_rates5[i,j])")
    end
end
p5 = contourf(input_rates2d,input_rates2d,output_rates5, levels = 20)
display(p5)
savefig(p5, "examples/rate_coding/experiment_5.svg")

@save "examples/rate_coding/examples_data.jld2" input_rates output_rates output_rates2 output_rates3 output_rates4 output_rates5
