using ADSP, Distributions, DataFrames, Plots

# Create Poisson spike-trains
population_size = 25
tspan = (0.0,500.0)

draw_spikes(population_firing_rates) = Dict(pop => [sort!(rand(rand(Poisson(firing_rate*(tspan[2]-tspan[1])))).*(tspan[2]-tspan[1]) .+ tspan[1]) for i ∈ 1:popsize] for (pop,(popsize,firing_rate)) in pairs(population_firing_rates))

function count_spikes(population_firing_rates)
    # generate spikes
    spikes = draw_spikes(population_firing_rates)
    input_spike_events = generateEventsFromSpikes(values(spikes), keys(spikes); spike_duration=net.neurons[:main].spike_duration)
    
    # simulate
    log_data=simulate!(net, input_spike_events; filter_events=ev-> ev isa Event{Float64, :spike_start} && ev.value==:main, show_progress=true)
    return size(log_data,1)
end

# experiment 1: single soma-segment
# load network
net,meta = load_yaml(Network{Float64}, "examples/rate_coding/cfg/example_1.yaml")

input_rates = LinRange(0,100,21)
output_rates = zeros(length(input_rates))
for (i,input_rate) ∈ enumerate(input_rates)
    output_rates[i] = count_spikes((A=(population_size, input_rate),)) / (tspan[2]-tspan[1])
    println("$(input_rates[i]): $(output_rates[i])")
end

# experiment 2: single plateau-segment
# load network
net,meta = load_yaml(Network{Float64}, "examples/rate_coding/cfg/example_2.yaml")

output_rates2 = zeros(length(input_rates))
for (i,input_rate) ∈ enumerate(input_rates)
    output_rates2[i] = count_spikes((A=(population_size, input_rate),)) / (tspan[2]-tspan[1])
    println("$(input_rates[i]): $(output_rates2[i])")
end
plot(input_rates, [output_rates output_rates2])


# experiment 3: three segments in a chain
# load network
net,meta = load_yaml(Network{Float64}, "examples/rate_coding/cfg/example_3.yaml")
input_rates2d = LinRange(0,50,11)
output_rates3 = zeros(length(input_rates2d),length(input_rates2d))
for (i,input_rate_B) ∈ enumerate(input_rates2d)
    for (j,input_rate_C) ∈ enumerate(input_rates2d)
        output_rates3[i,j] = count_spikes((A=(population_size, 100.0), B=(population_size, input_rate_B), C=(population_size, input_rate_C))) / (tspan[2]-tspan[1])
        println("($(input_rates2d[i]),$(input_rates2d[j])): $(output_rates3[i,j])")
    end
end
contourf(input_rates2d,input_rates2d,output_rates3, levels = 20)

# experiment 4: three segments in 2 branches, `or` connected
# load network
net,meta = load_yaml(Network{Float64}, "examples/rate_coding/cfg/example_4.yaml")
output_rates4 = zeros(length(input_rates2d),length(input_rates2d))
for (i,input_rate_B) ∈ enumerate(input_rates2d)
    for (j,input_rate_C) ∈ enumerate(input_rates2d)
        output_rates4[i,j] = count_spikes((A=(population_size, 100.0), B=(population_size, input_rate_B), C=(population_size, input_rate_C))) / (tspan[2]-tspan[1])
        println("($(input_rates2d[i]),$(input_rates2d[j])): $(output_rates4[i,j])")
    end
end
contourf(input_rates2d,input_rates2d,output_rates4, levels = 20)

# experiment 5: three segments in 2 branches, `and` connected
# load network
net,meta = load_yaml(Network{Float64}, "examples/rate_coding/cfg/example_5.yaml")
output_rates5 = zeros(length(input_rates2d),length(input_rates2d))
for (i,input_rate_B) ∈ enumerate(input_rates2d)
    for (j,input_rate_C) ∈ enumerate(input_rates2d)
        output_rates5[i,j] = count_spikes((A=(population_size, 100.0), B=(population_size, input_rate_B), C=(population_size, input_rate_C))) / (tspan[2]-tspan[1])
        println("($(input_rates2d[i]),$(input_rates2d[j])): $(output_rates5[i,j])")
    end
end
contourf(input_rates2d,input_rates2d,output_rates5, levels = 20)

@save "examples/rate_coding/examples_data.jld2" input_rates output_rates output_rates2 output_rates3 output_rates4 output_rates5
