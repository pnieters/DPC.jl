using ADSP, JLD2, Distributions, DataFrames, Plots


function smooth(spiketrains, kernel; scale=1.0/length(spiketrains))
    function smoothed(t)
        s = 0.0
        for spiketrain in spiketrains
            for spike in spiketrain
                s += kernel(t - spike)
            end
        end
        return s*scale
    end
end

function extract_population_sta(data, populations, segments; kernel=x->pdf(Normal(0,5e-3),x), spike_window=(-0.5,0.1), plateau_window=spike_window, selector=(event_type=:spike_start,event_source=:main))
    # group event by source and event type
    data=select(data, 
        :t, 
        :event=>(col->getindex.(getproperty.(typeof.(col),:parameters),2))=>:event_type, 
        :event=>(col->getproperty.(col, :value))=>:event_source
    )
    groups = groupby(data,[:event_type, :event_source])

    # find spikes to average by
    trigger_times = groups[selector].t

    sta_spikes = Dict{Symbol,Function}()
    for (population,ids) ∈ pairs(populations)
        spike_trains = Vector{Float64}[]
        for tt ∈ trigger_times
            # filter only events in the window
            spike_train = vcat((filter(:t=>t->tt+spike_window[1] ≤ t ≤ tt+spike_window[2], groups[(event_type=:spike_start,event_source=id)]).t .- tt for id ∈ ids)...)
            push!(spike_trains, spike_train)
        end
        sta_spikes[population] = smooth(spike_trains,kernel)
    end
    
    sta_plateaus = Dict{Symbol,Vector{Float64}}()
    for (segment,id) ∈ pairs(segments)
        plateaus = Float64[]
        for tt ∈ trigger_times
            p = filter(:t=>t->tt+plateau_window[1] ≤ t ≤ tt+plateau_window[2], groups[(event_type=:plateau_start,event_source=id)]).t .- tt
            if length(p) == 1
                push!(plateaus, p[])
            elseif length(p) > 1
                println("More than one plateau for segment $(segment) in the window $([tt+plateau_window[1] , tt+plateau_window[2]]): $(p). Taking last plateau time.")
                push!(plateaus, maximum(p))
            elseif isempty(p)
                println("No plateau at all for segment $(segment) in the window $([tt+plateau_window[1] , tt+plateau_window[2]]): $(p)")
                push!(plateaus, NaN)
            end
        end
        sta_plateaus[segment] = plateaus
    end

    return sta_spikes, sta_plateaus
end

println("Loading data...")
@load "examples/rate_coding/examples_data2.jld2" input_rates output_rates output_rates2 output_rates3 output_rates4 output_rates5 spikes_1 spikes_2 spikes_3 spikes_4 spikes_5

function interpolate(inp,out,r) 
    i=searchsortedlast(inp,r)
    c = (r-inp[i])/(inp[i+1]-inp[i])
    return c*out[i+1]+(1-c)*out[i]
end

"""Truncated exponential kernel with length c"""
dist_1(x,r;c=0.1) = r*exp(-r*x)/(1.0-exp(-r*c))*(0.0 ≤ x ≤ c)

"""Convolution of two exponential kernels with lengths c and d"""
dist_2(x,r1,r2;c=0.1,d=0.1) = r1*r2/(r2-r1)*(exp(-r1*x)-exp(-r2*x))/((1.0-exp(-r1*c))*(1.0-exp(-r2*d)))*(0≤x≤c+d)

"""

Derivation:
∫₀ᵗ exp(-ax)exp(-b*(t-x))*(0≤x≤c)*(0≤t-x≤d) dx
∫ exp(-ax) exp(-b*(t-x)) dx from max(t-d,0) to min(c,t)
∫ exp(-ax) exp(-b*(t-x)) dx = ∫ exp(-ax-b*(t-x)) dx = exp(-(a-b)x-b t)/(b-a)

WLOG: d≤c
"""
function dist_2(t, a, b; c=0.1, d=0.1)
    if d>c
        return dist_2(x,b,a;d=c,c=d)
    end
    range = (min(t,c),max(0.0,t-d))
    return if range[1]≤range[2]
        0.0
    else
        a*b/((b-a)*(1-exp(-a*c))*(1-exp(-b*d)))*exp(-b*t)*(exp(-(a-b)*range[1])-exp(-(a-b)*range[2]))
    end
end

println("Processing data...")
kernel(t) = float(0.0 ≤ t ≤ 5e-3)
populations = (A=Symbol.("A",1:25),B=Symbol.("B",1:25),C=Symbol.("C",1:25))
segments = (A=SegmentID(NeuronID(:main),:soma),B=SegmentID(NeuronID(:main),:B),C=SegmentID(NeuronID(:main),:C))
sta3, sta_plateaus3 = extract_population_sta(spikes_3, populations, segments; kernel=kernel, spike_window=(-0.25,0.1), plateau_window=(-0.2,0.0))
sta4, sta_plateaus4 = extract_population_sta(spikes_4, populations, segments; kernel=kernel, spike_window=(-0.15,0.1), plateau_window=(-0.1,0.0))
sta5, sta_plateaus5 = extract_population_sta(spikes_5, populations, segments; kernel=kernel, spike_window=(-0.15,0.1), plateau_window=(-0.1,0.0))

r_A = 50.0
r_B = 35.0
r_C = 35.0
r_trigger_A = 1.0/(1.0/interpolate(input_rates, output_rates2, r_A)-0.1)
r_trigger_B = 1.0/(1.0/interpolate(input_rates, output_rates2, r_B)-0.1)
r_trigger_C = 1.0/(1.0/interpolate(input_rates, output_rates2, r_C)-0.1)

# Figure for experiment 3
num_plateau_examples = 5
plateau_idx = sample(1:length(first(values(sta_plateaus3))),num_plateau_examples;replace=false)

xticks = collect(-0.1:0.02:0.0)
xlims = (-0.2,0.01)
x = LinRange(-0.2,0.01, 251)
p3A = plot(x,sta3[:C].(x), xlims=xlims, xticks=xticks, legend=false)
plot!(p3A, x,sta3[:B].(x), xlims=xlims, xticks=xticks, legend=false)
plot!(p3A, x,sta3[:A].(x), xlims=xlims, xticks=xticks, legend=false)

p3B = plot_spike_raster(xlims, [(x->[x]).(sta_plateaus3[pop][plateau_idx]) for pop ∈ (:C,:B,:A)], 0.1, legend=false)

p3C = histogram(sta_plateaus3[:B],bins=[collect(LinRange(-0.2,-5e-3,100));0.0+eps()], normalize=true, legend=false)
histogram!(sta_plateaus3[:C],bins=[collect(LinRange(-0.2,-5e-3,100));0.0+eps()], normalize=true)
plot!(x, dist_1.(-x, r_trigger_A), linestyle=:dash, linewidth=2, color=:black)
plot!(x, dist_2.(-x, r_trigger_A, r_trigger_B), linestyle=:dash, linewidth=2, color=:black)

p3 = plot(p3A, p3B, p3C, layout=grid(3,1), size=(500,1000))
savefig("examples/rate_coding/experiment_3_STA.svg")

xlims = (-0.1,0.01)
x = LinRange(-0.1,0.01, 251)
# Figure for experiment 4
p4A = plot(x,sta4[:C].(x), xlims=xlims, xticks=xticks, legend=false)
plot!(p4A, x,sta4[:B].(x), xlims=xlims, xticks=xticks, legend=false)
plot!(p4A, x,sta4[:A].(x), xlims=xlims, xticks=xticks, legend=false)

p4B = plot_spike_raster(xlims, [(x->[x]).(sta_plateaus4[pop][plateau_idx]) for pop ∈ (:C,:B,:A)], 0.1, legend=false)

p4C = histogram(sta_plateaus4[:B],bins=[collect(LinRange(-0.1,-5e-3,100));0.0+eps()], normalize=true, legend=false)
histogram!(sta_plateaus4[:C],bins=[collect(LinRange(-0.1,-5e-3,100));0.0+eps()], normalize=true)
plot!(x, dist_1.(-x, r_trigger_A), linestyle=:dash, linewidth=2, color=:black)

p4 = plot(p4A, p4B, p4C, layout=grid(3,1), size=(500,1000))
savefig("examples/rate_coding/experiment_4_STA.svg")


# Figure for experiment 5
p5A = plot(x,sta5[:C].(x), xlims=xlims, xticks=xticks, legend=false)
plot!(p5A, x,sta5[:B].(x), xlims=xlims, xticks=xticks, legend=false)
plot!(p5A, x,sta5[:A].(x), xlims=xlims, xticks=xticks, legend=false)

p5B = plot_spike_raster(xlims, [(x->[x]).(sta_plateaus5[pop][plateau_idx]) for pop ∈ (:C,:B,:A)], 0.1, legend=false)

p5C = histogram(sta_plateaus5[:B],bins=[collect(LinRange(-0.1,-5e-3,100));0.0+eps()], normalize=true, legend=false)
histogram!(sta_plateaus5[:C],bins=[collect(LinRange(-0.1,-5e-3,100));0.0+eps()], normalize=true)
# plot!(x, dist_1.(-x, r_trigger_A), linestyle=:dash, linewidth=2, color=:black)

p5 = plot(p5A, p5B, p5C, layout=grid(3,1), size=(500,1000))
savefig("examples/rate_coding/experiment_5_STA.svg")


#=
for (sta,sta_plateaus,x) ∈ [(sta3,sta_plateaus3,LinRange(-0.245,0.05, 250)),(sta4,sta_plateaus4,LinRange(-0.145,0.05, 250)),(sta5,sta_plateaus5,LinRange(-0.145,0.05, 250))]
    plts = []
    idx = sample(1:length(first(values(sta_plateaus))),5;replace=false)
    for pop ∈ [:C,:B,:A]
        p = plot(x,sta[pop].(x))
        plot!(x, exp.(.*x).*(x≤0), linestyle=:dash)
        vline!(p, sta_plateaus[pop][idx])
        push!(plts, p)
    end
    display(plot(plts...,layout=grid(3,1)))
end

#=
for k out of n neurons:
    n*p(neuron i fires at time [t,t+dt]) * p(K≥ k neurons fired at a time ∈ [max(0,t-τ),t])
=#

#=
A1 B1 C1 D1
A B C D
A B C D
A B C D



A1 = schema A
B1 = schema C->B->A
C1 = schema (C * B)->A
D1 = schema (C + B)->A
A2 = input/output function (2 "plateau durations")

B2,3,4 = 3 spike pattern examples 
B4 = STAs
B2 = input/output function (contours)
B3 = idealized input/output function (contours + "marginals")
=#

