using ADSP, JLD2, Distributions, DataFrames, LaTeXStrings, Plots, KernelDensity


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
    
    sta_plateaus = Dict{Symbol,Vector}()
    for (segment,id) ∈ pairs(segments)
        sta_plateaus[segment] = Union{Missing,Float64}[]
    end
    sta_plateaus[:is_unique] = Bool[]

    for tt ∈ trigger_times
        plat = Dict{Symbol,Vector}()
        # collect all relevant plateaus
        for (segment,id) ∈ pairs(segments)
            plat[segment] = Vector{Union{Missing,Float64}}(filter(:t=>t->tt+plateau_window[1] ≤ t ≤ tt+plateau_window[2], groups[(event_type=:plateau_start,event_source=id)]).t .- tt)
            if isempty(plat[segment])
                push!(plat[segment], missing)
            end
        end
        plat[:is_unique] = [all(length.(values(plat)) .== 1)]

        # record all combinations of plateaus
        for comb ∈ Iterators.product(values(plat)...)
            for (i,k) ∈ enumerate(keys(plat))
                push!(sta_plateaus[k], comb[i])
            end
        end
    end

    return sta_spikes, sta_plateaus
end

"""Truncated exponential kernel with length c"""
dist_1(x,r;c=0.1) = r*exp(-r*x)/(1.0-exp(-r*c))*(0.0 ≤ x ≤ c)

"""
Convolution of 2 truncated exponentials

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

println("Loading data...")
@load "examples/rate_coding/data/examples_data.jld2" input_rates output_rates2 log_6 log_7 log_8

println("Processing data...")
kernel(t) = float(0.0 ≤ t ≤ 5e-3)
populations = (A=Symbol.("A",1:25),B=Symbol.("B",1:25),C=Symbol.("C",1:25))
segments = (A=SegmentID(NeuronID(:main),:soma),B=SegmentID(NeuronID(:main),:B),C=SegmentID(NeuronID(:main),:C))
sta6, sta_plateaus6 = extract_population_sta(log_6, populations, segments; kernel=kernel, spike_window=(-0.25,0.1), plateau_window=(-0.25,0.0))
sta7, sta_plateaus7 = extract_population_sta(log_7, populations, segments; kernel=kernel, spike_window=(-0.25,0.1), plateau_window=(-0.25,0.0))
sta8, sta_plateaus8 = extract_population_sta(log_8, populations, segments; kernel=kernel, spike_window=(-0.25,0.1), plateau_window=(-0.25,0.0))

plots = Plots.Plot[]
for (i,example,sta,sta_plateaus) ∈ [(1,:sequential, sta6, sta_plateaus6), (2,:and, sta8, sta_plateaus8), (3,:or, sta7, sta_plateaus7)]
    # # plot the spike-triggered averages of the spiking inputs
    # x = LinRange(-0.2,0.01, 251)
    # p_sta = plot(
    #     xlims=(-0.2,0.01), 
    #     legend=false, 
    #     title="$(example) segments",
    #     xlabel=L"t-t_A",
    #     ylabel=L"B(t), C(t)"
    # )
    # # plot!(x,[sta6[:C].(x) sta6[:B].(x) sta6[:A].(x)], linewidth=2, labels=["segment C" "segment B" "segment A"])
    # plot!(x,[sta[:C].(x) sta[:B].(x)], linewidth=2, labels=["segment C" "segment B"])
    # vline!([0],linewidth=2, color=:black)
    # push!(plots, p_sta)

    # plot the spike-triggered plateaus
    p_sta_plateaus = plot(
        aspect_ratio=:equal, 
        legend=false, 
        xlims=(-0.2,0.0),
        ylims=(-0.2,0.0),
        xlabel=L"\Delta t_B = t_B-t_A",
        ylabel=L"\Delta t_C = t_C-t_A",
        framestyle=:semi,
        title=latexstring("P_{$i}(\\Delta t_B,\\Delta t_C|t_A)"),
    )
    # scatter!(sta_plateaus[:B],sta_plateaus[:C], alpha=0.5, marker=:o, markersize=2, label="all pairs")
    idx = Base.:~.(ismissing.(sta_plateaus[:B]) .| ismissing.(sta_plateaus[:C]))
    K = kde((collect(skipmissing(sta_plateaus[:B][idx])),collect(skipmissing(sta_plateaus[:C][idx]))))
    contourf!(K.x, K.y, K.density',levels=25, colorbar=false)
    scatter!(sta_plateaus[:B][sta_plateaus[:is_unique]],sta_plateaus[:C][sta_plateaus[:is_unique]], color=:white, marker=:., markersize=2, label="unambiguous pairs")
    if example== :sequential
        plot!([-0.25 -0.25 -0.15 -0.1; 0 0 0 -0.1],[-0.1 -0.25 -0.25 -0.25; -0.1 0.0 -0.1 0.0], linewidth=2, color=:gray, linestyle=:dash, label="")
        plot!([-0.1, 0.0, 0.0, -0.1, -0.1], [-0.2,-0.1,0.0,-0.1,-0.2], linewidth=3, color=:white, label="")
    elseif example == :and
        plot!([-0.25 -0.1; 0 -0.1],[-0.1 -0.25; -0.1 0.0], linewidth=2, color=:gray, linestyle=:dash, label="")
        plot!([-0.1, 0.0, 0.0, -0.1, -0.1], [-0.1, -0.1, 0.0, 0.0, -0.1], linewidth=3, color=:white, label="")
    else
        plot!([-0.25 -0.1; 0 -0.1],[-0.1 -0.25; -0.1 0.0], linewidth=2, color=:gray, linestyle=:dash, label="")
        plot!([-0.2, -0.1, -0.1, 0.0, 0.0, -0.2, -0.2], [-0.1, -0.1, -0.2, -0.2, 0.0, 0.0, -0.1], linewidth=3, color=:white, label="")
    end

    push!(plots, p_sta_plateaus)
end

w = 4.7747 * 150
h = 0.33*w
# plot(plots[[1,3,2,4]]..., link=:x, layout=grid(2,2, widths=[0.5, 0.5], heights=[0.25,0.75]), size=(w,h))
plot(plots..., link=:x, layout=grid(1,3), size=(w,h))
savefig("examples/rate_coding/figures/sta.svg")
