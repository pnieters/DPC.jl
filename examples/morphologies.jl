using ADSP, CairoMakie
include("utils.jl")

config_soma = """
refractory_duration: 100.01
plateau_duration: 100.0
spike_duration: 10.0
inputs:
  - id: iA1
  - id: iA2
  - id: iA3
  - id: iA4
  - id: iA5
  - id: iB1
  - id: iB2
  - id: iB3
  - id: iB4
  - id: iB5
neurons:
  - id: n
    θ_syn: 10
synapses:
  - {id: synA1, source: iA1, target: n}
  - {id: synA2, source: iA2, target: n}
  - {id: synA3, source: iA3, target: n}
  - {id: synA4, source: iA4, target: n}
  - {id: synA5, source: iA5, target: n}
  - {id: synB1, source: iB1, target: n}
  - {id: synB2, source: iB2, target: n}
  - {id: synB3, source: iB3, target: n}
  - {id: synB4, source: iB4, target: n}
  - {id: synB5, source: iB5, target: n}
"""

config_soma_segment = """
refractory_duration: 100.01
plateau_duration: 100
inputs:
  - id: iA1
  - id: iA2
  - id: iA3
  - id: iA4
  - id: iA5
  - id: iB1
  - id: iB2
  - id: iB3
  - id: iB4
  - id: iB5
neurons:
  - id: n
    θ_syn: 5
    branches:
    - id: seg
      θ_syn: 5
synapses:
  - {id: synA1, source: iA1, target: seg}
  - {id: synA2, source: iA2, target: seg}
  - {id: synA3, source: iA3, target: seg}
  - {id: synA4, source: iA4, target: seg}
  - {id: synA5, source: iA5, target: seg}
  - {id: synB1, source: iB1, target: n}
  - {id: synB2, source: iB2, target: n}
  - {id: synB3, source: iB3, target: n}
  - {id: synB4, source: iB4, target: n}
  - {id: synB5, source: iB5, target: n}
"""

config_soma_segment_r = """
refractory_duration: 100.01
plateau_duration: 100
inputs:
  - id: iA1
  - id: iA2
  - id: iA3
  - id: iA4
  - id: iA5
  - id: iB1
  - id: iB2
  - id: iB3
  - id: iB4
  - id: iB5
neurons:
  - id: n
    θ_syn: 5
    branches:
    - id: seg
      θ_syn: 5
synapses:
  - {id: synA1, source: iA1, target: n}
  - {id: synA2, source: iA2, target: n}
  - {id: synA3, source: iA3, target: n}
  - {id: synA4, source: iA4, target: n}
  - {id: synA5, source: iA5, target: n}
  - {id: synB1, source: iB1, target: seg}
  - {id: synB2, source: iB2, target: seg}
  - {id: synB3, source: iB3, target: seg}
  - {id: synB4, source: iB4, target: seg}
  - {id: synB5, source: iB5, target: seg}
"""

config_soma_parallel_segments = """
refractory_duration: 100.01
plateau_duration: 100
inputs:
  - id: iA1
  - id: iA2
  - id: iA3
  - id: iA4
  - id: iA5
  - id: iB1
  - id: iB2
  - id: iB3
  - id: iB4
  - id: iB5
neurons:
  - id: n
    θ_syn: 0
    θ_seg: 2
    branches:
    - id: seg1
      θ_syn: 5
    - id: seg2
      θ_syn: 5
synapses:
  - {id: synA1, source: iA1, target: seg1}
  - {id: synA2, source: iA2, target: seg1}
  - {id: synA3, source: iA3, target: seg1}
  - {id: synA4, source: iA4, target: seg1}
  - {id: synA5, source: iA5, target: seg1}
  - {id: synB1, source: iB1, target: seg2}
  - {id: synB2, source: iB2, target: seg2}
  - {id: synB3, source: iB3, target: seg2}
  - {id: synB4, source: iB4, target: seg2}
  - {id: synB5, source: iB5, target: seg2}
"""

config_soma_serial_segments = """
refractory_duration: 100.01
plateau_duration: 100
inputs:
  - id: iA1
  - id: iA2
  - id: iA3
  - id: iA4
  - id: iA5
  - id: iB1
  - id: iB2
  - id: iB3
  - id: iB4
  - id: iB5
neurons:
  - id: n
    θ_syn: 5
    branches:
    - id: seg2
      θ_syn: 5
      branches:
      - id: seg1
        θ_syn: 5
synapses:
  - {id: synA1, source: iA1, target: seg1}
  - {id: synA2, source: iA2, target: seg1}
  - {id: synA3, source: iA3, target: seg1}
  - {id: synA4, source: iA4, target: seg1}
  - {id: synA5, source: iA5, target: seg1}
  - {id: synB1, source: iB1, target: seg2}
  - {id: synB2, source: iB2, target: seg2}
  - {id: synB3, source: iB3, target: seg2}
  - {id: synB4, source: iB4, target: seg2}
  - {id: synB5, source: iB5, target: seg2}
  - {id: synA6, source: iA1, target: n}
  - {id: synA7, source: iA2, target: n}
  - {id: synA8, source: iA3, target: n}
  - {id: synA9, source: iA4, target: n}
  - {id: synA10, source: iA5, target: n}
"""

pop_A = [:iA1, :iA2, :iA3, :iA4, :iA5]
pop_B = [:iB1, :iB2, :iB3, :iB4, :iB5]
pop_A_volleys = [6.0,      185.0, 400.0, 550.0,  750.0, 752.0, 754.0,                   ]
pop_B_volleys = [ 14.0, 176.0,        475.0,                         775.0, 777.0, 779.0]
A_spikes = [(t+rand(), id) for t in pop_A_volleys for id in pop_A]
B_spikes = [(t+rand(), id) for t in pop_B_volleys for id in pop_B]
extra_spikes = []
all_spikes = [A_spikes; B_spikes; extra_spikes]

inputs(objects)=[Event(:input_spikes, 0.0, t, objects[id]) for (t,id) in all_spikes]

scenarios = [config_soma, config_soma_segment, config_soma_segment_r, config_soma_parallel_segments, config_soma_serial_segments]


tmax = 1000
old_ax = nothing
fig = Figure(resolution = (800, 600), show_axis=false)
for (i,scenario) in enumerate(scenarios)
    ax = fig[i, 1] = Axis(fig)
    ax_schema = fig[i, 2] = Axis(fig, width=Fixed(100), aspect=DataAspect(), backgroundcolor=:transparent, )
    local (net,objects) = load_network(YAML_source=scenario)
    local logger=simulate!(net, inputs(objects), tmax)

    # plot spikes for population A
    for (i,id) in enumerate(pop_A)
        spike_times = first.(filter(x->(x[2] == id), A_spikes))
        linesegments!(ax, 
            repeat(spike_times,inner=2), 
            -repeat(i .+ [0,0.9], outer=length(spike_times)), 
            linewidth=2,
            color=color_1
        )
    end

    # plot spikes for population B
    for (i,id) in enumerate(pop_B)
        spike_times = first.(filter(x->(x[2] == id), B_spikes))
        linesegments!(ax, 
            repeat(spike_times, inner=2), 
            -repeat(i+5 .+ [0,0.9], outer=length(spike_times)), 
            linewidth=2,
            color=color_2
        )
    end

    spike_times = filter(x->(x.object == :n && x.event == :spikes), logger.data).t
    linesegments!(ax, 
        repeat(spike_times, inner=2), 
        -repeat([0.0,11.0], outer=length(spike_times)), 
        linewidth=2,
        color=:gray50
    )

    
    if i==2
        seg = get_trace(:seg, logger.data)
        steps!(ax, [0;seg.t;tmax], 2.49 .* [0;Int.(seg.state);0] .- 6.05, color=:transparent, fill=color_1_25)
    elseif i==3
        seg = get_trace(:seg, logger.data)
        steps!(ax, [0;seg.t;tmax], 2.49 .* [0;Int.(seg.state);0] .- 11.05, color=:transparent, fill=color_2_25)
    elseif i==4
        seg1 = get_trace(:seg1, logger.data)
        seg2 = get_trace(:seg2, logger.data)
        steps!(ax, [0;seg1.t;tmax], 2.49 .* [0;Int.(seg1.state);0] .- 6.05, color=:transparent, fill=color_1_25)
        steps!(ax, [0;seg2.t;tmax], 2.49 .* [0;Int.(seg2.state);0] .- 11.05, color=:transparent, fill=color_2_25)
    elseif i==4 || i==5
        seg1 = get_trace(:seg1, logger.data)
        seg2 = get_trace(:seg2, logger.data)
        steps!(ax, [0;seg1.t;tmax], 2.49 .* [0;Int.(seg1.state);0] .- 6.05, color=:transparent, fill=color_1_25)
        steps!(ax, [0;seg2.t;tmax], 2.49 .* [0;Int.(seg2.state);0] .- 11.05, color=:transparent, fill=color_2_25)
    end

    if i==1
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>:gray))
    elseif i==2
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_2, :seg=>color_1))
    elseif i==3
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_1, :seg=>color_2))
    elseif i==4
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>:gray, :seg1=>color_2, :seg2=>color_1))
    elseif i==5
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_1, :seg1=>color_1, :seg2=>color_2))
    end

    ax.yticks = [-8.5, -3.5]
    ax.ytickformat = x->["B", "A"]

    limits!(ax_schema, Rect(-0.35,-0.25,0.7,2.5))
    hidedecorations!(ax_schema)
    hidespines!(ax_schema)

    limits!(ax, Rect(0,-11,tmax,10))
    if i<length(scenarios)
        hidexdecorations!(ax, grid=false, minorgrid=false)
    end

    if i==1
        ax.title = "Neural dynamics"
        ax_schema.title = "Morphology"
    end
end
fig

save(joinpath("figures","morphologies.pdf"), fig)
save(joinpath("figures","morphologies.svg"), fig)
save(joinpath("figures","morphologies.png"), fig)
fig