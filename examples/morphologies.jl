using ADSP, CairoMakie
include("utils.jl")

set_theme!(presentation_theme)

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
  - id: iC1
  - id: iC2
  - id: iC3
  - id: iC4
  - id: iC5
  - id: iD1
  - id: iD2
  - id: iD3
  - id: iD4
  - id: iD5
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
  - id: iC1
  - id: iC2
  - id: iC3
  - id: iC4
  - id: iC5
  - id: iD1
  - id: iD2
  - id: iD3
  - id: iD4
  - id: iD5
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
  - id: iC1
  - id: iC2
  - id: iC3
  - id: iC4
  - id: iC5
  - id: iD1
  - id: iD2
  - id: iD3
  - id: iD4
  - id: iD5
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
  - id: iC1
  - id: iC2
  - id: iC3
  - id: iC4
  - id: iC5
  - id: iD1
  - id: iD2
  - id: iD3
  - id: iD4
  - id: iD5
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
  - id: iC1
  - id: iC2
  - id: iC3
  - id: iC4
  - id: iC5
  - id: iD1
  - id: iD2
  - id: iD3
  - id: iD4
  - id: iD5
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
  - {id: synA6, source: iC1, target: n}
  - {id: synA7, source: iC2, target: n}
  - {id: synA8, source: iC3, target: n}
  - {id: synA9, source: iC4, target: n}
  - {id: synA10, source: iC5, target: n}
"""

config_soma_Y = """
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
  - id: iC1
  - id: iC2
  - id: iC3
  - id: iC4
  - id: iC5
  - id: iD1
  - id: iD2
  - id: iD3
  - id: iD4
  - id: iD5
neurons:
  - id: n
    θ_syn: 5
    branches:
    - id: seg2
      θ_syn: 5
      branches:
      - id: seg1A
        θ_syn: 5
      - id: seg1B
        θ_syn: 5
synapses:
  - {id: synA1, source: iA1, target: seg1A}
  - {id: synA2, source: iA2, target: seg1A}
  - {id: synA3, source: iA3, target: seg1A}
  - {id: synA4, source: iA4, target: seg1A}
  - {id: synA5, source: iA5, target: seg1A}
  - {id: synA1, source: iB1, target: seg1B}
  - {id: synA2, source: iB2, target: seg1B}
  - {id: synA3, source: iB3, target: seg1B}
  - {id: synA4, source: iB4, target: seg1B}
  - {id: synA5, source: iB5, target: seg1B}
  - {id: synB1, source: iC1, target: seg2}
  - {id: synB2, source: iC2, target: seg2}
  - {id: synB3, source: iC3, target: seg2}
  - {id: synB4, source: iC4, target: seg2}
  - {id: synB5, source: iC5, target: seg2}
  - {id: synA6, source: iD1, target: n}
  - {id: synA7, source: iD2, target: n}
  - {id: synA8, source: iD3, target: n}
  - {id: synA9, source: iD4, target: n}
  - {id: synA10, source: iD5, target: n}
"""

pop_A = [:iA1, :iA2, :iA3, :iA4, :iA5]
pop_B = [:iB1, :iB2, :iB3, :iB4, :iB5]
pop_C = [:iC1, :iC2, :iC3, :iC4, :iC5]
pop_D = [:iD1, :iD2, :iD3, :iD4, :iD5]
pop_A_volleys = [6.0,      185.0, 400.0, 550.0,  750.0, 752.0, 754.0,                   ]
pop_B_volleys = [ 14.0, 176.0,        475.0,                         775.0, 777.0, 779.0]
pop_C_volleys = [  20.0,          190.0,    555.0,                               800]
pop_D_volleys = [1.0,              200.0,      600.0,                                   ]
A_spikes = [(t+rand(), id) for t in pop_A_volleys for id in pop_A]
B_spikes = [(t+rand(), id) for t in pop_B_volleys for id in pop_B]
C_spikes = [(t+rand(), id) for t in pop_C_volleys for id in pop_C]
D_spikes = [(t+rand(), id) for t in pop_D_volleys for id in pop_D]
extra_spikes = []
all_spikes = [A_spikes; B_spikes; C_spikes; D_spikes; extra_spikes]

inputs(objects)=[Event(:input_spikes, 0.0, t, objects[id]) for (t,id) in all_spikes]

scenarios = [config_soma, config_soma_segment, config_soma_segment_r, config_soma_parallel_segments, config_soma_serial_segments, config_soma_Y]


tmax = 800
old_ax = nothing
fig = Figure(resolution = (800, 600), show_axis=false)
for (i,scenario) in enumerate(scenarios)
    ax = fig[i, 1] = Axis(fig)
    ax_schema = fig[i, 2] = Axis(fig, width=Fixed(100), aspect=DataAspect(), backgroundcolor=:transparent, )
    local (net,objects) = load_network(YAML_source=scenario)
    local logger=simulate!(net, inputs(objects), tmax)

    spike_times = filter(x->(x.object == :n && x.event == :spikes), logger.data).t
    linesegments!(ax, 
        repeat(spike_times, inner=2), 
        repeat([-3,2], outer=length(spike_times)), 
        linewidth=2,
        color=:gray50
    )

    # plot volleys for population A
    scatter!(ax, pop_A_volleys, ones(length(pop_A_volleys)), color=color_1)
    # plot volleys for population B
    scatter!(ax, pop_B_volleys, zeros(length(pop_B_volleys)), color=color_2)

    if i==1
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>:gray))
    elseif i==2
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_2, :seg=>color_1))
    elseif i==3
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_1, :seg=>color_2))
    elseif i==4
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>:gray, :seg1=>color_2, :seg2=>color_1))
    elseif i==5
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_3, :seg1=>color_1, :seg2=>color_2))
    elseif i==6
        plot!(ax_schema, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_4, :seg1A=>color_1, :seg1B=>color_2, :seg2=>color_3))
    end

    hidedecorations!(ax_schema)
    hidespines!(ax_schema)

    if i==1
        ax.title = "Neural dynamics"
        ax_schema.title = "Topology"
    end

    if scenario == config_soma_serial_segments
        ax.yticks = [-1, 0, 1]#[-8.5, -3.5]
        ax.ytickformat = x->["C", "B", "A"]
        scatter!(ax, pop_C_volleys, -ones(length(pop_C_volleys)), color=color_3)
        limits!(ax, Rect(-100,-1.5,tmax+200,3)) #limits!(ax, Rect(0,-11,tmax,10))
        rowsize!(fig.layout, i, Fixed(75))
        hidexdecorations!(ax, grid=false, minorgrid=false)
        limits!(ax_schema, Rect(-0.35,-0.25,0.7,2.5))
    elseif scenario == config_soma_Y
        ax.yticks = [-2, -1, 0, 1]#[-8.5, -3.5]
        ax.ytickformat = x->["D", "C", "B", "A"]
        scatter!(ax, pop_C_volleys, -ones(length(pop_C_volleys)), color=color_3)
        scatter!(ax, pop_D_volleys, -2*ones(length(pop_D_volleys)), color=color_4)
        limits!(ax, Rect(-100,-2.5,tmax+200,4)) #limits!(ax, Rect(0,-11,tmax,10))
        rowsize!(fig.layout, i, Fixed(100))
        limits!(ax_schema, Rect(-0.35,-0.25,0.7,2.5))
      else
        ax.yticks = [0, 1]#[-8.5, -3.5]
        ax.ytickformat = x->["B", "A"]
        limits!(ax, Rect(-100,-0.5,tmax+200,2)) #limits!(ax, Rect(0,-11,tmax,10))
        limits!(ax_schema, Rect(-0.35,-0.25,0.7,1.5))
        hidexdecorations!(ax, grid=false, minorgrid=false)
    end
end
fig

save(joinpath("figures","morphologies.pdf"), fig)
save(joinpath("figures","morphologies.svg"), fig)
save(joinpath("figures","morphologies.png"), fig)
fig