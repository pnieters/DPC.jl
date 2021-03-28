using ADSP, CairoMakie
#include(joinpath(@__DIR__, "utils.jl"))
include("utils.jl")

################################################################################

fig = Figure(resolution = (1200, 900))
ax11 = fig[1, 1] = Axis(fig, title = "Single neuron, 2 branches")
ax12 = fig[1, 2] = Axis(fig, aspect=DataAspect())


config = """
refractory_duration: 1.01
inputs:
- id: i1
- id: i2
- id: i3
neurons:
- id: n
  θ_syn: 5
  θ_seg: 1
  branches: 
    - id: seg1
    - id: seg2
synapses:
- {id: syn11, source: i1, target: seg1}
- {id: syn12, source: i1, target: seg1}
- {id: syn13, source: i1, target: seg1}
- {id: syn14, source: i1, target: seg1}
- {id: syn15, source: i1, target: seg1}
- {id: syn21, source: i2, target: seg2}
- {id: syn22, source: i2, target: seg2}
- {id: syn23, source: i2, target: seg2}
- {id: syn24, source: i2, target: seg2}
- {id: syn25, source: i2, target: seg2}
- {id: syn31, source: i3, target: n}
- {id: syn32, source: i3, target: n}
- {id: syn33, source: i3, target: n}
- {id: syn34, source: i3, target: n}
- {id: syn35, source: i3, target: n}
"""

(net,objects) = load_network(YAML_source=config)

volleys = [(5.0, :i1), (10.0, :i2), (25.0, :i3), (300.0, :i1), (325.0, :i3), (455.0, :i3), (700.0, :i2), (755.0, :i3)]
input=[ Event(:input_spikes, 0.0, t+τ-0.5, objects[pop]) for (t, pop) in volleys for τ in rand(5)]

logger=simulate!(net, input)

syn11 = get_trace(:syn11, logger.data)
syn12 = get_trace(:syn12, logger.data)
syn13 = get_trace(:syn13, logger.data)
syn14 = get_trace(:syn14, logger.data)
syn15 = get_trace(:syn15, logger.data)
syn21 = get_trace(:syn21, logger.data)
syn22 = get_trace(:syn22, logger.data)
syn23 = get_trace(:syn23, logger.data)
syn24 = get_trace(:syn24, logger.data)
syn25 = get_trace(:syn25, logger.data)
syn31 = get_trace(:syn31, logger.data)
syn32 = get_trace(:syn32, logger.data)
syn33 = get_trace(:syn33, logger.data)
syn34 = get_trace(:syn34, logger.data)
syn35 = get_trace(:syn35, logger.data)
seg1  = get_trace(:seg1, logger.data)
seg2  = get_trace(:seg2, logger.data)
n     = get_trace(:n, logger.data)

spike_times = filter(x->(x.object == :n && x.event == :spikes), logger.data).t

linesegments!(ax11, repeat(spike_times, inner=2), repeat([0,16], outer=length(spike_times)), linewidth=2, color=:gray10)
for (i, syn) in enumerate([syn11,syn12,syn13,syn14,syn15])
  steps!(ax11, [0;syn.t;900], 9 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_1)
end
steps!(ax11, [0;seg1.t;900], 10 .+ 2.45 .* [0;Int.(seg1.state);0], fill=color_1_50)

for (i, syn) in enumerate([syn21,syn22,syn23,syn24,syn25])
  steps!(ax11, [0;syn.t;900], 4 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_2)
end
steps!(ax11, [0;seg2.t;900], 5 .+ 2.45 .* [0;Int.(seg2.state);0], fill=color_2_50)

for (i, syn) in enumerate([syn31,syn32,syn33,syn34,syn35])
  steps!(ax11, [0;syn.t;900], -1 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3)
end
# steps!(ax11, [0;n.t;900], 0 .+ 2.45 .* [0;Int.(n.state);0], fill=color_3_50)
# steps!(ax11, [0;seg1.t;900], 1 .+ 0.45 .* [0;Int.(seg1.state);0], fill=color_1_50)

plot!(ax12, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_3, :seg1=>color_1, :seg2=>color_2))

################################################################################

ax211 = fig[2, 1][1,1] = Axis(fig, title = "Multiple neurons")
ax212 = fig[2, 1][2,1] = Axis(fig)
ax213 = fig[2, 1][3,1] = Axis(fig)
ax22 = fig[2, 2] = Axis(fig, aspect=DataAspect())

config = """
refractory_duration: 1.01
inputs:
- id: i1
- id: i2
- id: i3
neurons:
- id: n1
  θ_syn: 5
  branches: 
    - id: seg1
      θ_syn: 5
- id: n2
  θ_syn: 5
  branches: 
    - id: seg2
      θ_syn: 5
- id: n3
  θ_syn: 1
synapses:
- {id: syn11, source: i1, target: seg1}
- {id: syn12, source: i1, target: seg1}
- {id: syn13, source: i1, target: seg1}
- {id: syn14, source: i1, target: seg1}
- {id: syn15, source: i1, target: seg1}
- {id: syn21, source: i2, target: seg2}
- {id: syn22, source: i2, target: seg2}
- {id: syn23, source: i2, target: seg2}
- {id: syn24, source: i2, target: seg2}
- {id: syn25, source: i2, target: seg2}
- {id: syn311, source: i3, target: n1}
- {id: syn321, source: i3, target: n1}
- {id: syn331, source: i3, target: n1}
- {id: syn341, source: i3, target: n1}
- {id: syn351, source: i3, target: n1}
- {id: syn312, source: i3, target: n2}
- {id: syn322, source: i3, target: n2}
- {id: syn332, source: i3, target: n2}
- {id: syn342, source: i3, target: n2}
- {id: syn352, source: i3, target: n2}
- {id: syn41, source: n1, target: n3}
- {id: syn42, source: n2, target: n3}
"""


(net,objects) = load_network(YAML_source=config)

volleys = [(5.0, :i1), (10.0, :i2), (25.0, :i3), (300.0, :i1), (325.0, :i3), (455.0, :i3), (700.0, :i2), (755.0, :i3)]
input=[ Event(:input_spikes, 0.0, t+τ-0.5, objects[pop]) for (t, pop) in volleys for τ in rand(5)]

logger=simulate!(net, input)

syn11 = get_trace(:syn11, logger.data)
syn12 = get_trace(:syn12, logger.data)
syn13 = get_trace(:syn13, logger.data)
syn14 = get_trace(:syn14, logger.data)
syn15 = get_trace(:syn15, logger.data)
syn21 = get_trace(:syn21, logger.data)
syn22 = get_trace(:syn22, logger.data)
syn23 = get_trace(:syn23, logger.data)
syn24 = get_trace(:syn24, logger.data)
syn25 = get_trace(:syn25, logger.data)
syn311 = get_trace(:syn311, logger.data)
syn321 = get_trace(:syn321, logger.data)
syn331 = get_trace(:syn331, logger.data)
syn341 = get_trace(:syn341, logger.data)
syn351 = get_trace(:syn351, logger.data)
syn312 = get_trace(:syn312, logger.data)
syn322 = get_trace(:syn322, logger.data)
syn332 = get_trace(:syn332, logger.data)
syn342 = get_trace(:syn342, logger.data)
syn352 = get_trace(:syn352, logger.data)
syn41 = get_trace(:syn41, logger.data)
syn42 = get_trace(:syn42, logger.data)
seg1  = get_trace(:seg1, logger.data)
seg2  = get_trace(:seg2, logger.data)
n1    = get_trace(:n1, logger.data)
n2    = get_trace(:n2, logger.data)
n3    = get_trace(:n3, logger.data)

spike_times1 = filter(x->(x.object == :n1 && x.event == :spikes), logger.data).t
spike_times2 = filter(x->(x.object == :n2 && x.event == :spikes), logger.data).t
spike_times3 = filter(x->(x.object == :n3 && x.event == :spikes), logger.data).t


linesegments!(ax211, repeat(spike_times1, inner=2), repeat([0,11], outer=length(spike_times1)), linewidth=2, color=:gray10)
for (i, syn) in enumerate([syn11,syn12,syn13,syn14,syn15])
  steps!(ax211, [0;syn.t;900], 4 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_1)
end
steps!(ax211, [0;seg1.t;900], 5 .+ 2.45 .* [0;Int.(seg1.state);0], fill=color_1_50)
for (i, syn) in enumerate([syn311,syn321,syn331,syn341,syn351])
  steps!(ax211, [0;syn.t;900], -1 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3)
end
# steps!(ax211, [0;n1.t;900], 10 .+ 2.45 .* [0;Int.(n1.state);0], fill=color_3_50)


linesegments!(ax212, repeat(spike_times2, inner=2), repeat([0,11], outer=length(spike_times2)), linewidth=2, color=:gray10)
for (i, syn) in enumerate([syn21,syn22,syn23,syn24,syn25])
  steps!(ax212, [0;syn.t;900], 4 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_2)
end
steps!(ax212, [0;seg2.t;900], 5 .+ 2.45 .* [0;Int.(seg2.state);0], fill=color_2_50)
for (i, syn) in enumerate([syn312,syn322,syn332,syn342,syn352])
  steps!(ax212, [0;syn.t;900], -1 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3)
end
# steps!(ax212, [0;n2.t;900], 0 .+ 2.45 .* [0;Int.(n2.state);0], fill=color_3_50)

linesegments!(ax213, repeat(spike_times3, inner=2), repeat([0,4], outer=length(spike_times3)), linewidth=2, color=:gray10)
for (i, syn) in enumerate([syn41,syn42])
  steps!(ax213, [0;syn.t;900], i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3)
end
# steps!(ax213, [0;n3.t;900], -4 .+ 0.9 .* [0;Int.(n3.state);0], fill=color_3_50)

plot!(ax22, objects[:n1], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n1=>color_3, :n2=>color_3, :n3=>color_3, :seg1=>color_1, :seg2=>color_2))
plot!(ax22, objects[:n2], root_position=Point2f0(1,0), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n1=>color_3, :n2=>color_3, :n3=>color_3, :seg1=>color_1, :seg2=>color_2))
plot!(ax22, objects[:n3], root_position=Point2f0(0.5,-1), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n3=>RGBAf0(0.2,0.2,0.2,1.0)))

################################################################################
hidedecorations!(ax12)
hidespines!(ax12)
ax12.backgroundcolor = :transparent
hidedecorations!(ax22)
hidespines!(ax22)
ax22.backgroundcolor = :transparent

colsize!(fig.layout, 2, 200)
rowsize!(fig.layout, 1, Relative(1/3))
hidexdecorations!(ax211, grid=false)
hidexdecorations!(ax212, grid=false)

save(joinpath("figures","or_branches.pdf"), fig)
save(joinpath("figures","or_branches.svg"), fig)
save(joinpath("figures","or_branches.png"), fig)
fig

##