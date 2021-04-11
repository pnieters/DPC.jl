using ADSP, CairoMakie
#include(joinpath(@__DIR__, "utils.jl"))
include("utils.jl")


fig = Figure(resolution = (800, 800))
volleys = [(45.0, :i1), (75.0, :i2), (110.0, :i3), (250.0, :i1), (300.0, :i3), (450.0, :i3), (500.0, :i2), (550.0, :i3)]

################################################################################
config = """
refractory_duration: 5.01
inputs:
- id: i11
- id: i12
- id: i13
- id: i14
- id: i15
- id: i21
- id: i22
- id: i23
- id: i24
- id: i25
- id: i31
- id: i32
- id: i33
- id: i34
- id: i35
neurons:
- id: n
  θ_syn: 5
  θ_seg: 2
  branches: 
    - id: seg1
      θ_syn: 5
    - id: seg2
      θ_syn: 5
synapses:
- {id: syn11, source: i11, target: seg1}
- {id: syn12, source: i12, target: seg1}
- {id: syn13, source: i13, target: seg1}
- {id: syn14, source: i14, target: seg1}
- {id: syn15, source: i15, target: seg1}
- {id: syn21, source: i21, target: seg2}
- {id: syn22, source: i22, target: seg2}
- {id: syn23, source: i23, target: seg2}
- {id: syn24, source: i24, target: seg2}
- {id: syn25, source: i25, target: seg2}
- {id: syn31, source: i31, target: n}
- {id: syn32, source: i32, target: n}
- {id: syn33, source: i33, target: n}
- {id: syn34, source: i34, target: n}
- {id: syn35, source: i35, target: n}
"""

(net,objects) = load_network(YAML_source=config)
input=[ Event(:input_spikes, 0.0, t+τ-0.5, objects[Symbol("$(pop)$(i)")]) for (t, pop) in volleys for (i,τ) in enumerate(5*rand(5))]


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

plateau1_starts = filter(x->(x.object == :seg1 && x.event == :plateau_starts), logger.data).t
plateau2_starts = filter(x->(x.object == :seg2 && x.event == :plateau_starts), logger.data).t
spike_times = filter(x->(x.object == :n && x.event == :spikes), logger.data).t



xticks,xtickformat = make_manual_ticks([0;200;400;600;first.(volleys) .+ 5], ["0ms";"200ms";"400ms";"600ms";["τ$(i)" for i ∈ ['₁','₂','₃','₄','₅','₆','₇','₈']]])
yticks,ytickformat = make_manual_ticks(0.5:14.5, reverse!(["$(grp)$(sub)" for grp in ["A","B","C"] for sub in ['₁','₂','₃','₄','₅']]))

ax11 = fig[1, 1] = Axis(fig, title = "Single neuron, 2 branches (AND-configuration)"; xticks, xtickformat,  yticks, ytickformat)
ax12 = fig[1, 2] = Axis(fig, aspect=DataAspect())



# linesegments!(ax11, repeat(spike_times, inner=2), repeat([0,15], outer=length(spike_times)), linewidth=2, linestyle=:dash, color=:gray10)
for (i, syn) in enumerate([syn11,syn12,syn13,syn14,syn15])
  steps!(ax11, [0;syn.t;650], 9 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_1, color=:transparent)
end
steps!(ax11, [0;seg1.t;650], 10 .+ 2.45 .* [0;Int.(seg1.state);0], fill=color_1_50, color=:transparent)

for (i, syn) in enumerate([syn21,syn22,syn23,syn24,syn25])
  steps!(ax11, [0;syn.t;650], 4 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_2, color=:transparent)
end
steps!(ax11, [0;seg2.t;650], 5 .+ 2.45 .* [0;Int.(seg2.state);0], fill=color_2_50, color=:transparent)

for (i, syn) in enumerate([syn31,syn32,syn33,syn34,syn35])
  steps!(ax11, [0;syn.t;650], -1 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3, color=:transparent)
end


lines!.(ax11, Rect.(plateau1_starts .- 5.5, 9.95, 11, 5.1), color=:gray10, linewidth=2)
lines!.(ax11, Rect.(plateau2_starts .- 5.5, 4.95, 11, 5.1), color=:gray10, linewidth=2)
lines!.(ax11, Rect.(spike_times .- 5.5, -0.05, 11, 5.1), color=:gray10, linewidth=2)

pn=plot!(ax12, objects[:n], ports=Dict(:n=>[:C],:seg1=>[:A],:seg2=>[:B]), angle_between=30/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_3, :seg1=>color_1, :seg2=>color_2))
ports = Dict(pn.attributes[:ports][])

text!(ax12, " A ", position=ports[:A] - Point2f0(-0.3,0), align=(:left, :center), textsize=0.25, color=:black)
text!(ax12, " B ", position=ports[:B] - Point2f0( 0.3,0), align=(:right, :center), textsize=0.25, color=:black)
text!(ax12, " C ", position=ports[:C] - Point2f0(-0.3,0), align=(:left, :center), textsize=0.25, color=:black)
arrows!(ax12, [ports[:A]-Point2f0(-0.3,0), ports[:B]-Point2f0(0.3,0), ports[:C]-Point2f0(-0.3,0)] , [Point2f0(-0.2,0.0),Point2f0(0.2,0.0),Point2f0(-0.2,0.0)], linewidth=2, arrowsize = [-20, 20, -20], color=:black, arrowcolor=:black)


xlims!(ax12, [-1,1])
################################################################################


config = """
refractory_duration: 5.01
inputs:
- id: i11
- id: i12
- id: i13
- id: i14
- id: i15
- id: i21
- id: i22
- id: i23
- id: i24
- id: i25
- id: i31
- id: i32
- id: i33
- id: i34
- id: i35
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
  θ_syn: 2
synapses:
- {id: syn11, source: i11, target: seg1}
- {id: syn12, source: i12, target: seg1}
- {id: syn13, source: i13, target: seg1}
- {id: syn14, source: i14, target: seg1}
- {id: syn15, source: i15, target: seg1}
- {id: syn21, source: i21, target: seg2}
- {id: syn22, source: i22, target: seg2}
- {id: syn23, source: i23, target: seg2}
- {id: syn24, source: i24, target: seg2}
- {id: syn25, source: i25, target: seg2}
- {id: syn311, source: i31, target: n1}
- {id: syn321, source: i32, target: n1}
- {id: syn331, source: i33, target: n1}
- {id: syn341, source: i34, target: n1}
- {id: syn351, source: i35, target: n1}
- {id: syn312, source: i31, target: n2}
- {id: syn322, source: i32, target: n2}
- {id: syn332, source: i33, target: n2}
- {id: syn342, source: i34, target: n2}
- {id: syn352, source: i35, target: n2}
- {id: syn41, source: n1, target: n3}
- {id: syn42, source: n2, target: n3}
"""


(net,objects) = load_network(YAML_source=config)
input=[ Event(:input_spikes, 0.0, t+τ-0.5, objects[Symbol("$(pop)$(i)")]) for (t, pop) in volleys for (i,τ) in enumerate(5*rand(5))]

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

plateau1_starts = filter(x->(x.object == :seg1 && x.event == :plateau_starts), logger.data).t
plateau2_starts = filter(x->(x.object == :seg2 && x.event == :plateau_starts), logger.data).t
spike_times1 = filter(x->(x.object == :n1 && x.event == :spikes), logger.data).t
spike_times2 = filter(x->(x.object == :n2 && x.event == :spikes), logger.data).t
spike_times3 = filter(x->(x.object == :n3 && x.event == :spikes), logger.data).t

yticks1,ytickformat1 = make_manual_ticks(0.5:9.5, reverse!([["A$(sub)" for sub in ['₁','₂','₃','₄','₅']];["C$(sub)" for sub in ['₁','₂','₃','₄','₅']]]))
yticks2,ytickformat2 = make_manual_ticks(0.5:9.5, reverse!([["B$(sub)" for sub in ['₁','₂','₃','₄','₅']];["C$(sub)" for sub in ['₁','₂','₃','₄','₅']]]))
yticks3,ytickformat3 = make_manual_ticks(0.5:1.5, reverse!(["N₁";"N₂"]))

lower_grid = fig[2, 1] = GridLayout()
ax211 = lower_grid[1,1] = Axis(fig, title = "Multiple neurons (AND-configuration)"; yticks=yticks1, ytickformat=ytickformat1)
ax212 = lower_grid[2,1] = Axis(fig; yticks=yticks2, ytickformat=ytickformat2)
ax213 = lower_grid[3,1] = Axis(fig; yticks=yticks3, ytickformat=ytickformat3)
ax22 = fig[2, 2] = Axis(fig, aspect=DataAspect())

steps!(ax211, [0;seg1.t;650], 5 .+ 2.45 .* [0;Int.(seg1.state);0], fill=color_1_50, color=:transparent)
# linesegments!(ax211, repeat(spike_times1, inner=2), repeat([0,10], outer=length(spike_times1)), linestyle=:dash, linewidth=2, color=:gray10)
for (i, syn) in enumerate([syn11,syn12,syn13,syn14,syn15])
  steps!(ax211, [0;syn.t;650], 4 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_1, color=:transparent)
end
for (i, syn) in enumerate([syn311,syn321,syn331,syn341,syn351])
  steps!(ax211, [0;syn.t;650], -1 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3, color=:transparent)
end

steps!(ax212, [0;seg2.t;650], 5 .+ 2.45 .* [0;Int.(seg2.state);0], fill=color_2_50, color=:transparent)
# linesegments!(ax212, repeat(spike_times2, inner=2), repeat([0,10], outer=length(spike_times2)), linestyle=:dash, linewidth=2, color=:gray10)
for (i, syn) in enumerate([syn21,syn22,syn23,syn24,syn25])
  steps!(ax212, [0;syn.t;650], 4 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_2, color=:transparent)
end
for (i, syn) in enumerate([syn312,syn322,syn332,syn342,syn352])
  steps!(ax212, [0;syn.t;650], -1 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3, color=:transparent)
end

for (i, syn) in enumerate([syn41,syn42])
  steps!(ax213, [0;syn.t;650], -1.0 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=:gray50, color=:transparent)
end

lines!.(ax211, Rect.(spike_times1 .- 5.5, -0.05, 11, 5.1), color=:gray10, linewidth=2)
lines!.(ax211, Rect.(plateau1_starts .- 5.5,  4.95, 11, 5.1), color=:gray10, linewidth=2)
lines!.(ax212, Rect.(spike_times2 .- 5.5, -0.05, 11, 5.1), color=:gray10, linewidth=2)
lines!.(ax212, Rect.(plateau2_starts .- 5.5,  4.95, 11, 5.1), color=:gray10, linewidth=2)
lines!.(ax213, Rect.(spike_times3 .- 5.5, -0.05, 11, 2.1), color=:gray10, linewidth=2)

pn=plot!(ax22, objects[:n3], ports=Dict(:n3=>[:dummy, :A, :B]), root_position=Point2f0(0.5,-3.5), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n3=>:gray50))
# plot synapses
ports = Dict(pn.attributes[:ports][])
# arrows!(ax22, [ports[x] - Point2f0(0.4, 0) for x in [:C,:B,:A]] , fill(Point2f0(0.3,0),3), linewidth=2, arrowsize = 20)
arrows!(ax22, [ports[:A]-Point2f0(0.25,0), ports[:B]-Point2f0(-0.25,0)] , [Point2f0(0.1,0.0),Point2f0(-0.1,0.0)], linewidth=2, arrowsize = [20, -20], color=:black, arrowcolor=:black)
lines!(ax22, [Point2f0(0.0, 0.0), ports[:A] - Point2f0(0.5,-0.25), ports[:A] - Point2f0(0.25,0.0)], color=:black, linewidth=2)
lines!(ax22, [Point2f0(1.0, -2), ports[:B] - Point2f0(-0.5,-0.25), ports[:B] - Point2f0(-0.25,0.0)], color=:black, linewidth=2)

pn2=plot!(ax22, objects[:n1], ports=Dict(:n1=>[:C],:seg1=>[:A]), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n1=>color_3, :seg1=>color_1), linewidth=2)
pn3=plot!(ax22, objects[:n2], ports=Dict(:n2=>[:C],:seg2=>[:B]), root_position=Point2f0(1,-2.5), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n2=>color_3, :seg2=>color_2), linewidth=2)
ports2 = Dict(pn2.attributes[:ports][])
ports3 = Dict(pn3.attributes[:ports][])
arrows!(ax22, [ports2[:A]-Point2f0(-0.3,0), ports3[:B]-Point2f0(0.3,0), ports2[:C]-Point2f0(-0.3,0), ports3[:C]-Point2f0(0.3,0)] , [Point2f0(-0.2,0.0),Point2f0(0.2,0.0),Point2f0(-0.2,0.0),Point2f0(0.2,0.0)], linewidth=2, arrowsize = [-20, 20, -20, 20], color=:black, arrowcolor=:black)

text!(ax22, " A ", position=ports2[:A] - Point2f0(-0.3,0), align=(:left, :center), textsize=0.25, color=:black)
text!(ax22, " B ", position=ports3[:B] - Point2f0( 0.3,0), align=(:right, :center), textsize=0.25, color=:black)
text!(ax22, " C ", position=ports2[:C] - Point2f0(-0.3,0), align=(:left, :center), textsize=0.25, color=:black)
text!(ax22, " C ", position=ports3[:C] - Point2f0( 0.3,0), align=(:right, :center), textsize=0.25, color=:black)



ylims!(ax213, [-0.5, 2.5])
ylims!(ax22, [-3.75, 1.5])
################################################################################

hidedecorations!(ax12)
hidespines!(ax12)
ax12.backgroundcolor = :transparent
hidedecorations!(ax22)
hidespines!(ax22)
ax22.backgroundcolor = :transparent

colsize!(fig.layout, 2, Fixed(125))
rowsize!(fig.layout, 1, Relative(15/37))
rowsize!(lower_grid, 3, Relative(2/22))
hidexdecorations!(ax211, grid=false)
hidexdecorations!(ax212, grid=false)

save(joinpath("figures","multineuron_and.pdf"), fig)
save(joinpath("figures","multineuron_and.svg"), fig)
save(joinpath("figures","multineuron_and.png"), fig)
fig

##