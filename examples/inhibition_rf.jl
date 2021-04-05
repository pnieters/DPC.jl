using ADSP, CairoMakie
#include(joinpath(@__DIR__, "utils.jl"))
include("utils.jl")


fig = Figure(resolution = (800, 600))
## First setup: without inhibition

config1 = """
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
  θ_seg: 1
  branches: 
    - id: seg2
      θ_seg: 1
      θ_syn: 5
      branches: 
        - id: seg1
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

(net1,objects1) = load_network(YAML_source=config1)

volleys = [
    (70.0, :i1), (120.0, :i2), (170.0, :i3),
    (315.0, :i3), (325.0, :i2), (335.0, :i1), 
    (370.0, :i3), (380.0, :i2), (390.0, :i1), 
    (425.0, :i3), (435.0, :i2), (445.0, :i1) 
]
input=[ Event(:input_spikes, 0.0, t+τ-0.5, objects1[Symbol("$(pop)$(i)")]) for (t, pop) in volleys for (i,τ) in enumerate(5*rand(5))]

logger=simulate!(net1, input)

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
plateau1_ends = filter(x->(x.object == :seg1 && x.event == :plateau_ends), logger.data).t
plateau2_ends = filter(x->(x.object == :seg2 && x.event == :plateau_ends), logger.data).t
spike_times = filter(x->(x.object == :n && x.event == :spikes), logger.data).t
xticks,xtickformat = make_manual_ticks([0;250;500;plateau1_starts;plateau2_starts;spike_times], ["0ms";"250ms";"500ms";fill("t₀",length(plateau1_starts));fill("t₁",length(plateau2_starts));fill("t₂",length(spike_times))])
yticks,ytickformat = make_manual_ticks(0.5:14.5, reverse!(["$(grp)$(sub)" for grp in ["A","B","C"] for sub in ['₁','₂','₃','₄','₅']]))

ax11 = fig[1, 1] = Axis(fig, title = "Without inhibition", 
    height=Fixed(180),
    xticks=xticks, 
    xtickformat=xtickformat, 
    yticks=yticks, 
    ytickformat=ytickformat, 
    yticklabelsize=14
    )
ax12 = fig[1, 2] = Axis(fig, aspect=DataAspect(), backgroundcolor=:transparent)


for (i, syn) in enumerate([syn11,syn12,syn13,syn14,syn15])
    steps!(ax11, [0;syn.t;150], 9 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_1, color=:transparent)
end
steps!(ax11, [0;seg1.t;150], 10 .+ 2.45 .* [0;Int.(seg1.state);0], fill=color_1_50, color=:transparent)

for (i, syn) in enumerate([syn21,syn22,syn23,syn24,syn25])
    steps!(ax11, [0;syn.t;150], 4 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_2, color=:transparent)
end
steps!(ax11, [0;seg2.t;150], 5 .+ 5 .* ([0;Int.(seg2.state);0].>1), fill=color_2_50, color=:transparent)

for (i, syn) in enumerate([syn31,syn32,syn33,syn34,syn35])
    steps!(ax11, [0;syn.t;150], -1 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3, color=:transparent)
end
# steps!(ax11, [0;n.t;900], 0 .+ 2.45 .* [0;Int.(n.state);0], fill=color_3_50)
# steps!(ax11, [0;seg1.t;900], 1 .+ 0.45 .* [0;Int.(seg1.state);0], fill=color_1_50)


lines!.(ax11, Rect.(plateau1_starts .- 5.5, 9.95, 11, 5.1), color=:darkgray, linewidth=2)
lines!.(ax11, Rect.(plateau2_starts .- 5.5, 4.95, 11, 5.1), color=:darkgray, linewidth=2)
lines!.(ax11, Rect.(spike_times .- 5.5, -0.05, 11, 5.1), color=:darkgray, linewidth=2)

# arrows!(ax11, Point2f0[(spike_times[]-5.5,2.5)], Point2f0[(-59.5,0)], linewidth=2, arrowsize = 20, linecolor=:darkgray, arrowcolor=:darkgray)
# arrows!(ax11, Point2f0[(spike_times[]+5.5,2.5)], Point2f0[(59.5,0)], linewidth=2, arrowsize = -20, linecolor=:darkgray, arrowcolor=:darkgray)

linesegments!(ax11, repeat(spike_times, inner=2), repeat([-0.5,15.5], outer=length(spike_times)), linewidth=3, color=:gray10)
ylims!(ax11, (-0.5,15.5))

pn1 = plot!(ax12, objects1[:n], ports=Dict(:n=>[:C],:seg2=>[:dummy, :B],:seg1=>[:dummy, :A]), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_3, :seg1=>color_1, :seg2=>color_2))
hidedecorations!(ax12)

ports = Dict(pn1.attributes[:ports][])
arrows!(ax12, [ports[x] - Point2f0(0.4, 0) for x in [:C,:B,:A]] , fill(Point2f0(0.3,0),3), linewidth=2, arrowsize = 20)
for (label,pos)  in zip(["C","B","A"], [ports[x] - Point2f0(0.45, 0) for x in [:C,:B,:A]])
    text!(ax12, label, position=pos, align=(:right, :center), textsize=0.25, color=:black)
end

fig


## Second setup: with inhibition

config2 = """
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
  θ_seg: 1
  branches: 
    - id: seg2
      θ_seg: 1
      θ_syn: 5
      branches: 
        - id: seg1
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
- {id: syn41, source: i31, target: seg1, weight: -1, delay: 10.1}
- {id: syn42, source: i32, target: seg1, weight: -1, delay: 10.1}
- {id: syn43, source: i33, target: seg1, weight: -1, delay: 10.1}
- {id: syn44, source: i34, target: seg1, weight: -1, delay: 10.1}
- {id: syn45, source: i35, target: seg1, weight: -1, delay: 10.1}
- {id: syn51, source: i31, target: seg2, weight: -1, delay: 10.1}
- {id: syn52, source: i32, target: seg2, weight: -1, delay: 10.1}
- {id: syn53, source: i33, target: seg2, weight: -1, delay: 10.1}
- {id: syn54, source: i34, target: seg2, weight: -1, delay: 10.1}
- {id: syn55, source: i35, target: seg2, weight: -1, delay: 10.1}
"""

(net2,objects2) = load_network(YAML_source=config2)

volleys = [
    (70.0, :i1), (120.0, :i2), (170.0, :i3),
    (315.0, :i3), (325.0, :i2), (335.0, :i1), 
    (370.0, :i3), (380.0, :i2), (390.0, :i1), 
    (425.0, :i3), (435.0, :i2), (445.0, :i1) 
]
input=[ Event(:input_spikes, 0.0, t+τ-0.5, objects2[Symbol("$(pop)$(i)")]) for (t, pop) in volleys for (i,τ) in enumerate(5*rand(5))]

logger=simulate!(net2, input)

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
syn41 = get_trace(:syn41, logger.data)
syn42 = get_trace(:syn42, logger.data)
syn43 = get_trace(:syn43, logger.data)
syn44 = get_trace(:syn44, logger.data)
syn45 = get_trace(:syn45, logger.data)
syn51 = get_trace(:syn51, logger.data)
syn52 = get_trace(:syn52, logger.data)
syn53 = get_trace(:syn53, logger.data)
syn54 = get_trace(:syn54, logger.data)
syn55 = get_trace(:syn55, logger.data)
seg1  = get_trace(:seg1, logger.data)
seg2  = get_trace(:seg2, logger.data)
n     = get_trace(:n, logger.data)

plateau1_starts = filter(x->(x.object == :seg1 && x.event == :plateau_starts), logger.data).t
plateau2_starts = filter(x->(x.object == :seg2 && x.event == :plateau_starts), logger.data).t
plateau1_ends = filter(x->(x.object == :seg1 && x.event == :plateau_ends), logger.data).t
plateau2_ends = filter(x->(x.object == :seg2 && x.event == :plateau_ends), logger.data).t
spike_times = filter(x->(x.object == :n && x.event == :spikes), logger.data).t
xticks,xtickformat = make_manual_ticks([0;250;500;plateau1_starts;plateau2_starts;spike_times], ["0ms";"250ms";"500ms";fill("t₀",length(plateau1_starts));fill("t₁",length(plateau2_starts));fill("t₂",length(spike_times))])
yticks,ytickformat = make_manual_ticks(0.5:14.5, reverse!(["$(grp)$(sub)" for grp in ["A","B","C"] for sub in ['₁','₂','₃','₄','₅']]))

ax21 = fig[2, 1] = Axis(fig, title = "With inhibition", 
    height=Fixed(180),
    xticks=xticks, 
    xtickformat=xtickformat, 
    yticks=yticks, 
    ytickformat=ytickformat, 
    yticklabelsize=14
    )
ax22 = fig[2, 2] = Axis(fig, aspect=DataAspect(), backgroundcolor=:transparent)


for (i, syn) in enumerate([syn11,syn12,syn13,syn14,syn15])
    steps!(ax21, [0;syn.t;150], 9 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_1, color=:transparent)
end
steps!(ax21, [0;seg1.t;150], 10 .+ 2.45 .* [0;Int.(seg1.state);0], fill=color_1_50, color=:transparent)

for (i, syn) in enumerate([syn21,syn22,syn23,syn24,syn25])
    steps!(ax21, [0;syn.t;150], 4 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_2, color=:transparent)
end
steps!(ax21, [0;seg2.t;150], 5 .+ 5 .* ([0;Int.(seg2.state);0].>1), fill=color_2_50, color=:transparent)

for (i, syn) in enumerate([syn31,syn32,syn33,syn34,syn35])
    steps!(ax21, [0;syn.t;150], -1 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3, color=:transparent)
end
# steps!(ax21, [0;n.t;900], 0 .+ 2.45 .* [0;Int.(n.state);0], fill=color_3_50)
# steps!(ax21, [0;seg1.t;900], 1 .+ 0.45 .* [0;Int.(seg1.state);0], fill=color_1_50)


lines!.(ax21, Rect.(plateau1_starts .- 5.5, 9.95, 11, 5.1), color=:darkgray, linewidth=2)
lines!.(ax21, Rect.(plateau2_starts .- 5.5, 4.95, 11, 5.1), color=:darkgray, linewidth=2)
lines!.(ax21, Rect.(spike_times .- 5.5, -0.05, 11, 5.1), color=:darkgray, linewidth=2)

# arrows!(ax21, Point2f0[(spike_times[]-5.5,2.5)], Point2f0[(-59.5,0)], linewidth=2, arrowsize = 20, linecolor=:darkgray, arrowcolor=:darkgray)
# arrows!(ax21, Point2f0[(spike_times[]+5.5,2.5)], Point2f0[(59.5,0)], linewidth=2, arrowsize = -20, linecolor=:darkgray, arrowcolor=:darkgray)

linesegments!(ax21, repeat(spike_times, inner=2), repeat([-0.5,15.5], outer=length(spike_times)), linewidth=3, color=:gray10)
ylims!(ax21, (-0.5,15.5))

pn2 = plot!(ax22, objects2[:n], ports=Dict(:n=>[:C],:seg2=>[:dummy, :B, :nC2],:seg1=>[:dummy, :A, :nC1]), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_3, :seg1=>color_1, :seg2=>color_2))
hidedecorations!(ax22)

ports = Dict(pn2.attributes[:ports][])
arrows!(ax22, [ports[x] - Point2f0(0.4, 0) for x in [:C,:B,:A]] , fill(Point2f0(0.3,0),3), linewidth=2, arrowsize = 20)
arrows!(ax22, [ports[x] + Point2f0(0.4, 0) for x in [:nC2,:nC1]] , fill(Point2f0(-0.3,0),2), linewidth=2, arrowsize = -20)
for (label,pos)  in zip(["C","B","A"], [ports[x] - Point2f0(0.45, 0) for x in [:C,:B,:A]])
    text!(ax22, label, position=pos, align=(:right, :center), textsize=0.25, color=:black)
end
for (label,pos)  in zip(["-C","-C"], [ports[x] + Point2f0(0.45, 0) for x in [:nC2,:nC1]])
    text!(ax22, label, position=pos, align=(:left, :center), textsize=0.25, color=:black)
end

fig

## Save
colsize!(fig.layout, 2, Fixed(150))

save(joinpath("figures","inhibition_rf.pdf"), fig)
save(joinpath("figures","inhibition_rf.svg"), fig)
save(joinpath("figures","inhibition_rf.png"), fig)
fig
################################################################################