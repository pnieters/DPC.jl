using ADSP, CairoMakie
#include(joinpath(@__DIR__, "utils.jl"))
include("utils.jl")

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
  θ_seg: 1
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

volleys = [(15.0, :i1), (45.0, :i2), (80.0, :i3)]
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
plateau1_ends = filter(x->(x.object == :seg1 && x.event == :plateau_ends), logger.data).t
plateau2_ends = filter(x->(x.object == :seg2 && x.event == :plateau_ends), logger.data).t
spike_times = filter(x->(x.object == :n && x.event == :spikes), logger.data).t
xticks,xtickformat = make_manual_ticks([0;100;plateau1_starts;plateau2_starts;plateau1_ends;plateau2_ends;spike_times], ["0ms","100ms","t₀","t₁","t₀+τ","t₁+τ","t₂"])
yticks,ytickformat = make_manual_ticks(0.5:14.5, reverse!(["$(grp)$(sub)" for grp in ["A","B","C"] for sub in ['₁','₂','₃','₄','₅']]))

fig = Figure(resolution = (800, 400))
grd = fig[1, 1] = GridLayout()
ax11 = grd[1, 1] = Axis(fig, title = "Neuron with two active dendrite segments", 
    height=Fixed(180),
    xticks=xticks, 
    xtickformat=xtickformat, 
    yticks=yticks, 
    ytickformat=ytickformat, 
    yticklabelsize=14)
ax12 = grd[1, 2] = Axis(fig, width=Fixed(50), aspect=DataAspect(), backgroundcolor=:transparent)


for (i, syn) in enumerate([syn11,syn12,syn13,syn14,syn15])
    steps!(ax11, [0;syn.t;150], 9 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_1, color=:transparent)
end
steps!(ax11, [0;seg1.t;150], 10 .+ 2.45 .* [0;Int.(seg1.state);0], fill=color_1_50, color=:transparent)

for (i, syn) in enumerate([syn21,syn22,syn23,syn24,syn25])
    steps!(ax11, [0;syn.t;150], 4 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_2, color=:transparent)
end
steps!(ax11, [0;seg2.t;150], 5 .+ 2.45 .* [0;Int.(seg2.state);0], fill=color_2_50, color=:transparent)

for (i, syn) in enumerate([syn31,syn32,syn33,syn34,syn35])
    steps!(ax11, [0;syn.t;150], -1 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3, color=:transparent)
end
# steps!(ax11, [0;n.t;900], 0 .+ 2.45 .* [0;Int.(n.state);0], fill=color_3_50)
# steps!(ax11, [0;seg1.t;900], 1 .+ 0.45 .* [0;Int.(seg1.state);0], fill=color_1_50)


lines!(ax11, Rect(plateau1_starts[]-5.5, 9.95, 11, 5.1), color=:gray10, linewidth=2)
lines!(ax11, Rect(plateau2_starts[]-5.5, 4.95, 11, 5.1), color=:gray10, linewidth=2)
lines!(ax11, Rect(spike_times[]-5.5, -0.05, 11, 5.1), color=:gray10, linewidth=2)

arrows!(ax11, Point2f0[(spike_times[]-5.5,2.5)], Point2f0[(-59.5,0)], linewidth=2, arrowsize = 20, linecolor=:gray10, arrowcolor=:gray10)
arrows!(ax11, Point2f0[(spike_times[]+5.5,2.5)], Point2f0[(59.5,0)], linewidth=2, arrowsize = -20, linecolor=:gray10, arrowcolor=:gray10)

linesegments!(ax11, repeat(spike_times, inner=2), repeat([0,16], outer=length(spike_times)), linewidth=3, linestyle=:dash, color=:gray10)
ylims!(ax11, (-0.5,15.5))

plot!(ax12, objects[:n], angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:n=>color_3, :seg1=>color_1, :seg2=>color_2))


hidedecorations!(ax12)


################################################################################
## Draw receptive field                                                       ##
################################################################################
#=
grd_rf = fig[1:2, 2] = GridLayout()

ax_rf_1 = grd_rf[1,1] = Axis(fig, aspect=DataAspect())
ax_rf_2 = grd_rf[1,2] = Axis(fig, aspect=DataAspect())
ax_rf_3 = grd_rf[1,3] = Axis(fig, aspect=DataAspect())

config_3 = """
neurons:
  - id: n
    branches:
      - id: seg1
      - id: seg2
      - id: seg3
"""
(net_3,objects_3) = load_network(YAML_source=config_3)




function make_rf((x₀,y₀), σ)
  function rf(x,y)
    r = ((x-x₀)^2 + (y-y₀)^2)/(2*σ^2)
    inv(π*σ^4)*(1-r)*exp(-r)
  end
end

σ = 1.0
r₀ = √(2)*σ # zero-crossing of receptive field

A_centers = Point2f0[(-2,-4),(-1,-3),( 0,-2),( 1,-1),( 2, 0)]
B_centers = Point2f0[(-2,-2),(-1,-1),( 0, 0),( 1, 1),( 2, 2)]
C_centers = Point2f0[(-2, 0),(-1, 1),( 0, 2),( 1, 3),( 2, 4)]

A_rfs = make_rf.(A_centers,  σ)
B_rfs = make_rf.(B_centers,  σ)
C_rfs = make_rf.(C_centers,  σ)

A_rf = (x,y) -> [0.5, 1.0, 1.2, 1.0, 0.5]' * map(f->f(x,y), A_rfs)
B_rf = (x,y) -> [0.5, 1.0, 1.2, 1.0, 0.5]' * map(f->f(x,y), B_rfs)
C_rf = (x,y) -> [0.5, 1.0, 1.2, 1.0, 0.5]' * map(f->f(x,y), C_rfs)

x = LinRange(-5,5,100)
y = LinRange(-5,5,100)
xx = repeat(-2:1:2, outer=5)
yy = repeat(-2:1:2, inner=5)

heatmap!(ax_rf_1, x, y, A_rf.(x', y), 
    colormap=:diverging_bwr_40_95_c42_n256, 
    interpolate=true, colorrange=(-A_rf(0,-2),A_rf(0,-2))
)
lines!(ax_rf_1, Point2f0[(-4,-2),(0,2)], linewidth=3)
lines!(ax_rf_1, Point2f0[(-2,-2),(2,2)], linewidth=2, linestyle=:dash, color=:gray50)
lines!(ax_rf_1, Point2f0[( 0,-2),(4,2)], linewidth=2, linestyle=:dash, color=:gray50)
hidedecorations!(ax_rf_1)

heatmap!(ax_rf_2, x, y, B_rf.(x', y), 
    colormap=:diverging_bwr_40_95_c42_n256, 
    interpolate=true, colorrange=(-B_rf(0,0),B_rf(0,0))
)
lines!(ax_rf_2, Point2f0[(-4,-2),(0,2)], linewidth=2, linestyle=:dash, color=:gray50)
lines!(ax_rf_2, Point2f0[(-2,-2),(2,2)], linewidth=3)
lines!(ax_rf_2, Point2f0[( 0,-2),(4,2)], linewidth=2, linestyle=:dash, color=:gray50)
hidedecorations!(ax_rf_2)

heatmap!(ax_rf_3, x, y, C_rf.(x', y), 
    colormap=:diverging_bwr_40_95_c42_n256, 
    interpolate=true, colorrange=(-C_rf(0,2),C_rf(0,2))
)
lines!(ax_rf_3, Point2f0[(-4,-2),(0,2)], linewidth=2, linestyle=:dash, color=:gray50)
lines!(ax_rf_3, Point2f0[(-2,-2),(2,2)], linewidth=2, linestyle=:dash, color=:gray50)
lines!(ax_rf_3, Point2f0[( 0,-2),(4,2)], linewidth=3)
hidedecorations!(ax_rf_3)

ax_schema = grd_rf[2,:] = Axis(fig,  backgroundcolor=:transparent, aspect=DataAspect())
plot!(ax_schema, objects_3[:n],  
    branch_width=1, 
    branch_length=8.0, 
    angle_between=30.0/180*π,
    color=Dict(:n=>:gray50, :seg3=>color_3, :seg1=>color_2, :seg2=>color_1)
)
hidedecorations!(ax_schema)
xlims!(ax_schema, -10,10)
ylims!(ax_schema, -2,10)

# colsize!(grd_rf, 1, Aspect(1,1))
# colsize!(grd_rf, 2, Aspect(1,1))
colsize!(fig.layout, 2, Relative(0.45))
colgap!(fig.layout, 1, Fixed(50))
rowsize!(grd_rf, 1, Aspect(1,1))
=#
## Save

save(joinpath("figures","invariant_rf.pdf"), fig)
save(joinpath("figures","invariant_rf.svg"), fig)
save(joinpath("figures","invariant_rf.png"), fig)
fig
################################################################################