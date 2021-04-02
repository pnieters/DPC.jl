using ADSP, CairoMakie
include("utils.jl")

# draw_rf = true

draw_rf = true #false
set_theme!(mytheme)

config_1 = """
plateau_duration: 100.0
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
neurons:
- id: n
  θ_syn: 10
synapses:
- {id: syn11, source: i11, target: n}
- {id: syn12, source: i12, target: n}
- {id: syn13, source: i13, target: n}
- {id: syn14, source: i14, target: n}
- {id: syn15, source: i15, target: n}
- {id: syn21, source: i21, target: n}
- {id: syn22, source: i22, target: n}
- {id: syn23, source: i23, target: n}
- {id: syn24, source: i24, target: n}
- {id: syn25, source: i25, target: n}
"""

config_2 = """
plateau_duration: 100.0
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
neurons:
- id: n
  θ_syn: 5
  θ_seg: 1
  branches: 
    - id: seg
      θ_syn: 5
synapses:
- {id: syn11, source: i11, target: seg}
- {id: syn12, source: i12, target: seg}
- {id: syn13, source: i13, target: seg}
- {id: syn14, source: i14, target: seg}
- {id: syn15, source: i15, target: seg}
- {id: syn21, source: i21, target: n}
- {id: syn22, source: i22, target: n}
- {id: syn23, source: i23, target: n}
- {id: syn24, source: i24, target: n}
- {id: syn25, source: i25, target: n}
"""

(net_1,objects_1) = load_network(YAML_source=config_1)

volleys_1 = [(25.0, :i1), (25.0, :i2)]
input_1=[ Event(:input_spikes, 0.0, t+τ-0.5, objects_1[Symbol("$(pop)$(i)")]) for (t, pop) in volleys_1 for (i,τ) in enumerate(5*rand(5))]

volleys_2 = [(25.0, :i1), (75.0, :i2)]
input_2=[ Event(:input_spikes, 0.0, t+τ-0.5, objects_2[Symbol("$(pop)$(i)")]) for (t, pop) in volleys_2 for (i,τ) in enumerate(5*rand(5))]

logger_1=simulate!(net_1, input_1, 200.0)
spikes_1 = filter(x->(x.object==:n && x.event==:spikes), logger_1.data)

(net_2,objects_2) = load_network(YAML_source=config_2)


logger_2=simulate!(net_2, input_2, 200.0)
spikes_2 = filter(x->(x.object==:n && x.event==:spikes), logger_2.data)

plateau_starts = filter(x->(x.object==:seg && x.event==:plateau_starts), logger_2.data)


## Start plotting



A_starts = [(-4,-i+0.5) for i  in 1:5]
B_starts = [(-4,-5-i+0.5) for i  in 1:5]
A_joints = [(-2,-1.25-i*0.5) for i  in 1:5]
B_joints = [(-2,-6.25-i*0.5) for i  in 1:5]
A_ends = [(0,-1.25-i*0.5) for i  in 1:5]
B_ends = [(0,-6.25-i*0.5) for i  in 1:5]


################################################################################
## Neuron without active dendrite segment                                     ##
################################################################################
fig = Figure(resolution = (800, 400), show_axis=false)
xticks,xtickformat = make_manual_ticks([0;50;100;150;200;spikes_1.t], ["0ms","50ms","100ms","150ms","200ms","t₀"])
yticks,ytickformat = make_manual_ticks(-9.5:-0.5, reverse!(["$(grp)$(sub)" for grp in ["A","B"] for sub in ['₁','₂','₃','₄','₅']]))

# Draw response
grd = fig[1, 1] = GridLayout()
ax_res_1 = grd[1,1] = Axis(fig, title = "Point-neuron", height=Fixed(120), xticks=xticks, xtickformat=xtickformat, yticks=yticks, ytickformat=ytickformat, yticklabelsize=14)
ax_morph_1 = grd[1,2] = Axis(fig, width=Fixed(100), aspect=DataAspect(), backgroundcolor=:transparent)

seg = get_trace(:seg, logger_1.data)
steps!(ax_res_1, [0;seg.t;150.0], 2.49*[0;Int.(seg.state);0] .- 5.05, color=:transparent, fill=color_1_50)
for (i,syn) in enumerate([:syn11, :syn12, :syn13, :syn14, :syn15])
    epsp = get_trace(syn, logger_1.data)
    steps!(ax_res_1, [0;epsp.t;200.0], -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=color_1)
end

for (i,syn) in enumerate([:syn21, :syn22, :syn23, :syn24, :syn25])
  epsp = get_trace(syn, logger_1.data)
  steps!(ax_res_1, [0;epsp.t;200.0], -5.005 + -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=color_2)
end

lines!(ax_res_1, Rect(spikes_1.t[]-5.5, -10.05, 11, 10.1), color=:darkgray, linewidth=2)
lines!(ax_res_1, spikes_1.t, [-10.1, 0.1], linewidth=3, color=:black)

# arrows!(ax_res_1, Point2f0[(16,-5)], Point2f0[(-1,0)], linewidth=2, arrowsize = 20, color=:darkgray)
# arrows!(ax_res_1, Point2f0[(27,-5)], Point2f0[(1,0)], linewidth=2, arrowsize = -20, color=:darkgray)

ylims!(ax_res_1, (-10.5,0.5))
xlims!(ax_res_1, -5, 210)


# Draw morphology

for i in 1:5
  lines!(ax_morph_1, Point2f0[A_starts[i], A_joints[i], A_ends[i]], linewidth=2, color=color_1)
  arc!(ax_morph_1, (0,-5.25), abs(A_ends[i][2]+5.25), 3π/2, 5π/2, linewidth=2, color=color_1)
  lines!(ax_morph_1, Point2f0[B_starts[i], B_joints[i], B_ends[i]], linewidth=2, color=color_2)
end

plot!(ax_morph_1, objects_1[:n], root_position=Point2f0(0,-8), branch_width=1, branch_length=5.0, color=Dict(:n=>color_3))
hidedecorations!(ax_morph_1)

hidedecorations!(ax_morph_1)
hidespines!(ax_morph_1)
xlims!(ax_morph_1, (-4,4))
ylims!(ax_morph_1, (-10.5,0.5))
colgap!(grd, 1, 0)


################################################################################
## Neuron with active dendrite segment                                        ##
################################################################################
xticks,xtickformat = make_manual_ticks([0;50;100;150;200;plateau_starts.t;spikes_2.t;plateau_starts.t.+objects_2[:seg].plateau_duration],["0ms","50ms","100ms","150ms","200ms","t₀","t₁","t₀+τ"])

ax_res_2 = grd[2,1] = Axis(fig, title = "Neuron with active dendrite segment", height=Fixed(120), xticks=xticks, xtickformat=xtickformat, yticks=yticks, ytickformat=ytickformat, yticklabelsize=14)

seg = get_trace(:seg, logger_2.data)
steps!(ax_res_2, [0;seg.t;150.0], 2.49*[0;Int.(seg.state);0] .- 5.05, color=:transparent, fill=color_1_50)
for (i,syn) in enumerate([:syn11, :syn12, :syn13, :syn14, :syn15])
    epsp = get_trace(syn, logger_2.data)
    steps!(ax_res_2, [0;epsp.t;150.0], -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=color_1)
end

for (i,syn) in enumerate([:syn21, :syn22, :syn23, :syn24, :syn25])
  epsp = get_trace(syn, logger_2.data)
  steps!(ax_res_2, [0;epsp.t;150.0], -5.005 + -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=color_2)
end

lines!(ax_res_2, Rect(plateau_starts.t[]-5.5, -5.25, 11, 5.5), color=:darkgray, linewidth=2)
lines!(ax_res_2, Rect(spikes_2.t[]-5.5, -10.25, 11, 5.5), color=:darkgray, linewidth=2)
lines!(spikes_2.t, [-10.1, 0.1], linewidth=3, color=:black)
arrows!(ax_res_2, Point2f0[(spikes_2.t[]-5.5,-7.5)], Point2f0[(-44,0)], linewidth=2, arrowsize = 20, linecolor=:darkgray, arrowcolor=:darkgray)
arrows!(ax_res_2, Point2f0[(spikes_2.t[]+5.5,-7.5)], Point2f0[(44,0)], linewidth=2, arrowsize = -20, linecolor=:darkgray, arrowcolor=:darkgray)



# Draw morphology

ax_morph_2 = grd[2,2] = Axis(fig, width=Fixed(100), aspect=DataAspect(), backgroundcolor=:transparent)

for i in 1:5
  lines!(ax_morph_2, Point2f0[A_starts[i], A_joints[i], A_ends[i]], linewidth=2, color=color_1)
  lines!(ax_morph_2, Point2f0[B_starts[i], B_joints[i], B_ends[i]], linewidth=2, color=color_2)
end

plot!(ax_morph_2, objects_2[:n], root_position=Point2f0(0,-8), branch_width=1, branch_length=8.0, color=Dict(:n=>color_2, :seg=>color_1))
hidedecorations!(ax_morph_2)

hidedecorations!(ax_morph_2)
hidespines!(ax_morph_2)
xlims!(ax_morph_2, (-4,4))
ylims!(ax_morph_2, (-10.5,0.5))

ylims!(ax_res_2, (-10.5,0.5))
xlims!(ax_res_2, -5, 210)

################################################################################
## Draw receptive field                                                       ##
################################################################################
grd_rf = fig[1, 2] = GridLayout()
ax_schema = grd_rf[1:3,2] = Axis(fig,  backgroundcolor=:transparent, aspect=DataAspect())
ax_rf_1 = grd_rf[1,1] = Axis(fig)
ax_rf_2 = grd_rf[2,1] = Axis(fig)
ax_rf_3 = grd_rf[3,1] = Axis(fig)

config_3 = """
neurons:
  - id: n
    θ_syn: 5
    branches:
      - id: seg
        θ_syn: 5
        branches:
          - id: seg2
            θ_syn: 5
"""
(net_3,objects_3) = load_network(YAML_source=config_3)

plot!(ax_schema, objects_3[:n],  
    branch_width=1, 
    branch_length=8.0, 
    color=Dict(:n=>color_3, :seg=>color_2, :seg2=>color_1)
)
hidedecorations!(ax_schema)



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
hidedecorations!(ax_rf_1)

heatmap!(ax_rf_2, x, y, B_rf.(x', y), 
    colormap=:diverging_bwr_40_95_c42_n256, 
    interpolate=true, colorrange=(-B_rf(0,0),B_rf(0,0))
)
lines!(ax_rf_2, Point2f0[(-4,-2),(0,2)], linewidth=3, color=RGBAf0(0,0,0,0.5))
lines!(ax_rf_2, Point2f0[(-2,-2),(2,2)], linewidth=3)
hidedecorations!(ax_rf_2)

heatmap!(ax_rf_3, x, y, C_rf.(x', y), 
    colormap=:diverging_bwr_40_95_c42_n256, 
    interpolate=true, colorrange=(-C_rf(0,2),C_rf(0,2))
)
lines!(ax_rf_3, Point2f0[(-4,-2),(0,2)], linewidth=3, color=RGBAf0(0,0,0,0.25))
lines!(ax_rf_3, Point2f0[(-2,-2),(2,2)], linewidth=3, color=RGBAf0(0,0,0,0.5))
lines!(ax_rf_3, Point2f0[( 0,-2),(4,2)], linewidth=3)
hidedecorations!(ax_rf_3)

colsize!(grd_rf, 1, Aspect(1,1))
colsize!(grd_rf, 2, Aspect(1,1))
colsize!(fig.layout, 2, Relative(0.3))

# Save figures
save(joinpath("figures","spatio_temporal_rf.pdf"), fig)
save(joinpath("figures","spatio_temporal_rf.svg"), fig)
save(joinpath("figures","spatio_temporal_rf.png"), fig)
fig
