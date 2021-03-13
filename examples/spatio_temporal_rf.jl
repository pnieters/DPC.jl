using ADSP, CairoMakie
include("utils.jl")


config_1 = """
refractory_duration: 5.01
plateau_duration: 100
inputs:
  - id: i1
  - id: i2
  - id: i3
  - id: i4
  - id: i5
  - id: i6
  - id: i7
  - id: i8
  - id: i9
  - id: i10
neurons:
  - id: n
    θ_syn: 10
synapses:
  - {id: syn1, source: i1, target: n}
  - {id: syn2, source: i3, target: n}
  - {id: syn3, source: i8, target: n}
  - {id: syn4, source: i6, target: n}
  - {id: syn5, source: i5, target: n}
  - {id: syn6, source: i10, target: n}
  - {id: syn7, source: i2, target: n}
  - {id: syn8, source: i4, target: n}
  - {id: syn9, source: i7, target: n}
  - {id: syn10, source: i9, target: n}
"""

config_2 = """
refractory_duration: 5.01
plateau_duration: 100
inputs:
  - id: i1
  - id: i2
  - id: i3
  - id: i4
  - id: i5
  - id: i6
  - id: i7
  - id: i8
  - id: i9
  - id: i10
neurons:
  - id: n
    θ_syn: 5
    branches:
    - id: seg
      θ_syn: 5
synapses:
  - {id: syn1, source: i1, target: seg}
  - {id: syn2, source: i3, target: seg}
  - {id: syn3, source: i8, target: seg}
  - {id: syn4, source: i6, target: seg}
  - {id: syn5, source: i5, target: seg}
  - {id: syn6, source: i10, target: n}
  - {id: syn7, source: i2, target: n}
  - {id: syn8, source: i4, target: n}
  - {id: syn9, source: i7, target: n}
  - {id: syn10, source: i9, target: n}
"""

(net_1,objects_1) = load_network(YAML_source=config_1)

inp_1=[
    Event(:input_spikes, 0.0, 16.0, objects_1[:i1]),
    Event(:input_spikes, 0.0, 17.0, objects_1[:i5]),
    Event(:input_spikes, 0.0, 18.0, objects_1[:i6]),
    Event(:input_spikes, 0.0, 19.0, objects_1[:i3]),
    Event(:input_spikes, 0.0, 20.0, objects_1[:i8]),
    Event(:input_spikes, 0.0, 16.0, objects_1[:i9]),
    Event(:input_spikes, 0.0, 17.0, objects_1[:i7]),
    Event(:input_spikes, 0.0, 18.0, objects_1[:i4]),
    Event(:input_spikes, 0.0, 19.0, objects_1[:i2]),
    Event(:input_spikes, 0.0, 20.0, objects_1[:i10]),
]

logger_1=simulate!(net_1, inp_1, 150.0)
spikes_1 = filter(x->(x.object==:n && x.event==:spikes), logger_1.data)

(net_2,objects_2) = load_network(YAML_source=config_2)

inp_2=[
    Event(:input_spikes, 0.0, 16.0, objects_2[:i1]),
    Event(:input_spikes, 0.0, 17.0, objects_2[:i5]),
    Event(:input_spikes, 0.0, 18.0, objects_2[:i6]),
    Event(:input_spikes, 0.0, 19.0, objects_2[:i3]),
    Event(:input_spikes, 0.0, 20.0, objects_2[:i8]),
    Event(:input_spikes, 0.0, 66.0, objects_2[:i9]),
    Event(:input_spikes, 0.0, 67.0, objects_2[:i7]),
    Event(:input_spikes, 0.0, 68.0, objects_2[:i4]),
    Event(:input_spikes, 0.0, 69.0, objects_2[:i2]),
    Event(:input_spikes, 0.0, 70.0, objects_2[:i10]),
]

logger_2=simulate!(net_2, inp_2, 150.0)
spikes_2 = filter(x->(x.object==:n && x.event==:spikes), logger_2.data)

plateau_starts = filter(x->(x.object==:seg && x.event==:plateau_starts), logger_2.data)


## Start plotting

fig = Figure(resolution = (800, 600), show_axis=false)


A_starts = [(-4,-i+0.5) for i  in 1:5]
B_starts = [(-4,-5-i+0.5) for i  in 1:5]
A_joints = [(-2,-1.25-i*0.5) for i  in 1:5]
B_joints = [(-2,-6.25-i*0.5) for i  in 1:5]
A_ends = [(0,-1.25-i*0.5) for i  in 1:5]
B_ends = [(0,-6.25-i*0.5) for i  in 1:5]

function force_ticks(manual_ticks, manual_labels)
  @assert length(manual_ticks) == length(manual_labels)
  idx = sortperm(manual_ticks)
  t = copy(manual_ticks)
  l = copy(manual_labels)
  t[idx], x->l[idx]
end

################################################################################
## Neuron without active dendrite segment                                     ##
################################################################################
xticks,xtickformat = force_ticks([0;50;100;150;spikes_1.t], ["0ms","50ms","100ms","150ms","t₀"])
yticks,ytickformat = force_ticks(-9.5:-0.5, reverse!(["$(grp)$(sub)" for grp in ["A","B"] for sub in ['₁','₂','₃','₄','₅']]))

# Draw response
grd = fig[1, 1] = GridLayout()
ax_res_1 = grd[1,1] = Axis(fig, title = "Point-neuron", height=Fixed(120), xticks=xticks, xtickformat=xtickformat, yticks=yticks, ytickformat=ytickformat)
ax_morph_1 = grd[1,2] = Axis(fig, width=Fixed(100), aspect=DataAspect(), backgroundcolor=:transparent)

seg = get_trace(:seg, logger_1.data)
steps!(ax_res_1, [0;seg.t;150.0], 2.49*[0;Int.(seg.state);0] .- 5.05, color=:transparent, fill=color_1_25)
for (i,syn) in enumerate([:syn1, :syn2, :syn3, :syn4, :syn5])
    epsp = get_trace(syn, logger_1.data)
    steps!(ax_res_1, [0;epsp.t;150.0], -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=color_1_50)
end

for (i,syn) in enumerate([:syn6, :syn7, :syn8, :syn9, :syn10])
  epsp = get_trace(syn, logger_1.data)
  steps!(ax_res_1, [0;epsp.t;150.0], -5.005 + -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=color_2_50)
end

lines!(ax_res_1, Rect(16, -10.05, 11, 10.1), color=:darkgray, linewidth=2)
lines!(ax_res_1, spikes_1.t, [-10.1, 0.1], linewidth=3, color=:black)

# arrows!(ax_res_1, Point2f0[(16,-5)], Point2f0[(-1,0)], linewidth=2, arrowsize = 20, color=:darkgray)
# arrows!(ax_res_1, Point2f0[(27,-5)], Point2f0[(1,0)], linewidth=2, arrowsize = -20, color=:darkgray)

ylims!(ax_res_1, (-10.5,0.5))



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
xticks,xtickformat = force_ticks([0;50;100;150;plateau_starts.t;spikes_2.t;plateau_starts.t.+objects_2[:seg].plateau_duration],["0ms","50ms","100ms","150ms","t₀","t₁","t₀+τ"])

ax_res_2 = grd[2,1] = Axis(fig, title = "Neuron with active dendrite", height=Fixed(120), xticks=xticks, xtickformat=xtickformat, yticks=yticks, ytickformat=ytickformat)

seg = get_trace(:seg, logger_2.data)
steps!(ax_res_2, [0;seg.t;150.0], 2.49*[0;Int.(seg.state);0] .- 5.05, color=:transparent, fill=color_1_25)
for (i,syn) in enumerate([:syn1, :syn2, :syn3, :syn4, :syn5])
    epsp = get_trace(syn, logger_2.data)
    steps!(ax_res_2, [0;epsp.t;150.0], -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=color_1_50)
end

for (i,syn) in enumerate([:syn6, :syn7, :syn8, :syn9, :syn10])
  epsp = get_trace(syn, logger_2.data)
  steps!(ax_res_2, [0;epsp.t;150.0], -5.005 + -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=color_2_50)
end

lines!(ax_res_2, Rect(16, -5.25, 11, 5.5), color=:darkgray, linewidth=2)
lines!(ax_res_2, Rect(66, -10.25, 11, 5.5), color=:darkgray, linewidth=2)
lines!(spikes_2.t, [-10.1, 0.1], linewidth=3, color=:black)

arrows!(ax_res_2, Point2f0[(66,-7.5)], Point2f0[(-44,0)], linewidth=2, arrowsize = 20, linecolor=:darkgray, arrowcolor=:darkgray)
arrows!(ax_res_2, Point2f0[(77,-7.5)], Point2f0[(44,0)], linewidth=2, arrowsize = -20, linecolor=:darkgray, arrowcolor=:darkgray)



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

################################################################################
## Draw receptive field example                                               ##
################################################################################

function make_rf((x₀,y₀), σ)
    function rf(x,y)
        r = ((x-x₀)^2 + (y-y₀)^2)/(2*σ^2)
        inv(π*σ^4)*(1-r)*exp(-r)
    end
end


σ = 1.0
r₀ = √(2)*σ # zero-crossing of receptive field

A_centers = Point2f0[(-2,-2),(-1,-1),(0,0),(1,1),(2,2)]
B_centers = Point2f0[(-2,2),(-1,1),(0,0),(1,-1),(2,-2)]

A_rfs = make_rf.(A_centers,  σ)
B_rfs = make_rf.(B_centers,  σ)

A_rf = (x,y) -> [0.5, 1.0, 1.2, 1.0, 0.5]' * map(f->f(x,y), A_rfs)
B_rf = (x,y) -> [0.5, 1.0, 1.2, 1.0, 0.5]' * map(f->f(x,y), B_rfs)

x = LinRange(-5,5,100)
y = LinRange(-5,5,100)
xx = repeat(-2:1:2, outer=5)
yy = repeat(-2:1:2, inner=5)


ax2 = fig[2, 1][1,1] = Axis(fig, title="Early RF", limits=Rect(-5,5,10,10), aspect = DataAspect(), height=Fixed(150))
heatmap!(ax2, x, y, A_rf.(x', y), colormap=:diverging_bwr_40_95_c42_n256, interpolate=true, colorrange=(-A_rf(0,0),A_rf(0,0)))
scatter!(ax2, xx, yy, marker=:x, color=:gray)
lines!.(ax2, decompose.(Point2f0, Circle.(A_centers,σ)), linewidth=2, color=:black)
hidedecorations!(ax2)

ax3 = fig[2, 1][1,2] = Axis(fig, title="Late RF", limits=Rect(-5,5,10,10), aspect = DataAspect(), height=Fixed(150))
heatmap!(ax3, x, y, B_rf.(x', y), colormap=:diverging_bwr_40_95_c42_n256, interpolate=true, colorrange=(-B_rf(0,0),B_rf(0,0)))
scatter!(ax3, xx, yy, marker=:x, color=:gray)
lines!.(ax3, decompose.(Point2f0, Circle.(B_centers,σ)), linewidth=2, color=:black)
hidedecorations!(ax3)


A_y = [ 1.2,  1.4,  1.6,  1.8,  2.0]
B_y = [-1.2, -1.4, -1.6, -1.8, -2.0]
A_nodes_1 = collect(zip([4,4,4,4,4], A_y))
A_nodes_2 = collect(zip([7,7,7,7,7], A_y))
B_nodes_1 = collect(zip([4,4,4,4,4], B_y))
B_nodes_2 = collect(zip([7,7,7,7,7], B_y))

ax4 = fig[2,1][1,3] = Axis(fig, title="RF connectome", height=Fixed(150), aspect = DataAspect(), backgroundcolor=:transparent)
lines!(ax4, Rect(-5,-5,10,10), color=:black)
plot!(ax4, objects_2[:n], branch_width=0.75, branch_length=5.0, root_position=Point2f0(7.0, -2.0), color=Dict(:seg=>color_1,:n=>color_2))
hidedecorations!(ax4)
hidespines!(ax4)

poly!.(ax4, decompose.(Point2f0, Circle.(A_centers,σ)), linestyle=:dash, linewidth=2, color=color_1_50)
poly!.(ax4, decompose.(Point2f0, Circle.(B_centers,σ)), linestyle=:dot, linewidth=2, color=color_2_50)
for i in 1:5
  lines!(ax4, Point2f0[A_centers[i], A_nodes_1[i], A_nodes_2[i]], linewidth=2, color=color_1)
  lines!(ax4, Point2f0[B_centers[i], B_nodes_1[i], B_nodes_2[i]], linewidth=2, color=color_2)
end
scatter!(ax4, xx, yy, marker=:x, color=:gray)

hidedecorations!(ax4)
hidespines!(ax4)
xlims!(ax4, (-5,9))
ylims!(ax4, (-5,5))

# Save figures
save(joinpath("figures","spatio_temporal_rf.pdf"), fig)
save(joinpath("figures","spatio_temporal_rf.svg"), fig)
save(joinpath("figures","spatio_temporal_rf.png"), fig)
fig