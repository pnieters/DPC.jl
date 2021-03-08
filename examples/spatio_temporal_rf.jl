using ADSP, CairoMakie
include("utils.jl")


color_1 = pal.colors[1]
color_1_50 = RGBAf0(color_1.r,color_1.g,color_1.b,0.5)
color_1_25 = RGBAf0(color_1.r,color_1.g,color_1.b,0.25)
color_2 = pal.colors[2]
color_2_50 = RGBAf0(color_2.r,color_2.g,color_2.b,0.5)
color_2_25 = RGBAf0(color_2.r,color_2.g,color_2.b,0.25)

config = """
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
  - {id: syn6, source: i5, target: n}
  - {id: syn7, source: i2, target: n}
  - {id: syn8, source: i4, target: n}
  - {id: syn9, source: i7, target: n}
  - {id: syn10, source: i9, target: n}
"""

(net,objects) = load_network(YAML_source=config)

inp=[
    Event(:input_spikes, 0.0, 16.0, objects[:i1]),
    Event(:input_spikes, 0.0, 17.0, objects[:i5]),
    Event(:input_spikes, 0.0, 18.0, objects[:i6]),
    Event(:input_spikes, 0.0, 19.0, objects[:i3]),
    Event(:input_spikes, 0.0, 20.0, objects[:i8]),
    Event(:input_spikes, 0.0, 66.0, objects[:i9]),
    Event(:input_spikes, 0.0, 67.0, objects[:i7]),
    Event(:input_spikes, 0.0, 68.0, objects[:i4]),
    Event(:input_spikes, 0.0, 69.0, objects[:i2]),
    Event(:input_spikes, 0.0, 70.0, objects[:i5]),
]

logger=simulate!(net, inp, 150.0)


fig = Figure(resolution = (800, 400), show_axis=false)
ax1 = fig[1, 1] = Axis(fig, title = "Neural spatio-temporal response", height=Fixed(100), axis=(xticks=[spikes.t;plateau_starts.t]))

seg = filter(x->(x.object==:seg && x.event==:backprop), logger.data)
steps!(ax1, [0;seg.t;150.0], 2.49*[0;Int.(seg.state);0] .- 5.05, color=:transparent, fill=color_1_25)
for (i,syn) in enumerate([:syn1, :syn2, :syn3, :syn4, :syn5])
    epsp = filter(x->(x.object==syn), logger.data)
    steps!(ax1, [0;epsp.t;150.0], -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=color_1_50)
end
spikes = filter(x->(x.object==:n && x.event==:spikes), logger.data)
lines!(spikes.t, [-10.1, -5.1], linewidth=3, color=color_1)

for (i,syn) in enumerate([:syn6, :syn7, :syn8, :syn9, :syn10])
    epsp = filter(x->(x.object==syn), logger.data)
    steps!(ax1, [0;epsp.t;150.0], -5.005 + -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=color_2_50)
end

plateau_starts = filter(x->(x.object==:seg && x.event==:plateau_starts), logger.data)
# lines!(spikes.t, [-10.1, -5.1], linewidth=5, color=:black)
# lines!(plateau_starts.t, [-5.1, -0.1], linewidth=5, color=:black)
lines!(ax1, Rect(16, -5.25, 11, 5.5), color=:darkgray, linewidth=2)
lines!(ax1, Rect(66, -10.25, 11, 5.5), color=:darkgray, linewidth=2)

arrows!(ax1, Point2f0[(66,-7.5)], Point2f0[(-44,0)], linewidth=2, arrowsize = 20, color=:darkgray)
arrows!(ax1, Point2f0[(77,-7.5)], Point2f0[(44,0)], linewidth=2, arrowsize = -20, color=:darkgray)

manual_ticks = [0;50;100;150;plateau_starts.t;spikes.t;plateau_starts.t.+objects[:seg].plateau_duration]
manual_labels= ["0","50","100","150","t₀","t₁","t₀+τ"]
idx = sortperm(manual_ticks)
ax1.xticks = manual_ticks[idx]
ax1.xtickformat = x->manual_labels[idx]
ax1.yticks = [-7.5, -2.5]
ax1.ytickformat = x->["soma", "seg."]

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

ax2 = fig[2, 1][1,1] = Axis(fig, title="Early RF", limits=Rect(-5,5,10,10), aspect = DataAspect(), height=Fixed(150))
heatmap!(ax2, x, y, A_rf.(x', y), colormap=:diverging_bwr_40_95_c42_n256, interpolate=true, colorrange=(-A_rf(0,0),A_rf(0,0)))
lines!.(ax2, Circle.(A_centers,σ), linewidth=2, color=:black)
hidedecorations!(ax2)

ax3 = fig[2, 1][1,2] = Axis(fig, title="Late RF", limits=Rect(-5,5,10,10), aspect = DataAspect(), height=Fixed(150))
heatmap!(ax3, x, y, B_rf.(x', y), colormap=:diverging_bwr_40_95_c42_n256, interpolate=true, colorrange=(-B_rf(0,0),B_rf(0,0)))
lines!.(ax3, Circle.(B_centers,σ), linewidth=2, color=:black)
hidedecorations!(ax3)


A_y = [ 1.2,  1.4,  1.6,  1.8,  2.0]
B_y = [-1.2, -1.4, -1.6, -1.8, -2.0]
A_nodes_1 = collect(zip([4,4,4,4,4], A_y))
A_nodes_2 = collect(zip([7,7,7,7,7], A_y))
B_nodes_1 = collect(zip([4,4,4,4,4], B_y))
B_nodes_2 = collect(zip([7,7,7,7,7], B_y))

ax4 = fig[2,1][1,3] = Axis(fig, title="Morphology", limits=Rect(-5,5,14,14), height=Fixed(150), aspect = DataAspect())
lines!(ax4, Rect(-5,-5,10,10), color=:black)
plot!(ax4, objects[:n], Dict{Symbol,Vector{Port}}(), branch_width=0.75, branch_length=5.0, root_position=Point2f0(7.0, -2.0), color=Dict(:seg=>color_1,:n=>color_2))
hidedecorations!(ax4)
hidespines!(ax4)

x = repeat(-2:1:2, outer=5)
y = repeat(-2:1:2, inner=5)

poly!.(ax4, Circle.(A_centers,σ), linestyle=:dash, linewidth=2, color=color_1_50)
poly!.(ax4, Circle.(B_centers,σ), linestyle=:dot, linewidth=2, color=color_2_50)
for i in 1:5
  lines!(ax4, Point2f0[A_centers[i], A_nodes_1[i], A_nodes_2[i]], linewidth=2, color=color_1)
  lines!(ax4, Point2f0[B_centers[i], B_nodes_1[i], B_nodes_2[i]], linewidth=2, color=color_2)
end
scatter!(ax4, x, y, marker=:x)

hidedecorations!(ax4)
hidespines!(ax4)


# fig
save("spatio_temporal_rf.svg", fig)