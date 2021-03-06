using ADSP, GLMakie
include("utils.jl")

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
    Event(:input_spikes, 0.0, 96.0, objects[:i9]),
    Event(:input_spikes, 0.0, 97.0, objects[:i7]),
    Event(:input_spikes, 0.0, 98.0, objects[:i4]),
    Event(:input_spikes, 0.0, 99.0, objects[:i2]),
    Event(:input_spikes, 0.0, 100.0, objects[:i5]),
]

logger=simulate!(net, inp, 150.0)


fig = Figure(resolution = (1200, 900), show_axis=false)
ax1 = fig[1, 1:2] = Axis(fig, title = "Neuron with spatio-temporal receptive field", aspect = AxisAspect(10))
hideydecorations!(ax1)

seg = filter(x->(x.object==:seg && x.event==:backprop), logger.data)
steps!(ax1, [0;seg.t;150.0], 2.49*[0;Int.(seg.state);0] .- 5.05, color=:transparent, fill=RGBAf0(1,0,0,0.25))
for (i,syn) in enumerate([:syn1, :syn2, :syn3, :syn4, :syn5])
    epsp = filter(x->(x.object==syn), logger.data)
    steps!(ax1, [0;epsp.t;150.0], -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=RGBAf0(1,0,0,0.5))
end
lines!(spikes.t, [-10.1, -5.1], linewidth=3, color=RGBAf0(1,0,0,1))

for (i,syn) in enumerate([:syn6, :syn7, :syn8, :syn9, :syn10])
    epsp = filter(x->(x.object==syn), logger.data)
    steps!(ax1, [0;epsp.t;150.0], -5.005 + -i .+ 0.9 .* [0;Int.(epsp.state);0], color=:transparent, fill=RGBAf0(0,0,1,0.5))
end

spikes = filter(x->(x.object==:n && x.event==:spikes), logger.data)
lines!(spikes.t, [-10.1, -5.1], linewidth=5, color=:black)


lines!(ax1, Rect(16, -5.25, 11, 5.5), color=:darkgray)
lines!(ax1, Rect(96, -10.25, 11, 5.5), color=:darkgray)

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

ax2 = fig[2, 1] = Axis(fig, title="Early receptive field", limits=Rect(-5,5,10,10), aspect = AxisAspect(1))
heatmap!(ax2, x, y, A_rf.(x', y), colormap=:diverging_bwr_40_95_c42_n256, interpolate=true, colorrange=(-A_rf(0,0),A_rf(0,0)))
lines!.(ax2, Circle.(A_centers,σ), linewidth=2, color=:black)
hidedecorations!(ax2)

ax3 = fig[2, 2] = Axis(fig, title="Late receptive field", limits=Rect(-5,5,10,10), aspect = AxisAspect(1))
heatmap!(ax3, x, y, B_rf.(x', y), colormap=:diverging_bwr_40_95_c42_n256, interpolate=true, colorrange=(-B_rf(0,0),B_rf(0,0)))
lines!.(ax3, Circle.(B_centers,σ), linewidth=2, color=:black)
hidedecorations!(ax3)

#scene = Scene(camera = cam2d!, limits=Rect(-5,5,-5,5,-5,5), show_axis = false)

# subscene1 = Scene(show_axis = false, camera=cam3d_cad!, limits=Rect(-5,5,-5,5,-5,5), backgroundcolor=:blue)
# heatmap!(subscene1, x, y, A_rf.(x', y), colormap=:diverging_bwr_40_95_c42_n256, colorrange=(-A_rf(0,0),A_rf(0,0)))
# rotate!(subscene1, 1,0,0)
# update_cam!(subscene1)
# cam1 = cameracontrols(subscene1)
# cam1.lookat[] = Float32[2,2,-0.5]
# cam1.eyeposition[] = Float32[7.5,-15,0]
# update_cam!(subscene1)
# lines!.(subscene1, Circle.(A_centers,σ), linewidth=2, color=:red)
# # lines!.(subscene1, Circle.(A_centers,r₂), opacity=0.5, linewidth=2, linestyle=:dash)
#
# update_cam!(subscene1)

# subscene2 = Scene(ax2.scene, show_axis = false)
# heatmap!(subscene2, x, y, B_rf.(x', y), colormap=:diverging_bwr_40_95_c42_n256, colorrange=(-B_rf(0,0),B_rf(0,0)))
# lines!.(subscene2, Circle.(B_centers,σ), linewidth=2)
# # lines!.(subscene2, Circle.(B_centers,r₂), opacity=0.5, linewidth=2, linestyle=:dash)
# #rotate!(subscene2, 0,1,0)
# translate!(subscene2, 10, 0)

# # cam = cameracontrols(scene)
ax4 = fig[:, 3] = Axis(fig)
plot!(ax4, objects[:n])

fig
