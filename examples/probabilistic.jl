using ADSP, CairoMakie
include("utils.jl")

config0 = """
neurons:
  - id: C
    θ_syn: 0
    branches:
      - id: A
        θ_syn: 5
"""

config1 = """
neurons:
  - id: C
    θ_syn: 0
    θ_seg: 2
    branches:
      - id: A
        θ_syn: 5
      - id: B
        θ_syn: 5
"""

config2 = """
neurons:
  - id: C
    θ_syn: 0
    θ_seg: 1
    branches:
      - id: A
        θ_syn: 5
      - id: B
        θ_syn: 5
"""

num_synapses = 10
prob_A = fill(0.0, num_synapses)
prob_B = fill(0.0, num_synapses)
volleys = Dict(:A => (prob_A, ones(num_synapses)), :B => (prob_B, ones(num_synapses)))

probabilities = 0:0.01:1

firing_probabilities0 = Vector{Vector{Float64}}(undef, num_synapses)
(net0,objects0) = load_network(YAML_source=config0)
for θ in 1:num_synapses
    println("Computing for single segment, θ=$(θ)")
    objects0[:A].θ_syn = θ
    firing_probabilities0[θ] = [(prob_A .= a; P_fire(objects0[:C], volleys)) for a ∈ probabilities]
end

println("Computing for AND neuron")
(net1,objects1) = load_network(YAML_source=config1)
firing_probabilities1 = [(prob_A .= a; prob_B .= b; P_fire(objects1[:C], volleys)) for a ∈ probabilities, b ∈ probabilities]

println("Computing for OR neuron")
(net2,objects2) = load_network(YAML_source=config2)
firing_probabilities2 = [(prob_A .= a; prob_B .= b; P_fire(objects2[:C], volleys)) for a ∈ probabilities, b ∈ probabilities]


fig = Figure(resolution=(600,600))
ax0 = fig[1,2:3] = Axis(fig, title="Single segment", xlabel="Synaptic transmission probability P₀", ylabel="Plateau probability")#, aspect=DataAspect()
cell_col_1 = fig[1,4]
ax0b = fig[1,1] = Axis(fig, aspect=DataAspect(), backgroundcolor=:transparent)
ax3 = fig[2,1] = Axis(fig, aspect=DataAspect(), backgroundcolor=:transparent)
ax1 = fig[2,2] = Axis(fig, xlabel="P₁", ylabel="P₂", title="'AND' (Θ=2)")
ax2 = fig[2,3] = Axis(fig, xlabel="P₁", ylabel="P₂", title="'OR' (Θ=1)")
cell_col_2 = fig[2,4]

colors = cgrad(:viridis,length(firing_probabilities0), categorical=true)
for (prob,col) in zip(firing_probabilities0, colors)
    lines!(ax0, probabilities, prob, color=col, linewidth=3)
end
ax_col_1 = Colorbar(cell_col_1, width=10, colormap=colors, limits=(0.5,10.5), ticks=1:10, label="Threshold (out of $(num_synapses) synapses)")

contourf!(ax1, probabilities, probabilities, firing_probabilities1, colormap=:viridis, colorrange=(0,1))

contourf!(ax2, probabilities, probabilities, firing_probabilities2, colormap=:viridis, colorrange=(0,1))

p_n=plot!(ax3, objects1[:C], ports=Dict(:C=>[:C],:A=>[:A],:B=>[:B]), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:C=>:gray, :A=>color_1, :B=>color_2))
ports = lift(x->last.(x) .- Point2f0[(-0.4,0),(0.4,0)], p_n.attributes[:ports])
arrows!(ax3, ports , Node([Point2f0(-0.3,0),Point2f0(0.3,0)]), linewidth=2, arrowsize = [-20,20])
text!.(ax3, "P₁", position=(@lift $ports[1]), align=(:right, :bottom), textsize=0.25, color=:black)
text!.(ax3, "P₂", position=(@lift $ports[2]), align=(:left, :bottom), textsize=0.25, color=:black)

p_n=plot!(ax0b, objects0[:C], ports=Dict(:C=>[:C],:A=>[:A]), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:C=>:gray, :A=>color_1, :B=>color_2))
port = lift(x->[x[1][2] - Point2f0(-0.4,0)], p_n.attributes[:ports])
arrows!(ax0b, port , [Point2f0(-0.3,0)], linewidth=2, arrowsize = -20)
text!.(ax0b, "P₀", position=(@lift $ports[1]), align=(:right, :bottom), textsize=0.25, color=:black)


ax_col_2 = Colorbar(cell_col_2, width=10, colormap=:viridis, limits=(0,1), label="Plateau probability")

hidedecorations!(ax3)
hidespines!(ax3)
hidedecorations!(ax0b)
hidespines!(ax0b)

hideydecorations!(ax2)
ax1.xticks = [0,0.5,1.0]
ax2.xticks=[0,0.5,1.0]
ax1.yticks=[0,0.5,1.0]
# ax3.ticks = []
rowsize!(fig.layout, 2, Aspect(2,1)) 
colsize!(fig.layout, 1, Fixed(75)) 

ax0.xticks = [0,0.5,1.0]
ax0.yticks = [0,0.5,1.0]
ax1.xticks = [0,0.5,1.0]
ax1.yticks = [0,0.5,1.0]
ax2.xticks = [0,0.5,1.0]
ax2.yticks = [0,0.5,1.0]

ax0.limits[] = Rect(0.0f0,0.0f0,1.0f0,1.0f0)
ax1.limits[] = Rect(0.0f0,0.0f0,1.0f0,1.0f0)
ax2.limits[] = Rect(0.0f0,0.0f0,1.0f0,1.0f0)
ax1.yticks[] = [0,0.5,1.0]



save(joinpath("figures","probabilistic.pdf"), fig)
save(joinpath("figures","probabilistic.svg"), fig)
save(joinpath("figures","probabilistic.png"), fig)
fig
