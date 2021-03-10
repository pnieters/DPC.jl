using ADSP, CairoMakie
include("utils.jl")

config0 = """
neurons:
  - id: C
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
prob_C = fill(0.0, num_synapses)
volleys = Dict(:A => (prob_A, ones(num_synapses)), :B => (prob_B, ones(num_synapses)), :C => (prob_C, ones(num_synapses)))

probabilities = 0:0.01:1

println("Computing for single segment")
(net0,objects0) = load_network(YAML_source=config0)
firing_probabilities0 = [(prob_C .= a; P_fire(objects0[:C], volleys)) for a ∈ probabilities]

println("Computing for AND neuron")
(net1,objects1) = load_network(YAML_source=config1)
firing_probabilities1 = [(prob_A .= a; prob_B .= b; P_fire(objects1[:C], volleys)) for a ∈ probabilities, b ∈ probabilities]

println("Computing for OR neuron")
(net2,objects2) = load_network(YAML_source=config2)
firing_probabilities2 = [(prob_A .= a; prob_B .= b; P_fire(objects2[:C], volleys)) for a ∈ probabilities, b ∈ probabilities]

fig = Figure(resolution=(800,600))
ax0 = fig[1,1] = Axis(fig, title="$(num_synapses) synapses, θ = $(objects0[:C].θ_syn)", aspect=DataAspect(), xlabel="Trans. probability", ylabel="Firing probability")
ax1a = fig[1,2][2,1] = Axis(fig, aspect=DataAspect(), xlabel="P₁", ylabel="P₂")
ax1b = fig[1,2][2,2] = Axis(fig, aspect=DataAspect(), title="Θ=2")
ax2a = fig[1,2][3,1] = Axis(fig, aspect=DataAspect(), xlabel="P₁", ylabel="P₂")
ax2b = fig[1,2][3,2] = Axis(fig, aspect=DataAspect(), title="Θ=1")

lines!(ax0, probabilities, firing_probabilities0, colormap=:viridis, colorrange=(0,1), linewidth=3)

contourf!(ax1a, probabilities, probabilities, firing_probabilities1, colormap=:viridis, colorrange=(0,1))
p1=plot!(ax1b, objects1[:C], ports=Dict(:C=>[:C],:A=>[:A],:B=>[:B]), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:C=>:gray, :A=>color_1, :B=>color_2))
ports1 = lift(x->last.(x) .- Point2f0[(-0.4,0),(0.4,0)], p1.attributes[:ports])
arrows!(ax1b, ports1 , Node([Point2f0(-0.35,0),Point2f0(0.35,0)]), linewidth=2, arrowsize = [-20,20])
text!.(ax1b, "P₁", position=(@lift $ports1[1]), align=(:right, :bottom), textsize=0.25, color=:black)
text!.(ax1b, "P₂", position=(@lift $ports1[2]), align=(:left, :bottom), textsize=0.25, color=:black)

contourf!(ax2a, probabilities, probabilities, firing_probabilities2, colormap=:viridis, colorrange=(0,1))
p2=plot!(ax2b, objects2[:C], ports=Dict(:C=>[:C],:A=>[:A],:B=>[:B]), angle_between=20/180*π, branch_width=0.2, branch_length=1.0, color=Dict(:C=>:gray, :A=>color_1, :B=>color_2))
ports2 = lift(x->last.(x) .- Point2f0[(-0.4,0),(0.4,0)], p2.attributes[:ports])
arrows!(ax2b, ports2 , Node([Point2f0(-0.35,0),Point2f0(0.35,0)]), linewidth=2, arrowsize = [-20,20])
text!.(ax2b, "P₁", position=(@lift $ports2[1]), align=(:right, :bottom), textsize=0.25, color=:black)
text!.(ax2b, "P₂", position=(@lift $ports2[2]), align=(:left, :bottom), textsize=0.25, color=:black)

ax3 = Colorbar(fig[1,2][1,1], vertical=false, label="Firing Probability", height=10, limits=(0,1), colormap=:viridis)

hidedecorations!(ax1b)
hidedecorations!(ax2b)
hidespines!(ax1b)
hidespines!(ax2b)

hidexdecorations!(ax1a)
ax2a.xticks = [0,0.5,1.0]
ax1a.yticks=[0,0.5,1.0]
ax2a.yticks=[0,0.5,1.0]
ax3.ticks = [0,0.5,1.0]

ax0.xticks = [0,0.5,1.0]
ax0.yticks = [0,0.5,1.0]

ax0.limits[] = Rect(0.0f0,0.0f0,1.0f0,1.0f0)
ax1a.limits[] = Rect(0.0f0,0.0f0,1.0f0,1.0f0)
ax2a.limits[] = Rect(0.0f0,0.0f0,1.0f0,1.0f0)


save("probabilistic.svg", fig)