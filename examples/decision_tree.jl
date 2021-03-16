using ADSP, CairoMakie, Random, Optim
include("utils.jl")

## Network configurations
config_preamble = """
refractory_duration: 100.5
inputs:
  - id: X1
  - id: X2
  - id: X3
  - id: X4
  - id: X5
  - id: X6
  - id: X7
  - id: X8
  - id: X9
  - id: X10
  - id: Y1
  - id: Y2
  - id: Y3
  - id: Y4
  - id: Y5
  - id: Y6
  - id: Y7
  - id: Y8
  - id: Y9
  - id: Y10
  - id: bias
synapses:
  - {id: synXn1, source: X1, target: n, delay: 41}
  - {id: synXn2, source: X2, target: n, delay: 41}
  - {id: synXn3, source: X3, target: n, delay: 41}
  - {id: synXn4, source: X4, target: n, delay: 41}
  - {id: synXn5, source: X5, target: n, delay: 41}
  - {id: synXn6, source: X6, target: n, delay: 41}
  - {id: synXn7, source: X7, target: n, delay: 41}
  - {id: synXn8, source: X8, target: n, delay: 41}
  - {id: synXn9, source: X9, target: n, delay: 41}
  - {id: synXn10, source: X10, target: n, delay: 41}
  - {id: synYn1, source: Y1, target: n, delay: 41}
  - {id: synYn2, source: Y2, target: n, delay: 41}
  - {id: synYn3, source: Y3, target: n, delay: 41}
  - {id: synYn4, source: Y4, target: n, delay: 41}
  - {id: synYn5, source: Y5, target: n, delay: 41}
  - {id: synYn6, source: Y6, target: n, delay: 41}
  - {id: synYn7, source: Y7, target: n, delay: 41}
  - {id: synYn8, source: Y8, target: n, delay: 41}
  - {id: synYn9, source: Y9, target: n, delay: 41}
  - {id: synYn10, source: Y10, target: n, delay: 41}
  - {id: synbiasn, source: bias, target: n, delay: 41}
"""

config_shallow = config_preamble * """
neurons:
  - id: n
    θ_syn: 1.0
"""

config_deep = config_preamble * """
  - {id: synX31, source: X1, target: seg3, delay: 1}
  - {id: synX32, source: X2, target: seg3, delay: 1}
  - {id: synX33, source: X3, target: seg3, delay: 1}
  - {id: synX34, source: X4, target: seg3, delay: 1}
  - {id: synX35, source: X5, target: seg3, delay: 1}
  - {id: synX36, source: X6, target: seg3, delay: 1}
  - {id: synX37, source: X7, target: seg3, delay: 1}
  - {id: synX38, source: X8, target: seg3, delay: 1}
  - {id: synX39, source: X9, target: seg3, delay: 1}
  - {id: synX310, source: X10, target: seg3, delay: 1}
  - {id: synY31, source: Y1, target: seg3, delay: 1}
  - {id: synY32, source: Y2, target: seg3, delay: 1}
  - {id: synY33, source: Y3, target: seg3, delay: 1}
  - {id: synY34, source: Y4, target: seg3, delay: 1}
  - {id: synY35, source: Y5, target: seg3, delay: 1}
  - {id: synY36, source: Y6, target: seg3, delay: 1}
  - {id: synY37, source: Y7, target: seg3, delay: 1}
  - {id: synY38, source: Y8, target: seg3, delay: 1}
  - {id: synY39, source: Y9, target: seg3, delay: 1}
  - {id: synY310, source: Y10, target: seg3, delay: 1}
  - {id: synbias3, source: bias, target: seg3, delay: 1}
  - {id: synX21, source: X1, target: seg2, delay: 11}
  - {id: synX22, source: X2, target: seg2, delay: 11}
  - {id: synX23, source: X3, target: seg2, delay: 11}
  - {id: synX24, source: X4, target: seg2, delay: 11}
  - {id: synX25, source: X5, target: seg2, delay: 11}
  - {id: synX26, source: X6, target: seg2, delay: 11}
  - {id: synX27, source: X7, target: seg2, delay: 11}
  - {id: synX28, source: X8, target: seg2, delay: 11}
  - {id: synX29, source: X9, target: seg2, delay: 11}
  - {id: synX210, source: X10, target: seg2, delay: 11}
  - {id: synY21, source: Y1, target: seg2, delay: 11}
  - {id: synY22, source: Y2, target: seg2, delay: 11}
  - {id: synY23, source: Y3, target: seg2, delay: 11}
  - {id: synY24, source: Y4, target: seg2, delay: 11}
  - {id: synY25, source: Y5, target: seg2, delay: 11}
  - {id: synY26, source: Y6, target: seg2, delay: 11}
  - {id: synY27, source: Y7, target: seg2, delay: 11}
  - {id: synY28, source: Y8, target: seg2, delay: 11}
  - {id: synY29, source: Y9, target: seg2, delay: 11}
  - {id: synY210, source: Y10, target: seg2, delay: 11}
  - {id: synbias2, source: bias, target: seg2, delay: 11}
  - {id: synX11, source: X1, target: seg1, delay: 21}
  - {id: synX12, source: X2, target: seg1, delay: 21}
  - {id: synX13, source: X3, target: seg1, delay: 21}
  - {id: synX14, source: X4, target: seg1, delay: 21}
  - {id: synX15, source: X5, target: seg1, delay: 21}
  - {id: synX16, source: X6, target: seg1, delay: 21}
  - {id: synX17, source: X7, target: seg1, delay: 21}
  - {id: synX18, source: X8, target: seg1, delay: 21}
  - {id: synX19, source: X9, target: seg1, delay: 21}
  - {id: synX110, source: X10, target: seg1, delay: 21}
  - {id: synY11, source: Y1, target: seg1, delay: 21}
  - {id: synY12, source: Y2, target: seg1, delay: 21}
  - {id: synY13, source: Y3, target: seg1, delay: 21}
  - {id: synY14, source: Y4, target: seg1, delay: 21}
  - {id: synY15, source: Y5, target: seg1, delay: 21}
  - {id: synY16, source: Y6, target: seg1, delay: 21}
  - {id: synY17, source: Y7, target: seg1, delay: 21}
  - {id: synY18, source: Y8, target: seg1, delay: 21}
  - {id: synY19, source: Y9, target: seg1, delay: 21}
  - {id: synY110, source: Y10, target: seg1, delay: 21}
  - {id: synbias1, source: bias, target: seg1, delay: 21}
neurons:
  - id: n
    θ_syn: 1.0
    branches:
      - id: seg1
        θ_syn: 1.0
        branches:
          - id: seg2
            θ_syn: 1.0
            branches:
              - id: seg3
                θ_syn: 1.0
"""

function check_for_activation(net,objects,obj_id, stimulus; num_synapses=10)
    inputs = [[Event(:input_spikes, 0.0, 1.0, objects[:bias])];[Event(:input_spikes, 0.0, 1.0, objects[Symbol(XY,i)]) for (XY,S) in pairs(stimulus) for i in randperm(num_synapses)[1:S]]]
    logger = simulate!(net, inputs, 150.0)
    spikes = filter(x->x.object==obj_id && (x.event==:spikes || x.event==:plateau_starts), logger.data)
    !isempty(spikes)
end

function find_optimal_split(xgrid, ygrid, counts_pos, counts_neg)
    all = counts_pos .+ counts_neg

    """
    Apply a soft split with offset `off` and parameters `wx`, `wy` to the `data`.
    Returns the respective probabilities `s1_prob,s2_prob` for potential data points to land on either side.
    """
    function soft_split(xgrid,ygrid, off, wx, wy; slope=1)
        lin = @. xgrid*wx + ygrid*wy + off-1
        s1_prob = @. 1.0/(1.0+exp(-slope*lin))
        return (s1_prob,1 .- s1_prob)
    end

    # initial parameters (each is a view into the parameter vector)
    x₀,y₀ = 5.5 .+ randn(2)*donut_hole_size/2
    α = rand()*2π
    c = 1.0
    
    function loss(my_params)
        p_pos_unsplit = sum(counts_pos)/sum(counts_pos.+counts_neg)
        info_unsplit = -p_pos_unsplit*log2(p_pos_unsplit)-(1-p_pos_unsplit)*log2(1-p_pos_unsplit)

        off,α = my_params
        (s1_prob,s2_prob) = soft_split(xgrid, ygrid', off, sin(α),cos(α))
        p_side(s) = s ? sum(s1_prob .* all) / sum(all) : 1-p_side(~s)
        p_class(c) = c ? sum(counts_pos) / sum(all) : 1-p_class(~c)
        p_side_class(s,c) = if (s,c) == (true,true)
            sum(s1_prob .* counts_pos) / sum(counts_pos)
        elseif (s,c) == (true,false)
            sum(s1_prob .* counts_neg) / sum(counts_neg)
        else
            1-p_side_class(~s,c)
        end
        p_class_and_side(c,s) = p_side_class(s,c)*p_class(c)
        p_class_side(c,s) = p_class_and_side(c,s)/p_side(s)
        info(s) = sum(c -> -p_class_side(c,s)*log2(p_class_side(c,s)), (true,false))
        info_split = sum(s->p_side(s)*info(s), (true,false))
        
        # flip if the wrong side is more likely to contain the desired class
        flip = p_side_class(false,true) > p_side_class(true,false)
        return (loss=info_split-info_unsplit, flip=flip)
    end

    opt_params = optimize(x->loss(x).loss, [c,α], NelderMead()).minimizer
    off,wx,wy = opt_params[1],sin(opt_params[2]),cos(opt_params[2])
    
    if loss(opt_params).flip
        # x*wx+y*wy+c > 1 ==invert=> x*wx+y*wy+c < 1 ==> -x*wx+-y*wy-c > -1 ==> -x*wx-y*wy-c+2 > 1
        println("Flipped!")
        off,wx,wy = (-off+2,-wx,-wy) 
    end

    return off,wx,wy
end

function set_parameters(obj, offset, wx, wy)
    for input in obj.net.inputs
        w = if startswith(String(input.id),"X") 
            wx 
        elseif startswith(String(input.id),"Y")
            wy
        elseif input.id == :bias
            offset
        else
            error("Don't know what to do with $(input.id)")
        end

        for synapse in input.synapses
            if synapse.target == obj
                synapse.weight = w
                # hack: make sure the inhibitory spikes arrive first and lasts longer, otherwise we could spike before inhibition
                synapse.delay += synapse.weight > 0 ? 1.0 : 0.0
                synapse.spike_duration = synapse.weight > 0 ? 5.0 : 7.0
            end
        end
    end
end

function plot_classifier!(ax, off, wx, wy; kwargs...)
    lines!(ax, [Point2f0(-1,(1-off+wx)/wy), Point2f0(11,(1-off-wx*11)/wy)]; kwargs...)
end

## Generate data

num_samples = 500
donut_hole_size = 4
donut_width = 0.5
center = [5.5, 5.5]

xgrid = 0:10
ygrid = 0:10

quantize(ar, grid) = clamp.(searchsortedlast.(Ref((grid[2:end] .+ grid[1:end-1])/2), ar), 0, 10)

pos_real = center .+ donut_hole_size/2.5*randn(2, num_samples)
neg_real = center .+ [fun(ϕ)*r for fun ∈ (sin, cos), (ϕ,r) ∈ zip(2π.*rand(num_samples), donut_hole_size .+ donut_width*randn(num_samples))]

pos = [quantize(pos_real[1,:], xgrid)'; quantize(pos_real[2,:], ygrid)']
neg = [quantize(neg_real[1,:], xgrid)'; quantize(neg_real[2,:], ygrid)']

counts_pos = zeros(Int, length(xgrid), length(ygrid))
counts_neg = zeros(Int, length(xgrid), length(ygrid))

for i in 1:num_samples
    counts_pos[1 .+ pos[:,i]...] += 1
    counts_neg[1 .+ neg[:,i]...] += 1
end

## Compute optimal splits

params = Tuple[]
counts_pos_tmp = copy(counts_pos)
counts_neg_tmp = copy(counts_neg)
effective = fill(0, size(counts_pos)...)
effective_shallow = fill(0, size(counts_pos)...)

# shallow neuron
(net_shallow,objects_shallow) = load_network(;YAML_source=config_shallow, weight_type=Float32, synaptic_input_type=Float32)
off,wx,wy = find_optimal_split(xgrid, ygrid, counts_pos_tmp, counts_neg_tmp)
params_shallow = (off,wx,wy)
set_parameters(objects_shallow[:n], off, wx, wy)

for (i,x) in enumerate(xgrid), (j,y) in enumerate(ygrid)
    did_fire = check_for_activation(net_shallow, objects_shallow, :n, (X=x,Y=y); num_synapses=10)
    if did_fire
        effective_shallow[i,j] = 1
    end
end

# deep neuron
(net,objects) = load_network(;YAML_source=config_deep, weight_type=Float32, synaptic_input_type=Float32)
for (level,obj_id) in zip(1:4,[:seg3,:seg2,:seg1,:n])
    off,wx,wy = find_optimal_split(xgrid, ygrid, counts_pos_tmp, counts_neg_tmp)
    push!(params, (off,wx,wy))
    set_parameters(objects[obj_id], off, wx, wy)

    for (i,x) in enumerate(xgrid), (j,y) in enumerate(ygrid)
        did_fire = check_for_activation(net, objects, obj_id, (X=x,Y=y); num_synapses=10)
        if ~did_fire
            counts_pos_tmp[i,j] = 0
            counts_neg_tmp[i,j] = 0
        else
            effective[i,j] = level
        end
    end
end

function select_colors(effective, level)
    class_colors_full=[color_1, color_2]
    class_colors_shaded=[color_1_25, color_2_25]
    colors = zeros(RGBAf0,2,length(xgrid),length(ygrid))
    for (i,x) in enumerate(xgrid), (j,y) in enumerate(ygrid)
        colors[:, i, j] = if effective[i,j] ≥ level
            class_colors_full
        else 
            class_colors_shaded
        end
    end
    return colors
end


## Figure 1: show raw data
fig1 = Figure(resolution=(700,700))
ax = fig1[1,1] = Axis(fig1, xticks=xgrid, yticks=ygrid, 
    title="Example: 2D classification - raw data", 
    xlabel="feature encoded by population X",
    ylabel="feature encoded by populytion Y")
scatter!(ax, pos_real, marker=:x, strokecolor=color_1)
scatter!(ax, neg_real, marker=:x, strokecolor=color_2)

ylims!(ax, -0.5,10.5)
xlims!(ax, -0.5,10.5)
save(joinpath("figures","classifier_walkthrough","fig1.png"), fig1)
fig1

## Figure 2: show encoding
fig2 = Figure(resolution=(700,700))
ax = fig2[1,1] = Axis(fig2, xticks=xgrid, yticks=ygrid, 
    title="Example: encoding into a spike volley", 
    xlabel="#spikes from population X",
    ylabel="#spikes from population Y")

colors = select_colors(effective, 0)
plot_pie_grid(ax,xgrid,ygrid, counts_pos, counts_neg, color=colors)
    
ylims!(ax, -0.5,10.5)
xlims!(ax, -0.5,10.5)
save(joinpath("figures","classifier_walkthrough","fig2.png"), fig2)
fig2



## Figure 3:
fig3 = Figure(resolution=(800,700))
ax = fig3[1,1] = Axis(fig3, xticks=xgrid, yticks=ygrid, 
    title="Example: What a linear classifier can see & do", 
    xlabel="#spikes from population X",
    ylabel="#spikes from population Y")


ax_schema = fig3[1,2] = Axis(fig3, width=100, aspect=DataAspect(), backgroundcolor=:transparent)
plt=plot!(ax_schema, objects_shallow[:n], branch_width=1.0, branch_length=5.0, 
    color=DefaultDict(to_color(:gray40), neuron_colors[end-level+1:end]...),
    ports=Dict(:n=>[:nX,:nY]))

ports = Dict(plt.attributes[:ports][])

colors = select_colors(effective_shallow, 1)
plot_pie_grid(ax,xgrid,ygrid, counts_pos, counts_neg, color=colors)
plot_classifier!(ax, params_shallow...; linewidth=3.0)

for XY in [:X,:Y]
    pos = ports[Symbol("n"*String(XY))]
    arrows!(ax_schema, [pos+Point2f0(-2.0,0)], [Point2f0(1.25,0)], linewidth=2, arrowsize = 20)
    text!.(ax_schema, String(XY), position=pos+Point2f0(-2.0,0), align=(:right, :center), textsize=0.75, color=:black)
end

hidedecorations!(ax_schema)
hidespines!(ax_schema)
ylims!(ax, -0.5,10.5)
xlims!(ax, -0.5,10.5)
save(joinpath("figures","classifier_walkthrough","fig3.png"), fig3)
fig3

## Figures 4-6:
cols = cgrad(:viridis, 5, categorical=true, rev=true)
dists = get_root_distances(objects[:n])
neuron_colors = [k=>cols[2+v] for (k,v) in pairs(dists)]

for (level,name) in zip(0:4, ["first segment", "first segment", "second segment", "third segment", "soma"])
    figᵢ = Figure(resolution=(800,700))
    ax = figᵢ[1,1] = Axis(figᵢ, xticks=xgrid, yticks=ygrid, 
        title="Example: What the $(name) can see & do", 
        xlabel="#spikes from population X",
        ylabel="#spikes from population Y")

    ax_schema = figᵢ[1,2] = Axis(figᵢ, width=100, aspect=DataAspect(), backgroundcolor=:transparent)
    plt=plot!(ax_schema, objects[:n], branch_width=1.0, branch_length=5.0, 
        color=DefaultDict(to_color(:gray40), neuron_colors[end-level+1:end]...),
        ports=Dict(:n=>[:nX,:nY],:seg1=>[:seg1X,:seg1Y],:seg2=>[:seg2X,:seg2Y],:seg3=>[:seg3X,:seg3Y]))
    
    ports = Dict(plt.attributes[:ports][])
    for seg in [:n, :seg1, :seg2, :seg3], XY in [:X,:Y]
        pos = ports[Symbol(String(seg)*String(XY))]
        arrows!(ax_schema, [pos+Point2f0(-2.0,0)], [Point2f0(1.25,0)], linewidth=2, arrowsize = 20)
        text!.(ax_schema, String(XY), position=pos+Point2f0(-2.0,0), align=(:right, :center), textsize=0.75, color=:black)
    end

    colors = select_colors(effective, level)
    plot_pie_grid(ax,xgrid,ygrid, counts_pos, counts_neg, color=colors)
    for (i,p) in enumerate(params[1:level])
        plot_classifier!(ax, p...; linewidth= i==level ? 3.0 : 1.0)
    end

    hidedecorations!(ax_schema)
    hidespines!(ax_schema)
    ylims!(ax, -0.5,10.5)
    xlims!(ax, -0.5,10.5)
    display(figᵢ)
    save(joinpath("figures","classifier_walkthrough","fig$(level+4).png"), figᵢ)
end
