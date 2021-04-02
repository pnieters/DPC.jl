using ADSP
using AbstractPlotting
import ColorTypes: RGB
import DataStructures: DefaultDict, OrderedDict
CairoMakie.activate!()


function make_manual_ticks(manual_ticks, manual_labels)
    @assert length(manual_ticks) == length(manual_labels)
    idx = sortperm(manual_ticks)
    t = copy(manual_ticks)
    l = copy(manual_labels)
    t[idx], x->l[idx]
end


pal = AbstractPlotting.Palette(:Dark2)

mytheme = Theme(
    fontsize = 18,
    font = "Linux Libertine O",
    Axis = (
        backgroundcolor = :gray90,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        xgridcolor = :white,
        ygridcolor = :white,
        titlesize = 20,
        xticklabelsize = 16,
        yticklabelsize = 16,
        xlabelsize = 18,
        ylabelsize = 18,
        xlabelpadding = 5f0,
        ylabelpadding = 5f0,
        xticklabelcolor = :gray10,
        yticklabelcolor = :gray10,
        xtickcolor = :gray10,
        ytickcolor = :gray10,
    )
)

set_theme!(mytheme)

presentation_theme = copy(mytheme)
presentation_theme[:Axis][:titlesize][] = 28
presentation_theme[:Axis][:xticklabelsize][] = 24
presentation_theme[:Axis][:yticklabelsize][] = 24
presentation_theme[:Axis][:xlabelsize][] = 26
presentation_theme[:Axis][:ylabelsize][] = 26
presentation_theme[:fontsize][] = 26


color_1 = pal.colors[1]
color_1_50 = RGBAf0(color_1.r,color_1.g,color_1.b,0.5)
color_1_25 = RGBAf0(color_1.r,color_1.g,color_1.b,0.25)
color_2 = pal.colors[2]
color_2_50 = RGBAf0(color_2.r,color_2.g,color_2.b,0.5)
color_2_25 = RGBAf0(color_2.r,color_2.g,color_2.b,0.25)
color_3 = pal.colors[3]
color_3_50 = RGBAf0(color_3.r,color_3.g,color_3.b,0.5)
color_3_25 = RGBAf0(color_3.r,color_3.g,color_3.b,0.25)
color_4 = pal.colors[4]
color_4_50 = RGBAf0(color_3.r,color_3.g,color_3.b,0.5)
color_4_25 = RGBAf0(color_3.r,color_3.g,color_3.b,0.25)

@recipe(Steps, t, y) do scene
    Theme(
        color = :black,
        fill = nothing,
    )
end

function AbstractPlotting.plot!(myplot::Steps)
    t = myplot[1]
    y = myplot[2]

    line = Node(Point2f0[])
    function update!(t, y)
        empty!(line[])
        y_old = y[1]
        for (_t,_y) in zip(t,y)
            push!(line[], Point2f0(_t,y_old))
            push!(line[], Point2f0(_t,_y))
            y_old = _y
        end
    end

    AbstractPlotting.Observables.onany(update!, t, y)
    update!(t[],y[])

    poly!(myplot, line, color=myplot[:fill])
    lines!(myplot, line, color=myplot[:color])
    myplot
end



struct Port{W,I}
    name::I
end
Port{W}(name::I) where {W,I} = Port{W,I}(name)


function AbstractPlotting.default_theme(scene::AbstractPlotting.SceneLike, ::Type{<: AbstractPlotting.Plot(Neuron)})
    Theme(
        color=RGBAf0(0.9,0.9,0.9,1.0),
        branch_width=0.1,
        branch_length=1.0,
        angle_between=10/180*pi,
        xticks = [],
        root_position = Point2f0(0,0),
        ports = Dict{Symbol,Vector{Symbol}}(),
        fxaa = false,
    )
end


function AbstractPlotting.plot!(treeplot::AbstractPlotting.Plot(Neuron))
    # get the actual object to plot
    tree = to_value(treeplot[1])
    # get position of the root = offset of all nodes
    offset = treeplot[:root_position]
    
    # if obj is an observable dict-type holding named values, get named value ...
    maybe_get(obj::Node, key) = if isa(obj[],Union{DefaultDict,Dict})
        @lift $obj[key]
    else
        obj
    end

    # if obj is a dict-type holding named values, get named value ...
    maybe_get(obj::Union{DefaultDict,Dict}, key) = obj[key]
    # ... otherwise return the obj itself
    maybe_get(obj, key) = obj

    
    function serialize_tree(tree::NeuronOrSegment{ID,T,WT,IT}, l, w, ω, given_ports) where {ID,T,WT,IT}
        sectors = Float64[]
        depths = Float64[]
        all_points = Vector{Pair{ID,Point2f0}}[]
        all_ports = Vector{Pair{ID,Point2f0}}[]
        all_points_flat = Pair{ID,Point2f0}[]
        all_ports_flat = Pair{ID,Point2f0}[]
        all_parents_flat = Pair{ID,ID}[]

        # go through all branches once to collect information
        for subtree in tree.next_upstream
            # get angle and depth of sector
            α,d,points,ports,parents = serialize_tree(subtree, l, w, ω, given_ports)
            # compute depth of new sector
            ll = maybe_get(l, subtree.id)
            c = √(ll^2+d^2-2*ll*d*cos(π-α/2))
            # compute angle of new sector
            β=2acos((ll^2+c^2-d^2)/(2ll*c))

            push!(sectors, β+ maybe_get(ω, subtree.id))
            push!(depths, c)
            push!(all_points, points)
            push!(all_ports, ports)
            append!(all_parents_flat, parents)
            push!(all_parents_flat, subtree.id=>tree.id)
        end

        # calculate total angle and depth for current (sub-)tree
        sector = sum(sectors)
        depth = isempty(depths) ? 0.0 : maximum(depths)

        # calculate angles for all branches
        branch_angles = cumsum(sectors) .- sectors./2 .- sector/2

        for (branch, branch_angle, points, ports) in zip(tree.next_upstream, branch_angles, all_points, all_ports)
            # shift by ll and rotate by `branch_angle`
            ll = maybe_get(l, branch.id)
            c=cos(branch_angle)
            s=sin(branch_angle)
            transform = ((x,y),) -> Point2f0(c*x-s*(y+ll),s*x+c*(y+ll))
    
            branch_end = transform(Point2f0(0.0,0.0))
            # add branch for branch
            push!(all_points_flat, branch.id=>branch_end)
    
            # add ports for branch
            branch_ports = get(given_ports, branch.id, ID[])
            append!(all_ports_flat, Pair.(branch_ports, LinRange(Point2f0(0,0), branch_end, length(branch_ports)+2)[2:end-1]))
    
            # keep branches from branches' subtree
            if !isempty(points)            
                append!(all_points_flat, map(((name,pts),) -> name=>transform(pts), points))
            end
    
            # keep ports from branches' subtree
            if !isempty(ports)            
                append!(all_ports_flat, map(((prt,pos),) -> prt=>transform(pos), ports))
            end
        end

        if isa(tree, Neuron)
            neuron_ports = get(given_ports, tree.id, ID[])

            ww = maybe_get(w, tree.id)
            append!(all_ports_flat, Pair.(neuron_ports, LinRange(Point2f0(0,-ww), Point2f0(0,ww), length(neuron_ports)+2)[2:end-1]))
        end

        return (sector=sector, depth=depth, nodes=all_points_flat, ports=all_ports_flat, parents=all_parents_flat)
    end



    # get all nodes, their ports and parents in flat format
    serialized = lift(treeplot[:branch_length],treeplot[:branch_width],treeplot[:angle_between], treeplot[:ports]) do l,w,ω,prts
        ser=serialize_tree(tree,l,w,ω,prts)
        push!(ser.nodes, tree.id => Point2f0(0.0,0.0))
        ser
    end
    parent = DefaultDict(tree.id, serialized[].parents...)
    # get all the positions of the ports
    ports = lift(serialized,offset) do ser,off
        [k=>v+off for (k,v) in ser.ports]
    end

    # plot branches
    for (name,parent_name) in serialized[].parents
        c  = maybe_get(treeplot[:color], name)
        w1 = maybe_get(treeplot[:branch_width], name)
        w2 = maybe_get(treeplot[:branch_width], parent_name)

        # dynamically recompute polygon
        branch_poly = lift(w1, w2, serialized, offset) do w1,w2,ser,off
            node_dict = Dict(ser.nodes...)
            b1 = node_dict[name]
            b2 = node_dict[parent_name]
            normal = [0 1; -1 0]*(b2-b1)
            normal /= sqrt(normal'*normal)*2
            Ref(off) .+ Point2f0[
                -normal*w1+b1,
                 normal*w1+b1,
                 normal*w2+b2,
                -normal*w2+b2
            ]
        end
        
        # draw polygon for the branch
        poly!(treeplot, branch_poly, color=c, strokewidth = 1, strokecolor=c)

        # draw circles to cap off branch
        poly!(treeplot, lift(w1,serialized, offset) do w1,ser,off
            node_dict = Dict(ser.nodes...)
            b2 = node_dict[name]
            decompose(Point2f0, Circle(off + b2, w1/2))
        end, color=c, strokewidth = 1, strokecolor=c)

        poly!(treeplot, lift(w1,offset) do w1,off
            decompose(Point2f0, Circle(off, w1/2))
        end, color=c, strokewidth = 1, strokecolor=c)
    end

    # draw polygon for the root
    c = maybe_get(treeplot[:color], tree.id)
    w = maybe_get(treeplot[:branch_width], tree.id)
    root_poly = lift(w,offset) do w,offset
        Point2f0[(-2w/2*√(3), -2w/2), (0,2w), (2w/2*√(3), -2w/2)].+offset
    end
    poly!(treeplot, root_poly, color=c, strokewidth = 1, strokecolor=c)

    treeplot.attributes[:ports] = ports
    treeplot
end


"""
Computes diameters for branches to satisfy Rall's condition.
"""
function get_branch_diameters(tree::NeuronOrSegment{ID,T,WT,IT}; root_diam=1.0, diams=Dict{ID,Float64}()) where {ID,T,WT,IT}
    diams[tree.id] = root_diam
    child_diam = (root_diam^(3/2)/length(tree.next_upstream))^(2/3)
    for child in tree.next_upstream
        get_branch_diameters(child; root_diam=child_diam, diams=diams)
    end
    return diams
end



"""
Computes the depths in segments from the soma for each segment.
"""
function get_root_distances(tree::NeuronOrSegment{ID,T,WT,IT}; start=0, dists = OrderedDict{ID,Int}()) where {ID,T,WT,IT}
    dists[tree.id] = start
    for child in tree.next_upstream
        get_root_distances(child; start=start+1, dists=dists)
    end
    return dists
end


"""
Draw a grid of pie-plots to illustrate a discrete spatial distribution of class probabilities.
"""
function plot_pie_grid(ax, xgrid, ygrid, class_counts...; 
    class_colors=[RGBAf0(c,c,c) for c in LinRange(0,1,length(class_counts))], 
    color=repeat(class_colors, 1, length(xgrid), length(ygrid)))

    scale = min(minimum(diff(xgrid)),minimum(diff(ygrid)))/(2*√(maximum(sum(class_counts))))
    cnt_tmps = zeros(Int,length(class_counts))
    for (i,x) in  enumerate(xgrid), (j,y) in enumerate(ygrid)
        for (k,class) in enumerate(class_counts)
            cnt_tmps[k] = class[i,j]
        end

        frac_pos = cumsum(cnt_tmps) ./ sum(cnt_tmps)
        Φs = 2π .* frac_pos
        r = √(sum(cnt_tmps))*scale
        if sum(cnt_tmps) > 0
            Φ₀ = 0
            for (k,Φ) in enumerate(Φs)
                arc = Point2f0[[Point2f0(x,y)]; [Point2f0(x+r*sin(ϕ),y+r*cos(ϕ)) for ϕ in LinRange(Φ₀,Φ,round(Int,Φ*100))]; [Point2f0(x,y)]]
                
                poly!(ax, arc, color=color[k,i,j])
                Φ₀ = Φ
            end
        end
    end    
end