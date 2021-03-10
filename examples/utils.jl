using AbstractPlotting
import ColorTypes: RGB
import DataStructures: DefaultDict

pal = AbstractPlotting.Palette(:Dark2)
set_theme!(Theme(
    color = pal
))

color_1 = pal.colors[1]
color_1_50 = RGBAf0(color_1.r,color_1.g,color_1.b,0.5)
color_1_25 = RGBAf0(color_1.r,color_1.g,color_1.b,0.25)
color_2 = pal.colors[2]
color_2_50 = RGBAf0(color_2.r,color_2.g,color_2.b,0.5)
color_2_25 = RGBAf0(color_2.r,color_2.g,color_2.b,0.25)
color_3 = pal.colors[3]
color_3_50 = RGBAf0(color_3.r,color_3.g,color_3.b,0.5)
color_3_25 = RGBAf0(color_3.r,color_3.g,color_3.b,0.25)

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
        color=RGB(0.9,0.9,0.9),
        branch_width=0.1,
        branch_length=1.0,
        angle_between=10/180*pi,
        xticks = [],
        root_position = Point2f0(0,0),
        ports = Dict{Symbol,Vector{Symbol}}()
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

    
    function serialize_tree(tree::NeuronOrSegment{ID,T,WT,IT}, l, ω, given_ports) where {ID,T,WT,IT}
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
            α,d,points,ports,parents = serialize_tree(subtree, l, ω, given_ports)
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
            # shift by l and rotate by `branch_angle`
            c=cos(branch_angle)
            s=sin(branch_angle)
            transform = ((x,y),) -> Point2f0(c*x-s*(y+l),s*x+c*(y+l))
    
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

        return (sector=sector, depth=depth, nodes=all_points_flat, ports=all_ports_flat, parents=all_parents_flat)
    end



    # get all nodes, their ports and parents in flat format
    serialized = lift(treeplot[:branch_length],treeplot[:angle_between], treeplot[:ports]) do l,ω,prts
        ser=serialize_tree(tree,l,ω,prts)
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

        # draw circle to cap off branch
        poly!(treeplot, lift(w1,serialized, offset) do w1,ser,off
            node_dict = Dict(ser.nodes...)
            b2 = node_dict[name]
            Circle(off + b2, w1/2)
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


# """
# Computes diameters for branches to satisfy Rall's condition.
# """
# function get_branch_diameters(tree::Branch{I}; own_diam=1.0, diams=Dict{I,Float64}()) where {I}
#     diams[tree.name] = own_diam
#     child_diam = (own_diam^(3/2)/length(tree.branches))^(2/3)
#     for child in tree.branches
#         get_branch_diameters(child; own_diam=child_diam, diams=diams)
#     end
#     return diams
# end

