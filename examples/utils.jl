using AbstractPlotting
import ColorTypes: RGB
import DataStructures: DefaultDict

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
        show_port_labels=true,
        xticks = [],
        port_marker = x -> (marker=isa(x,Port{:pos}) ? "▲" : "o", color=RGB(0.9,0.9,0.9)),
        port_label = x -> (x.name, Point2f0(-0.1,-0.1), :textsize=>0.25, :align=>(:right,:center))
    )
end

function AbstractPlotting.plot!(treeplot::AbstractPlotting.Plot(Neuron{ID,T,WT,IT})) where {ID,T,WT,IT}
    # get the actual object to plot
    tree = to_value(treeplot[1])

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

    function serialize_tree(tree, l, ω)
        sectors = Float64[]
        depths = Float64[]
        all_points = Vector{Pair{ID,Point2f0}}[]
        all_points_flat = Pair{ID,Point2f0}[]
        all_ports_flat = Pair{ID,Vector{Port}}[tree.id=>[Port{syn.weight ≥ 0 ? :pos : :neg}(syn.id) for syn ∈ tree.synapses]]
        all_parents_flat = Pair{ID,ID}[]

        # go through all branches once to collect information
        for subtree in tree.next_upstream
            # get angle and depth of sector
            α,d,points,ports,parents = serialize_tree(subtree, l, ω)
            # compute depth of new sector
            ll = maybe_get(l, subtree.id)
            c = √(ll^2+d^2-2*ll*d*cos(π-α/2))
            # compute angle of new sector
            β=2acos((ll^2+c^2-d^2)/(2ll*c))

            push!(sectors, β+ maybe_get(ω, subtree.id))
            push!(depths, c)
            push!(all_points, points)
            append!(all_ports_flat, ports)
            append!(all_parents_flat, parents)
            push!(all_parents_flat, subtree.id=>tree.id)
        end

        # calculate total angle and depth for current (sub-)tree
        sector = sum(sectors)
        depth = isempty(depths) ? 0.0 : maximum(depths)

        # calculate angles for all branches
        branch_angles = cumsum(sectors) .- sectors./2 .- sector/2

        # go through all branches again and apply branch-specific transformations
        for (branch, branch_angle, points) in zip(tree.next_upstream, branch_angles, all_points)
            # shift by l and rotate by `branch_angle`
            c=cos(branch_angle)
            s=sin(branch_angle)
            transform = ((x,y),) -> (c*x-s*(y+l),s*x+c*(y+l))

            # add branch for branch
            push!(all_points_flat, branch.id=>transform(Point2f0(0.0,0.0)))
            
            # keep branches from branches' subtree
            if !isempty(points)            
                append!(all_points_flat, map(((name,pts),) -> name=>transform(pts), points))
            end
        end

        return (sector=sector, depth=depth, nodes=all_points_flat, ports=all_ports_flat, parents=all_parents_flat)
    end



    # get all nodes, their ports and parents in flat format
    serialized = lift(treeplot[:branch_length],treeplot[:angle_between]) do l,ω
        ser=serialize_tree(tree,l,ω)
        push!(ser.nodes, tree.name => Point2f0(0.0,0.0))
        ser
    end
    parent = DefaultDict(tree.name, serialized[].parents...)

    # plot branches
    for (name,parent_name) in serialized[].parents
        c  = maybe_get(treeplot[:color], name)
        w1 = maybe_get(treeplot[:branch_width], name)
        w2 = maybe_get(treeplot[:branch_width], parent_name)

        # dynamically recompute polygon
        branch_poly = lift(w1, w2, serialized) do w1,w2,ser
            node_dict = Dict(ser.nodes...)
            b1 = node_dict[name]
            b2 = node_dict[parent_name]
            normal = [0 1; -1 0]*(b2-b1)
            normal /= sqrt(normal'*normal)*2
            Point2f0[
                -normal*w1+b1,
                 normal*w1+b1,
                 normal*w2+b2,
                -normal*w2+b2
            ]
        end
        
        # draw polygon for the branch
        poly!(treeplot, branch_poly, color=c)

        # draw circle to cap off branch
        poly!(treeplot, lift(w2,serialized) do w2,ser
            node_dict = Dict(ser.nodes...)
            b2 = node_dict[name]
            Circle(b2, w2/2)
        end, color=c)
    end

    # draw polygon for the root
    c = maybe_get(treeplot[:color], tree.name)
    w = maybe_get(treeplot[:branch_width], tree.name)
    root_poly = lift(w) do w
        Point2f0[(-w/2*√(3), -w/2), (0,w), (w/2*√(3), -w/2)]
    end
    println(c)
    poly!(treeplot, root_poly, color=c)
   
    treeplot
end

