using AbstractPlotting

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
