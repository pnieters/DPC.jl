using ADSP, JLD2, Distributions, DataFrames, LaTeXStrings, Plots, KernelDensity

println("Loading data...")
@load "examples/rate_coding/data/examples_data.jld2" input_rates input_rates2d output_rates output_rates2 output_rates3 output_rates4 output_rates5 

r_A = 35.0

function interpolate(inp,out)
    function interpolated(r)
        i=searchsortedlast(inp,r)
        return if i ≤ 0
            out[1]
        elseif i ≥ length(inp)
            out[end]
        else
            c = (r-inp[i])/(inp[i+1]-inp[i])
            c*out[i+1]+(1-c)*out[i]
        end
    end
end

in_out = interpolate(input_rates, output_rates2)
out_rates2d = in_out.(input_rates2d)
ideal_and = out_rates2d.*0.1 * out_rates2d'*0.1 * in_out(r_A)
ideal_or = (1 .-(1 .-out_rates2d.*0.1) .* (1 .-out_rates2d'.*0.1)) * in_out(r_A)

"""
| schema | I/O curve |        |   AND   |        |   OR    |
| schema | contour   | schema | contour | schema | contour |
"""

schema_soma = plot(framestyle=:none)
schema_seq = plot(framestyle=:none)
schema_and = plot(framestyle=:none)
schema_or  = plot(framestyle=:none)

filler = plot(framestyle=:none)

io_curve = plot(legend=false, title="input/output", xlabel=L"r_A", ylabel=L"\varrho(r_A)")
plot!(input_rates, output_rates2, color=:black, linewidth=2)
# plot!(twinx(), input_rates, output_rates, color=:gray, linestyle=:dash, linewidth=2, legend=false)

contour_seq = plot(aspect_ratio=:equal, framesyle=:origin, colorbar=false,xlims=(0,50),ylims=(0,50), title="'then' neuron",xlabel=L"r_B", ylabel=L"r_C")
contourf!(input_rates2d,input_rates2d,output_rates3, levels=10)
contour_or  = plot(aspect_ratio=:equal, framesyle=:origin, colorbar=false,xlims=(0,50),ylims=(0,50), title="'or' neuron",xlabel=L"r_B", ylabel=L"r_C")
contourf!(input_rates2d,input_rates2d,output_rates4, levels=10)
contour_and = plot(aspect_ratio=:equal, framesyle=:origin, colorbar=false,xlims=(0,50),ylims=(0,50), title="'and' neuron",xlabel=L"r_B", ylabel=L"r_C")
contourf!(input_rates2d,input_rates2d,output_rates5, levels=10)

contour_ideal_and = plot(aspect_ratio=:equal, framesyle=:origin, colorbar=false,xlims=(0,50),ylims=(0,50), title="ideal 'and'",xlabel=L"r_B", ylabel=L"r_C")
contourf!(input_rates2d, input_rates2d, ideal_and,levels=10)
contour_ideal_or = plot(aspect_ratio=:equal, framesyle=:origin, colorbar=false,xlims=(0,50),ylims=(0,50), title="ideal 'or'",xlabel=L"r_B", ylabel=L"r_C")
contourf!(input_rates2d, input_rates2d, ideal_or, levels=10)


plts = [
    schema_soma io_curve    filler      contour_ideal_and   filler      contour_ideal_or;
    schema_seq  contour_seq schema_and  contour_and         schema_or   contour_or
]

widths = [1.0,4.0,1.0,4.0,1.0,4.0]
widths./=sum(widths)

w = 4.7747 * 150
h = 0.45*w
plot(permutedims(plts)..., layout=grid(2,6, widths=widths), size=(w,h))
savefig("examples/rate_coding/figures/functions.svg")
