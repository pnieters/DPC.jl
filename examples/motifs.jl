using ADSP, CairoMakie
#include(joinpath(@__DIR__, "utils.jl"))
include("utils.jl")

config_and = """
refractory_duration: 5.01
inputs:
- id: i11
- id: i12
- id: i13
- id: i14
- id: i15
- id: i21
- id: i22
- id: i23
- id: i24
- id: i25
- id: i31
- id: i32
- id: i33
- id: i34
- id: i35
neurons:
- id: n
  θ_syn: 5
  θ_seg: 1
  branches: 
    - id: seg0
      θ_syn: 5
      θ_seg: 2
      branches: 
          - id: seg1
            θ_syn: 5
          - id: seg2
            θ_syn: 5
synapses:
- {id: syn11, source: i11, target: seg1}
- {id: syn12, source: i12, target: seg1}
- {id: syn13, source: i13, target: seg1}
- {id: syn14, source: i14, target: seg1}
- {id: syn15, source: i15, target: seg1}
- {id: syn21, source: i21, target: seg2}
- {id: syn22, source: i22, target: seg2}
- {id: syn23, source: i23, target: seg2}
- {id: syn24, source: i24, target: seg2}
- {id: syn25, source: i25, target: seg2}
- {id: syn31, source: i31, target: seg0}
- {id: syn32, source: i32, target: seg0}
- {id: syn33, source: i33, target: seg0}
- {id: syn34, source: i34, target: seg0}
- {id: syn35, source: i35, target: seg0}
"""


config_or = """
refractory_duration: 5.01
inputs:
- id: i11
- id: i12
- id: i13
- id: i14
- id: i15
- id: i21
- id: i22
- id: i23
- id: i24
- id: i25
- id: i31
- id: i32
- id: i33
- id: i34
- id: i35
neurons:
- id: n
  θ_syn: 0
  θ_seg: 1
  branches: 
    - id: seg0
      θ_syn: 5
      θ_seg: 1
      branches: 
        - id: seg1
          θ_syn: 5
        - id: seg2
          θ_syn: 5
synapses:
- {id: syn11, source: i11, target: seg1}
- {id: syn12, source: i12, target: seg1}
- {id: syn13, source: i13, target: seg1}
- {id: syn14, source: i14, target: seg1}
- {id: syn15, source: i15, target: seg1}
- {id: syn21, source: i21, target: seg2}
- {id: syn22, source: i22, target: seg2}
- {id: syn23, source: i23, target: seg2}
- {id: syn24, source: i24, target: seg2}
- {id: syn25, source: i25, target: seg2}
- {id: syn31, source: i31, target: seg0}
- {id: syn32, source: i32, target: seg0}
- {id: syn33, source: i33, target: seg0}
- {id: syn34, source: i34, target: seg0}
- {id: syn35, source: i35, target: seg0}
"""


config_chain = """
refractory_duration: 5.01
inputs:
- id: i11
- id: i12
- id: i13
- id: i14
- id: i15
- id: i21
- id: i22
- id: i23
- id: i24
- id: i25
- id: i31
- id: i32
- id: i33
- id: i34
- id: i35
neurons:
- id: n
  θ_syn: 0
  θ_seg: 1
  branches: 
    - id: seg0
      θ_syn: 5
      θ_seg: 1
      branches: 
        - id: seg2
          θ_syn: 5
          branches: 
            - id: seg1
              θ_syn: 5
synapses:
- {id: syn11, source: i11, target: seg1}
- {id: syn12, source: i12, target: seg1}
- {id: syn13, source: i13, target: seg1}
- {id: syn14, source: i14, target: seg1}
- {id: syn15, source: i15, target: seg1}
- {id: syn21, source: i21, target: seg2}
- {id: syn22, source: i22, target: seg2}
- {id: syn23, source: i23, target: seg2}
- {id: syn24, source: i24, target: seg2}
- {id: syn25, source: i25, target: seg2}
- {id: syn31, source: i31, target: seg0}
- {id: syn32, source: i32, target: seg0}
- {id: syn33, source: i33, target: seg0}
- {id: syn34, source: i34, target: seg0}
- {id: syn35, source: i35, target: seg0}
"""

volleys = [(15.0, :i1), (75.0, :i2), (140.0, :i3), (215.0, :i2), (275.0, :i1), (300.0, :i3), (415.0, :i2), (475.0, :i3)]
spikes = [ (time=t+τ-0.5, name=Symbol("$(pop)$(i)")) for (t, pop) in volleys for (i,τ) in enumerate(5*rand(5))]

## Set up plot

fig = Figure(resolution = (0.75textwidth, 0.75textwidth))
ax_left  = fig[1,1] = Axis(fig; backgroundcolor=:transparent)
ax_right = fig[1,3] = Axis(fig; backgroundcolor=:transparent)
gl = fig[1,2] = GridLayout(alignmode=Outside())
ax_chain = gl[1,1] = Axis(fig; 
    yticks=(2.5:5:15, ["C","B","A"]), 
    yminorticks=AbstractPlotting.IntervalsBetween(5, true),
    yminorticksvisible = true, yminorgridvisible = true, 
    titlealign=:left, title = "a.    Sequential segments"
)
ax_or    = gl[2,1] = Axis(fig; 
    yticks=(2.5:5:15, ["C","B","A"]), 
    yminorticks=AbstractPlotting.IntervalsBetween(5, true),
    yminorticksvisible = true, yminorgridvisible = true, 
    titlealign=:left, title = "b.    Parallel segments (either required)"
)
ax_and   = gl[3,1] = Axis(fig; 
    yticks=(2.5:5:15, ["C","B","A"]),  
    yminorticks=AbstractPlotting.IntervalsBetween(5, true),
    yminorticksvisible = true, yminorgridvisible = true, 
    titlealign=:left, title = "c.    Parallel segments (both required)", xlabel="time [ms]"
)

linkaxes!(ax_chain,ax_and)
linkaxes!(ax_chain,ax_or)

hidexdecorations!(ax_chain, grid=false)
hidexdecorations!(ax_or, grid=false)
hidedecorations!(ax_left)
hidedecorations!(ax_right)

fig

## Simulate each motive
offset = Dict(:chain => Point2f0(0,-0.5), :and => Point2f0(0,-1.35) , :or => Point2f0(0,-2.7))

for (name, config, ax) in ((:chain,config_chain,ax_chain), (:or,config_or,ax_or), (:and,config_and,ax_and))
    (net,objects) = load_network(YAML_source=config)
    input=[ Event(:input_spikes, 0.0, t, objects[name]) for (t, name) in spikes]
    logger=simulate!(net, input)

    # get traces
    syn11 = get_trace(:syn11, logger.data)
    syn12 = get_trace(:syn12, logger.data)
    syn13 = get_trace(:syn13, logger.data)
    syn14 = get_trace(:syn14, logger.data)
    syn15 = get_trace(:syn15, logger.data)
    syn21 = get_trace(:syn21, logger.data)
    syn22 = get_trace(:syn22, logger.data)
    syn23 = get_trace(:syn23, logger.data)
    syn24 = get_trace(:syn24, logger.data)
    syn25 = get_trace(:syn25, logger.data)
    syn31 = get_trace(:syn31, logger.data)
    syn32 = get_trace(:syn32, logger.data)
    syn33 = get_trace(:syn33, logger.data)
    syn34 = get_trace(:syn34, logger.data)
    syn35 = get_trace(:syn35, logger.data)
    seg0  = get_trace(:seg0, logger.data)
    seg1  = get_trace(:seg1, logger.data)
    seg2  = get_trace(:seg2, logger.data)
    n     = get_trace(:n, logger.data)

    plateau0_starts = filter(x->(x.object == :seg0 && x.event == :plateau_starts), logger.data).t
    plateau1_starts = filter(x->(x.object == :seg1 && x.event == :plateau_starts), logger.data).t
    plateau2_starts = filter(x->(x.object == :seg2 && x.event == :plateau_starts), logger.data).t
    plateau1_ends = filter(x->(x.object == :seg1 && x.event == :plateau_ends), logger.data).t
    plateau2_ends = filter(x->(x.object == :seg2 && x.event == :plateau_ends), logger.data).t
    spike_times = filter(x->(x.object == :n && x.event == :spikes), logger.data).t
    # plot dynamics

    steps!(ax, [0;seg0.t;500], 0 .+ 4.9 .* [0;Int.(seg0.state .== ADSP.voltage_elevated);0], fill=color_3_25, color=:transparent)
    steps!(ax, [0;seg0.t;500], 0 .+ 4.9 .* [0;Int.(seg0.state .== ADSP.voltage_high);0], fill=color_3_50, color=:transparent)
    steps!(ax, [0;seg2.t;500], 5 .+  4.9 .* [0;Int.(seg2.state .== ADSP.voltage_elevated);0], fill=color_2_25, color=:transparent)
    steps!(ax, [0;seg2.t;500], 5 .+  4.9 .* [0;Int.(seg2.state .== ADSP.voltage_high);0], fill=color_2_50, color=:transparent)
    steps!(ax, [0;seg1.t;500], 10 .+ 4.9 .* [0;Int.(seg1.state .== ADSP.voltage_elevated);0], fill=color_1_25, color=:transparent)
    steps!(ax, [0;seg1.t;500], 10 .+ 4.9 .* [0;Int.(seg1.state .== ADSP.voltage_high);0], fill=color_1_50, color=:transparent)

    # plot spikes
    for (i, syn) in enumerate([syn11,syn12,syn13,syn14,syn15])
        steps!(ax, [0;syn.t;500], 9 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_1, color=:transparent)
    end
    
    for (i, syn) in enumerate([syn21,syn22,syn23,syn24,syn25])
        steps!(ax, [0;syn.t;500], 4 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_2, color=:transparent)
    end

    for (i, syn) in enumerate([syn31,syn32,syn33,syn34,syn35])
        steps!(ax, [0;syn.t;500], -1 .+ i .+ 0.9 .* [0;Int.(syn.state);0], fill=color_3, color=:transparent)
    end

    # plot output spikes
    lines!.(ax, Rect.(plateau0_starts.-5.5, -0.25, 11, 5.5), color=:gray10, linewidth=2)


    # plot neurons
    plot!(
        name == :and ? ax_right : ax_left, 
        objects[:n],
        angle_between = 20/180*pi,
        branch_length = 0.5,
        branch_width = 0.05,
        root_position = offset[name],
        color=Dict(:n=>color_3_25, :seg0=>color_3, :seg1=>color_1, :seg2=>color_2),
    )
end

fig

## Finish up and save
colsize!(fig.layout, 1, Relative(0.1))
colsize!(fig.layout, 3, Relative(0.1))
ylims!(ax_left, (-3.1,1.1))
ylims!(ax_right, (-3.1,1.1))
xlims!(ax_left, (-0.2,0.2))
xlims!(ax_right, (-0.2,0.2))
# xlims!(ax_input, (0,500))

save(joinpath("figures","motifs.pdf"), fig)
save(joinpath("figures","motifs.svg"), fig)
save(joinpath("figures","motifs.png"), fig)
fig
