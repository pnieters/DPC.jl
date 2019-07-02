import Base
using PlotTree


all_neurons = [
    # neuron 1
    Node("A"),

    # neuron 2
    Node([ Node("A"), Node("B") ], 1, "C" ),

    # neuron 3
    Node([ Node("A"), Node("B") ], 2, "C" ),

    # neuron 4
    Node([ Node("A"), Node("B"), Node("C") ], 1, "D" ),

    # neuron 5
    Node([ Node([ Node("A") ], 1, "B") ], 1, "C" ),
    
    # neuron 6
    Node([ Node([ Node("A"), Node("B") ], 2, "C"), Node("D") ], 1, "E" ),

    # neuron 7
    Node([ Node([ Node("A"), Node("B") ], 1, "E"), Node([Node("C"), Node("D") ], 2, "F") ], 2, "G" ),

    # neuron 8
    Node([ Node([ Node("A"), Node("B"), Node("C"), Node("D") ], 1, "I"), Node([ Node("E"), Node("F"), Node("G"), Node("H") ], 1, "J") ], 2, "K" ),

    # neuron 9
    Node([ Node([ Node([ Node("A"), ], 1, "C") ], 1, "E" ), Node([ Node([ Node("B"), ], 1, "D") ], 1, "F" )  ] ,1,"G"),
]

for (id, neuron) ∈ enumerate(all_neurons)
    save("examples/zoo/figures/neuron_$(id).svg", neuron; style=:tree, text_width=7)
    save("examples/zoo/figures/equation_$(id).svg", neuron; style=:equation, text_width=7)
end