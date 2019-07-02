# Active dendritic sequence processing (ADSP)

This package contains the accompanying code for the publication
"Event-based pattern detection in active dendrites"


# Files
This package is organized as follows:

* *model.jl*   defines the basic data structures used throughout the package
* *simulate.jl*   defines functions that implement basic simulation of neurons and networks
* *utils.jl*    defines additional utility functions for plotting, generating input etc.


# Examples

There are three scripts to be run to generate figures:

* `examples/grid_cells/neuron_paths.jl`
* `examples/grid_cells/ensemble_paths.jl` 
* `examples/zoo/zoo.jl` 

Executing `run_examples.jl` will execute all of them consecutively.
