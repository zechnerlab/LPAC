# LPAC - Local Polynomial Approximation for Compartmentalised reaction systems

This repository contains the code used to analyse the case studies and produce the figures for the paper _Bianucci, T., Zechner, C. (2023 - in preparation). A local polynomial approximation method for compartmentalised biochemical systems._

## Usage
### Dependencies
This repository already contains `Project.toml` and `Manifest.toml` files specifying
is dependencies in the form of a `Pkg` environment.  
To use it:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()
```
or, in the Julia REPL:
```julia
] activate .
] instantiate
] precompile
```
This will download, install and precompile all the dependencies.

### General usage and project structure
The code is structured as a main module `LPAC` which directly exposes
the main data structures and plotting features, together with the main 
additional submodules `Sim`, `Models` and `Figures`.

### Generating the paper figures
The functions generating the figures for the paper are available in the 
`LPAC.Figures` submodule.  
They can be generated with the following:
```
julia --project
```
and
```julia
using LPAC
Figures.generateAllFigures()
```
This uses the same serialized simulation results that we used for generating the paper figures,
that are stored in `*.jser` file within the `./Figures` subfolder.
This ensures that the output is identical to the figures in the paper.

In order to run a new set of SSA trajectories and solve the moment equations again, 
it is enough to disable the loading of data from the serialized dump:
```julia
Figures.generateAllFigures(; loadFromDump=false)
```
The new simulation results are again also serialized and saved into `*.jser` files, so that 
tweaking of the plots does not require running the SSAs all over again.  
To load the results from these saves, use the `loadFromDump=true` flag (which is the default):
```julia
Figures.generateAllFigures(; loadFromDump=true)
```

Please note that any run with `loadFromDump=false` overwrites any existing `*.jser` files.
It is possible to revert to the original data from this repository by using the following command:
```
git checkout Figures
```
