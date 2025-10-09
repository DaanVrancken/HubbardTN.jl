# HubbardTN

Welcome to the documentation for HubbardTN, a tool for implementing and solving general 1D multi-band Hubbard models using tensor networks. The framework is built upon the packages [MPSKit.jl](https://github.com/QuantumKitHub/MPSKit.jl) and [TensorKit.jl](https://github.com/jutho/TensorKit.jl). [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/stable/) is used to automatically store your results in the desired location.

## Installation

To add this package to your Julia environment, do
```
julia> using Pkg
julia> Pkg.add(url="https://github.com/DaanVrancken/HubbardTN.jl")
```
after which it can be used by loading
```
julia> using HubbardTN
```

## Usage

The overall outline of the simulations is always the same. First, choose the type of Hubbard model that you are interested in. The different options are defined by their symmetries: 
- One-band model
    - Spin symmetry and conserved number of electrons, $\mathbb{Z}_2\times SU(2)\times U(1)$.
    - Conserved number of spin up and down electrons, $\mathbb{Z}_2\times U(1)\times U(1)$.
    - Spin symmetry, $\mathbb{Z}_2\times SU(2)$.
- Multi-band model
    - Spin symmetry and conserved number of electrons, $\mathbb{Z}_2\times SU(2)\times U(1)$.
    - Conserved number of spin up and down electrons, $\mathbb{Z}_2\times U(1)\times U(1)$.
    - Spin symmetry, $\mathbb{Z}_2\times SU(2)$.
Proceed by inserting the parameters of the model. You are now ready to compute the ground and excited states and their properties!
