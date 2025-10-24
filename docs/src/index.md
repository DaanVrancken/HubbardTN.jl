# HubbardTN

Welcome to the documentation for **HubbardTN**, a Julia package for constructing and solving general **1D multi-band Hubbard models** using **tensor network** techniques.  
The framework builds on top of [MPSKit.jl](https://github.com/QuantumKitHub/MPSKit.jl) and [TensorKit.jl](https://github.com/jutho/TensorKit.jl).

---

## Installation

You can install `HubbardTN` directly from its GitHub repository:

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/DaanVrancken/HubbardTN.jl")
```
After installation, load the package with:
```
julia> using HubbardTN
```

## Usage Overview

A typical HubbardTN simulation follows a clear sequence of steps:

1. **Define the symmetries.**  
   Every model includes a built-in fermionic **ℤ₂ symmetry**.  
   For the particle and spin symmetries, you can choose among:
   - `Trivial`
   - `U1Irrep`
   - `SU2Irrep`  

   These are specified in a [`SymmetryConfig`](@ref) object.

2. **Specify the model parameters.**  
   Insert your model’s coupling constants and hopping terms into a [`ModelParams`](@ref) object.  
   This defines the Hamiltonian that will be constructed.

3. **Compute physical quantities.**  
   With the symmetry and model defined, you can:
   - Compute the **ground state** using [`compute_groundstate`](@ref).
   - Obtain **excitations** with [`compute_excitations`](@ref).
   - Investigate **domain walls** via [`compute_domainwall`](@ref).
