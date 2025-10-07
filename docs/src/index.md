# HubbardTN

Welcome to the documentation for HubbardTN, a tool for implementing and solving general 1D multi-band Hubbard models using tensor networks.

## Installation

To reproduce this project, do the following:

1. Download this code base.
   ```
   git clone https://github.com/DaanVrancken/HubbardTN
   ```
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

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