# Examples

This page demonstrates how to use **HubbardTN** to build and solve different 1D Hubbard models using tensor networks.

Two examples are provided:
1. A **minimal single-band Hubbard model**
2. A **general multi-band Hubbard model**

---

## 🧩 Minimal Example

```julia
using HubbardTN
using TensorKit

# Step 1: Define the symmetries
particle_symmetry = U1Irrep
spin_symmetry = U1Irrep
cell_width = 2
filling = (1, 1)

symm = SymmetryConfig(particle_symmetry, spin_symmetry, cell_width, filling)

# Step 2: Set up model parameters
t = [0.0, 1.0]   # [chemical_potential, nn_hopping, nnn_hopping, ...]
U = [4.0]        # [on-site interaction, nn_interaction, ...]

model = ModelParams(t, U)
calc = CalcConfig(symm, model)

# Step 3: Compute the ground state
gs = compute_groundstate(calc)
ψ = gs["groundstate"]
H = gs["ham"]

println("Ground-state energy density: ", expectation_value(ψ, H) / length(H))

# Step 4: Compute first excitation in fZ2(0) × U1Irrep(0) × U1Irrep(0) sector
momenta = collect(range(0, 2π, length = 10))
charges = [0.0, 0.0, 0.0]
ex = compute_excitations(gs, momenta, charges)
```

### Notes

- The `SymmetryConfig` object defines all symmetry information of the model:
  - The **particle** and **spin** symmetries (`Trivial`, `U1Irrep`, or `SU2Irrep`)
  - The **number of sites in the unit cell** (`cell_width`)
  - The **filling fraction**, defined by `N_electrons / N_sites` via the keyword `filling=(N_electrons, N_sites)`
- The `ModelParams` constructor shown above is the simplest form, suitable for single-band Hubbard models.

---

## ⚙️ General Multi-Band Model

The example below illustrates a **two-band Hubbard model** with custom hopping and interaction terms.

```julia
using HubbardTN
using TensorKit

# Step 1: Define the symmetries
particle_symmetry = U1Irrep
spin_symmetry = SU2Irrep
cell_width = 2
filling = (1, 1)

symm = SymmetryConfig(particle_symmetry, spin_symmetry, cell_width, filling)

# Step 2: Define model parameters
bands = 2

# Hopping amplitudes:
# (1,2) and (2,1): inter-band hopping
# (2,3) and (3,2): next-nearest-neighbor hopping across unit cells
t = Dict((1,2)=>1.0, (2,1)=>1.0, (2,3)=>0.5, (3,2)=>0.5)

# Interaction terms:
# (i,j,k,l) correspond to U_ijkl c⁺_i c⁺_j c_k c_l
U = Dict(
    (1,1,1,1) => 8.0,   # on-site band 1
    (2,2,2,2) => 8.0,   # on-site band 2
    (1,2,1,2) => 1.0,   # inter-orbital exchange
    (2,1,2,1) => 1.0
)

model = ModelParams(bands, t, U)
calc = CalcConfig(symm, model)

# Step 3: Compute the ground state
gs = compute_groundstate(calc)
ψ = gs["groundstate"]
H = gs["ham"]

println("Ground-state energy density: ", expectation_value(ψ, H) / length(H))

# Step 4: Compute first excitations
momenta = collect(range(0, 2π, length = 10))
charges = [0.0, 0.0, 0.0]
ex = compute_excitations(gs, momenta, charges)
```

### Notes

- **Dictionaries** are used to define hopping (`t`) and interaction (`U`) parameters:
  - Keys represent index tuples:
    - `(i, j)` for hopping → corresponds to c⁺ᵢ cⱼ
    - `(i, j, k, l)` for interactions → corresponds to c⁺ᵢ c⁺ⱼ cₖ cₗ
  - Indices larger than the number of bands refer to orbitals in neighboring unit cells.
- This flexible formulation allows the user to specify *arbitrary band connectivity* and *interaction structure*.
- Beyond the usual on-site interactions, exchange or other couplings can be included naturally.

---

📘 **Tip:**  
For advanced use cases (custom operators, constrained symmetry sectors, or DMRG sweep control), see the API reference for  
[`compute_groundstate`](@ref), [`compute_excitations`](@ref), and [`ModelParams`](@ref).
