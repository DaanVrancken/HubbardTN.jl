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