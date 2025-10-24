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