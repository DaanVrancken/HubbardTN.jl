using HubbardTN
using TensorKit, MPSKit

# Step 0: Simple Hubbard model
particle_symmetry = U1Irrep
spin_symmetry = U1Irrep
cell_width = 2
filling = 1//1

symm = SymmetryConfig(particle_symmetry, spin_symmetry, cell_width, filling)

t = [0.0, 1.0]     # [chemical_potential, nn_hopping, nnn_hopping, ...]
U = [-10.0]        # [on-site interaction, nn_interaction, ...]
svalue = 2.5

model = HubbardParams(t, U)
calc0 = CalcConfig(symm, model)
gs0 = compute_groundstate(calc0; svalue=svalue)
ψ0 = gs0["groundstate"]
H0 = gs0["ham"]

E0 = sum(real(expectation_value(ψ0, H0))) / length(H0)
println("Groundstate energy: ", E0)

dim0 = dim_state(ψ0)
println("Max bond dimension: ", maximum(dim0))

Ne0 = density_e(ψ0, calc0)
println("Number of electrons per site: ", Ne0)

spin_gap, spin_k = compute_spingap(gs0, calc0)
pairing_gap, pairing_k = compute_pairing_energy(gs0, calc0)
println("Spin gap   at k = $(spin_k): Δs = $(spin_gap)")
println("Pairing gap   at k = $(pairing_k): Δp = $(pairing_gap)")

# Step 1: Define the symmetries
particle_symmetry = Trivial
symm = SymmetryConfig(particle_symmetry, spin_symmetry, cell_width)

# Step 2: Set up model parameters
t = [-5.0, 1.0]    # [chemical_potential, nn_hopping, nnn_hopping, ...]
U = [-10.0]        # [on-site interaction, nn_interaction, ...]
alpha = [0.0,0.0]
beta = [0.0,0.0]
model = HubbardParams(t, U)

alpha_list = Vector{Vector{Float64}}()
beta_list  = Vector{Vector{Float64}}()
E_list     = Float64[]

for i in 1:5
    global alpha, beta
    println("Step $i: alpha = $alpha, beta = $beta")
    calc = CalcConfig(symm, model, PairGapMF(alpha, beta))

    # Step 3: Compute the ground state
    gs = compute_groundstate(calc; svalue=svalue)
    ψ = gs["groundstate"]
    H = gs["ham"]

    E = sum(real(expectation_value(ψ, H))) / length(H)

    dim = dim_state(ψ)
    println("Max bond dimension: ", maximum(dim))

    Ne = density_e(ψ, calc)
    println("Number of electrons per site: ", Ne)

    push!(alpha_list, copy(alpha))
    push!(beta_list, copy(beta))
    push!(E_list, E)

    ty = 0.5
    tz = 0.5
    alpha = get_alpha(ψ, calc, ty, tz, pairing_gap)
    beta = get_beta(ψ, calc, ty, tz, pairing_gap)
end