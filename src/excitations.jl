###############
# Excitations #
###############

"""
    compute_excitations(groundstate_dict, momenta, charges; nums=1, solver=Arnoldi(...))

Compute the low-lying quasiparticle excitations above a given ground state.

# Arguments
- `groundstate_dict::Dict{String,Any}`: A dictionary produced by `compute_groundstate`, 
   containing the keys `"groundstate"`, `"ham"`, and `"environments"`.
- `momenta::Union{Float64,Vector{Float64}}`: A collection of momentum values (in units of lattice sites) 
   at which excitations are evaluated. This argument is irrelevant for FiniteMPS.
- `charges::Vector{Float64}`: Target quantum numbers defining the excitation sector (one value per symmetry).
- `nums::Int64=1`: Number of excitations to compute per momentum.
- `solver`: The eigensolver used for diagonalization (default is `Arnoldi(; krylovdim=30, tol=1e-6, eager=true)`).

# Returns
A dictionary with the following keys:
- `"Es"`: Eigenenergies of the excitations.
- `"qps"`: Quasiparticle ansatz states.
- `"momenta"`: The input momentum values.
"""
function compute_excitations(
                    groundstate_dict::Dict{String,Any},
                    momenta::Union{Float64,Vector{Float64}},
                    charges::Union{Vector{Int64},Vector{Float64}};
                    nums::Int64=1,
                    solver=Arnoldi(;krylovdim=30,tol=1e-6,eager=true)
                )
    ψ = groundstate_dict["groundstate"]
    H = groundstate_dict["ham"]
    envs = groundstate_dict["environments"]

    trivial_sector = first(sectors(oneunit(physicalspace(H, 1))))
    @assert length(charges) == length(trivial_sector) "Number of charges must match number of symmetries ($(length(trivial_sector)))."
    sector = foldl(⊠, [typeof(f)(charges[i]) for (i, f) in enumerate((trivial_sector))])

    if typeof(ψ) <: FiniteMPS
        Es, qps = excitations(H, QuasiparticleAnsatz(solver, MPSKit.Defaults.alg_environments(;dynamic_tols=false)), 
                            ψ, envs; num=nums, sector=sector)
    else
        Es, qps = excitations(H, QuasiparticleAnsatz(solver, MPSKit.Defaults.alg_environments(;dynamic_tols=false)), 
                            momenta./length(H), ψ, envs; num=nums, sector=sector)
    end

    return Dict("Es" => Es, "qps" => qps, "momenta" => momenta)
end

"""
    compute_domainwall(groundstate_dict, momenta, charges; nums=1, shift=1, solver=Arnoldi(...))

Compute domain-wall excitations between a ground state and a spatially shifted version of itself.

# Arguments
- `groundstate_dict::Dict{String,Any}`: A dictionary produced by `compute_groundstate`,
   containing `"groundstate"`, `"ham"`, and `"environments"`.
- `momenta::Union{Float64,Vector{Float64}}`: A collection of momentum values (in units of lattice sites)
   at which domain-wall excitations are evaluated.
- `charges::Vector{Float64}`: Target quantum numbers defining the excitation sector (one value per symmetry).
- `nums::Int64=1`: Number of excitations to compute per momentum.
- `shift::Int64=1`: The number of lattice sites by which to shift the reference ground state to form the domain wall.
- `solver`: The eigensolver used for diagonalization (default is `Arnoldi(; krylovdim=30, tol=1e-6, eager=true)`).

# Returns
A dictionary with the following keys:
- `"Es"`: Eigenenergies of the domain-wall excitations.
- `"qps"`: Domain-wall quasiparticle ansatz states.
- `"momenta"`: The input momentum values.

# Notes
This function constructs the second “shifted” ground state and its environments automatically using `circshift` and `environments`.  
It then computes excitations between the two MPS states within the specified symmetry sector.
"""
function compute_domainwall(
                    groundstate_dict::Dict{String,Any},
                    momenta::Union{Float64,Vector{Float64}},
                    charges::Union{Vector{Int64},Vector{Float64}};
                    nums::Int64=1,
                    shift::Int64=1,
                    solver=Arnoldi(;krylovdim=30,tol=1e-6,eager=true)
                )
    ψ = groundstate_dict["groundstate"]
    H = groundstate_dict["ham"]
    envs = groundstate_dict["environments"]

    trivial_sector = first(sectors(oneunit(physicalspace(H, 1))))
    @assert length(charges) == length(trivial_sector) "Number of charges must match number of symmetries ($(length(trivial_sector)))."
    sector = foldl(⊠, [typeof(f)(charges[i]) for (i, f) in enumerate((trivial_sector))])

    ψ_s = circshift(ψ, shift)
    envs_s = environments(ψ_s, H);
    Es, qps = excitations(H, QuasiparticleAnsatz(solver, MPSKit.Defaults.alg_environments(;dynamic_tols=false)), 
                            momenta./length(H), ψ, envs, ψ_s, envs_s; num=nums, sector=sector)

    return Dict("Es" => Es, "qps" => qps, "momenta" => momenta)
end

"""
    compute_bandgap(gs, calc; resolution=5)

Compute the single-particle (charge) band gap from particle–addition and
particle–removal excitations over a uniform momentum grid.

# Arguments
- `gs::Dict{String,Any}`: Ground-state data as produced by `compute_groundstate`,
  containing the state, Hamiltonian, and environments.
- `calc::CalcConfig`: Calculation configuration. The particle symmetry must be
  `U1Irrep`.

# Keyword Arguments
- `resolution::Int=5`: Number of momentum points in the uniform grid between
  `0` and `π` (inclusive).

# Returns
- `(gap, kmin)`: The minimum value of `E_add(k) + E_remove(k)` over the sampled
  momenta and the corresponding momentum `k`.
"""
function compute_bandgap(gs::Dict{String,Any}, calc::CalcConfig; resolution::Int64=5)
    symm = calc.symmetries
    @assert symm.particle_symmetry==U1Irrep "Particle symmetry must be of type U1Irrep."
    d = denominator(symm.filling)
    if symm.spin_symmetry==Trivial
        charges_particle = [1, d]
    else
        charges_particle = [1, d, 1/2]
    end
    charges_hole = copy(charges_particle)
    charges_hole[2] *= -1

    momenta = collect(range(0, π, resolution))

    ex_add = compute_excitations(gs, momenta, charges_particle; nums=1)
    ex_rem = compute_excitations(gs, momenta, charges_hole; nums=1)

    Eb = ex_add["Es"] .+ ex_rem["Es"]
    gap, k = findmin(real.(Eb[:,1]))
    return gap, momenta[k]
end

"""
    compute_spingap(gs, calc; resolution=5)

Compute the spin gap from spin-flip excitations above the ground state over a
uniform momentum grid.

# Arguments
- `gs::Dict{String,Any}`: Ground-state data as produced by `compute_groundstate`,
  containing the state, Hamiltonian, and environments.
- `calc::CalcConfig`: Calculation configuration. The spin symmetry must not be
  `Trivial`.

# Keyword Arguments
- `resolution::Int=5`: Number of momentum points in the uniform grid between
  `0` and `π` (inclusive).

# Returns
- `(gap, kmin)`: The minimum spin excitation energy over the sampled momenta
  and the corresponding momentum `k`.
"""
function compute_spingap(gs::Dict{String,Any}, calc::CalcConfig; resolution::Int64=5)
    symm = calc.symmetries
    @assert symm.spin_symmetry!=Trivial "Spin symmetry must not be Trivial."
    if symm.particle_symmetry==Trivial
        charges = [0, 1]
    else
        charges = [0, 0, 1]
    end

    momenta = collect(range(0, π, resolution))
    ex = compute_excitations(gs, momenta, charges; nums=1)
    Es = ex["Es"]
    gap, k = findmin(real.(Es[:,1]))
    return gap, momenta[k]
end

"""
    compute_chargegap(gs, calc; resolution=5)

Compute the charge gap from two-particle addition excitations on a uniform momentum grid.

# Arguments
- `gs::Dict{String,Any}`: Ground-state data as produced by `compute_groundstate`,
  containing the state, Hamiltonian, and environments.
- `calc::CalcConfig`: Calculation configuration.

# Keyword Arguments
- `resolution::Int=5`: Number of momentum points in the uniform grid between
  `0` and `π` (inclusive).

# Returns
- `(gap, kmin)`: The minimum charge excitation energy over the sampled momenta
  and the corresponding momentum `k`.
"""
function compute_chargegap(gs::Dict{String,Any}, calc::CalcConfig; resolution::Int64=5)
    symm = calc.symmetries
    @assert symm.particle_symmetry==U1Irrep "Particle symmetry must be of type U1Irrep."
    d = denominator(symm.filling)
    if symm.spin_symmetry==Trivial
        charges_particle = [0, 2*d]
    else
        charges_particle = [0, 2*d, 0]
    end
    charges_hole = copy(charges_particle)
    charges_hole[2] *= -1

    momenta = collect(range(0, π, resolution))

    ex_add = compute_excitations(gs, momenta, charges_particle; nums=1)
    ex_rem = compute_excitations(gs, momenta, charges_hole; nums=1)

    Ec = 0.5 .* (ex_add["Es"] .+ ex_rem["Es"])
    gap, k = findmin(real.(Ec[:,1]))
    return gap, momenta[k]
end

"""
    compute_pairing_energy(gs, calc; resolution=5, tol::Float64=1e-6)

Compute the pairing energy from single-particle and two-particle addition excitations on a uniform momentum grid.

# Arguments
- `gs::Dict{String,Any}`: Ground-state data as produced by `compute_groundstate`, containing the state, Hamiltonian, and environments.
- `calc::CalcConfig`: Calculation configuration. The particle symmetry `calc.symmetries.particle_symmetry` must be `U1Irrep`.

# Keyword Arguments
- `resolution::Int=5`: Number of momentum points in the uniform grid between `0` and `π` (inclusive).

# Returns
- `(gap, kmin)`: The minimum value of `2E_add(k) - E_double(k)` over the sampled momenta and the corresponding momentum `k`.
"""
function compute_pairing_energy(gs::Dict{String,Any}, calc::CalcConfig; resolution::Int64=5, tol::Float64=1e-6)
    @assert calc.symmetries.particle_symmetry == U1Irrep "Particle symmetry must be of type U1Irrep."

    d = denominator(calc.symmetries.filling)

    if calc.symmetries.spin_symmetry == Trivial
        charges_particle = [1, d]
        charges_double = [0, 2*d]
    else
        charges_particle = [1, d, 1/2]
        charges_double = [0, 2*d, 0]
    end

    momenta = collect(range(0, π, resolution))

    ex_add = compute_excitations(gs, momenta, charges_particle; nums=1)
    ex_double = compute_excitations(gs, momenta, charges_double; nums=1)

    E1 = real.(vec(ex_add["Es"]))
    E2 = real.(vec(ex_double["Es"]))

    E1min, idx_E1min = findmin(E1)
    k_E1min = momenta[idx_E1min]

    E2min, idx_E2min = findmin(E2)
    k_E2min = momenta[idx_E2min]

    if !isapprox(k_E2min, 0.0; atol=tol)
        error("E2 minimum is not at k = 0.")
    end

    Ep = 2*E1min - E2min
    return Ep, k_E1min, k_E2min
end