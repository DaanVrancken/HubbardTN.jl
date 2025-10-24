###############
# Excitations #
###############

"""
    compute_excitations(groundstate_dict, momenta, charges; nums=1, solver=Arnoldi(...))

Compute the low-lying quasiparticle excitations above a given ground state.

# Arguments
- `groundstate_dict::Dict{String,Any}`: A dictionary produced by `compute_groundstate`, containing the keys `"groundstate"`, `"ham"`, and `"environments"`.
- `momenta::Union{Float64,Vector{Float64}}`: A collection of momentum values (in units of lattice sites) at which excitations are evaluated.
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
    @assert length(charges) == length(trivial_sector) "Number of charges must match number of symmetries."
    sector = foldl(⊠, [typeof(f)(charges[i]) for (i, f) in enumerate((trivial_sector))])

    Es, qps = excitations(H, QuasiparticleAnsatz(solver, MPSKit.Defaults.alg_environments(;dynamic_tols=false)), 
                            momenta./length(H), ψ, envs; num=nums, sector=sector)

    return Dict("Es" => Es, "qps" => qps, "momenta" => momenta)
end

"""
    compute_domainwall(groundstate_dict, momenta, charges; nums=1, shift=1, solver=Arnoldi(...))

Compute domain-wall excitations between a ground state and a spatially shifted version of itself.

# Arguments
- `groundstate_dict::Dict{String,Any}`: A dictionary produced by `compute_groundstate`, containing `"groundstate"`, `"ham"`, and `"environments"`.
- `momenta::Union{Float64,Vector{Float64}}`: A collection of momentum values (in units of lattice sites) at which domain-wall excitations are evaluated.
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
                    charges::Vector{Float64};
                    nums::Int64=1,
                    shift::Int64=1,
                    solver=Arnoldi(;krylovdim=30,tol=1e-6,eager=true)
                )
    ψ = groundstate_dict["groundstate"]
    H = groundstate_dict["ham"]
    envs = groundstate_dict["environments"]

    trivial_sector = first(sectors(oneunit(physicalspace(H, 1))))
    @assert length(charges) == length(trivial_sector) "Number of charges must match number of symmetries."
    sector = foldl(⊠, [typeof(f)(charges[i]) for (i, f) in enumerate((trivial_sector))])

    ψ_s = circshift(ψ, shift)
    envs_s = environments(ψ_s, H);
    Es, qps = excitations(H, QuasiparticleAnsatz(solver, MPSKit.Defaults.alg_environments(;dynamic_tols=false)), 
                            momenta./length(H), ψ, envs, ψ_s, envs_s; num=nums, sector=sector)

    return Dict("Es" => Es, "qps" => qps, "momenta" => momenta)
end
