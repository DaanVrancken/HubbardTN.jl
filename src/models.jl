##########################
# Symmetry configuration #
##########################

"""
    SymmetryConfig(particle_symmetry, spin_symmetry, cell_width, filling=nothing)

Represents the symmetry configuration of a lattice system, including particle and spin symmetries,
unit cell width, and optional filling information.

# Fields
- `particle_symmetry` : Union{Type{Trivial}, Type{U1Irrep}, Type{SU2Irrep}}
    - The symmetry type for particle number. Use `Trivial` for no symmetry, `U1Irrep` for U(1) symmetry, or `SU2Irrep` for SU(2) symmetry.
- `spin_symmetry` : Union{Type{Trivial}, Type{U1Irrep}, Type{SU2Irrep}}
    - The symmetry type for spin degrees of freedom.
- `cell_width` : Int64
    - Number of sites in the unit cell. Must be a positive integer.
- `filling` : Union{Nothing, Rational{Int}}
    - Optional particle filling specified as a rational number `P//Q` (numerator/denominator). Only allowed if `particle_symmetry` is `U1Irrep`.
    Otherwise the filling is determined by the chemical potential.

# Constructor Behavior
- If `filling` is provided, the constructor checks that `particle_symmetry == U1Irrep`.
- `cell_width` must be positive.
- If `particle_symmetry == U1Irrep` and `filling === nothing`, filling defaults to `1//1`.
"""
struct SymmetryConfig
    particle_symmetry::Union{Type{Trivial},Type{U1Irrep},Type{SU2Irrep}}
    spin_symmetry::Union{Type{Trivial},Type{U1Irrep},Type{SU2Irrep}}
    cell_width::Int64
    filling::Union{Nothing,Rational{Int}}

    function SymmetryConfig(
                particle_symmetry::Union{Type{Trivial},Type{U1Irrep},Type{SU2Irrep}},
                spin_symmetry::Union{Type{Trivial},Type{U1Irrep},Type{SU2Irrep}},
                cell_width::Int64,
                filling::Union{Nothing,Rational{Int}}=nothing
            )
        cell_width > 0 || throw(ArgumentError("Cell width must be a positive integer, got $cell_width."))

        if particle_symmetry == U1Irrep
            filling = filling === nothing ? 1//1 : filling
        elseif filling !== nothing
            throw(ArgumentError("Filling can only be specified when particle symmetry is U1Irrep, but got $(particle_symmetry)."))
        end

        return new(particle_symmetry, spin_symmetry, cell_width, filling)
    end
end


#################
# Hubbard model #
#################

# Convert hopping matrix to dictionary representation
function hopping_matrix2dict(t::Vector{T}) where {T<:AbstractFloat}
    hopping = Dict{NTuple{2,Int},T}()
    for (j, val) in enumerate(t)
        if val != 0.0
            hopping[(1,j)] = val
            hopping[(j,1)] = conj(val)
        end
    end

    return hopping
end
function hopping_matrix2dict(t::Matrix{T}) where {T<:AbstractFloat}
    hopping = Dict{NTuple{2,Int},T}()
    bands = size(t,1)
    size(t,2) % bands == 0 || throw(ArgumentError("Second dimension of t ($(size(t,2))) must be a multiple of number of bands ($bands)."))
    ishermitian(t[1:bands, 1:bands]) || throw(ArgumentError("t on-site matrix is not Hermitian."))

    for i in 1:bands
        for j in 1:size(t,2)
            if t[i, j] != 0.0
                hopping[(i, j)] = t[i, j]
                hopping[(j, i)] = conj(t[i, j])
            end
        end
    end

    return hopping
end
# Add 2-body interaction term and its Hermitian conjugate
function addU!(U::Dict{NTuple{4,Int},T}, key::NTuple{4,Int}, val::T) where {T}
    if val != 0
        U[key] = val
        i,j,k,l = key
        U[(l,k,j,i)] = conj(val)          # Hermitian conjugate term
    end
end
# Check Hermiticity of parameter dictionary
function check_hermitian_dict(d::Dict{NTuple{X, Int}, T}; atol::Real=1e-8, rtol::Real=1e-5) where {X, T}
    for (k, val) in d
        k_rev = reverse(k)
        
        # Check if the conjugate/inverted key exists
        if !haskey(d, k_rev)
            return false, k_rev
        end
        
        # Check if values match within tolerance
        val_rev = d[k_rev]
        if !isapprox(val, val_rev; atol=atol, rtol=rtol)
            return false, k_rev
        end
    end
    return true, nothing
end

"""
    HubbardParams{T<:AbstractFloat}

Represents the standard Hubbard Hamiltonian parameters for a lattice or
multi-orbital system.

# Fields
- `bands::Int64`  
    Number of electronic orbitals or bands per unit cell. Must be positive.
- `t::Dict{NTuple{2, Int64}, T}`  
    Hopping amplitudes. Convention: `t[(i,i)] = μ_i` is the on-site potential,
    `t[(i,j)]` for `i ≠ j` is the hopping amplitude from site i to j.
- `U::Dict{NTuple{4, Int64}, T}`  
    Two-body electronic interaction tensor. Entries `U[(i,j,k,l)]` correspond
    to the operator c⁺_i c⁺_j c_k c_l. Zero entries can be omitted.

# Constructors
- `HubbardParams(bands, t::Dict, U::Dict)` — standard constructor specifying bands, hopping, and interactions.
- `HubbardParams(t::Vector, U::Vector)` — single-band convenience constructor from vectors.
- `HubbardParams(t::Matrix, U::Matrix)` — multi-band constructor from matrices; automatically checks dimensions and Hermiticity.
"""
struct HubbardParams{T<:AbstractFloat}
    bands::Int64
    t::Dict{NTuple{2, Int64}, T}          # t_ii=µ_i, t_ij hopping i→j
    U::Dict{NTuple{4, Int64}, T}          # U_ijkl c⁺_i c⁺_j c_k c_l

    function HubbardParams(bands::Int64, t::Dict{NTuple{2,Int64}, T}, U::Dict{NTuple{4,Int},T}) where {T<:AbstractFloat}
        bands > 0 || throw(ArgumentError("Number of bands must be a positive integer, got $bands."))
        all(k -> all(>(0), k), keys(t)) || throw(ArgumentError("t has negative indices."))
        all(k -> all(>(0), k), keys(U)) || throw(ArgumentError("U has negative indices."))
        t_hermitian, key_t = check_hermitian_dict(t)
        t_hermitian || throw(ArgumentError("t is not Hermitian. Missing or inconsistent conjugate for key $(key_t)."))
        U_hermitian, key_U = check_hermitian_dict(U) 
        U_hermitian || throw(ArgumentError("U is not Hermitian. Missing or inconsistent conjugate for key $(key_U)."))
        return new{T}(bands, t, U)
    end
end
# Constructors
function HubbardParams(t::Union{Vector{T}, Matrix{T}}, U::Dict{NTuple{4,Int},T}) where {T<:AbstractFloat}
    bands = isa(t, Matrix) ? size(t,1) : 1
    return HubbardParams(bands, hopping_matrix2dict(t), U)
end
function HubbardParams(t::Vector{T}, U::Vector{T}) where {T<:AbstractFloat}
    interaction = Dict{NTuple{4,Int},T}()
    for (i, val) in enumerate(U)
        addU!(interaction, (1,i,i,1), val)
        addU!(interaction, (i,1,1,i), val)
    end
    return HubbardParams(1, hopping_matrix2dict(t), interaction)
end
function HubbardParams(t::Matrix{T}, U::Matrix{T}) where {T<:AbstractFloat}
    bands = size(t, 1)

    # --- basic checks ---
    size(U, 1) == bands || throw(ArgumentError("First dimension of U ($(size(U,1))) must be equal to number of bands ($bands)."))
    size(U, 2) % bands == 0 || throw(ArgumentError("Second dimension of U ($(size(U,2))) must be multiple of number of bands ($bands)."))
    ishermitian(U[1:bands, 1:bands]) || throw(ArgumentError("U on-site matrix is not Hermitian."))

    interaction = Dict{NTuple{4,Int},T}()
    for i in 1:bands
        for j in 1:size(U,2)
            addU!(interaction, (i,j,j,i), U[i,j])
            j>bands && addU!(interaction, (j,i,i,j), U[i,j])
        end
    end

    return HubbardParams(bands, hopping_matrix2dict(t), interaction)
end


###############
# Extra terms #
###############

abstract type AbstractHamiltonianTerm end

# Add 3-body interaction term and its Hermitian conjugate
function addV!(V::Dict{NTuple{6,Int},T}, key::NTuple{6,Int}, val::T) where {T<:AbstractFloat}
    if val != 0
        V[key] = val
        i,j,k,l,n,m = key
        V[(m,n,l,k,j,i)] = conj(val)          # Hermitian conjugate term
    end
end
# Add standard 3-body interaction terms to dictionary
function dominant_threebody_term(i::Int64, j::Int64, value::T, V::Dict{NTuple{6,Int},T}=Dict{NTuple{6,Int},T}()) where {T<:AbstractFloat}
    addV!(V, (i,i,j,j,i,i), value)
    addV!(V, (i,j,i,i,j,i), value)
    addV!(V, (j,i,i,i,i,j), value)
    addV!(V, (i,j,j,j,j,i), value)
    addV!(V, (j,i,j,j,i,j), value)
    addV!(V, (j,j,i,i,j,j), value)
    return V
end

"""
    ThreeBodyTerm{T<:AbstractFloat} <: AbstractHamiltonianTerm

Represents three-body interactions in the Hamiltonian.

# Fields
- `bands::Int64`  
    Number of bands (orbitals) in the system.
- `V::Dict{NTuple{6,Int}, T}`  
    Three-body interaction amplitudes `c⁺_i c⁺_j c⁺_k c_l c_m c_n`. Zero
    entries can be omitted.

# Constructors
- `ThreeBodyTerm(V::Vector{T})` — single-band constructor from a vector.
- `ThreeBodyTerm(V::Matrix{T})` — multi-band constructor from a matrix.
"""
struct ThreeBodyTerm{T<:AbstractFloat} <: AbstractHamiltonianTerm
    bands::Int64
    V::Dict{NTuple{6,Int}, T}
    function ThreeBodyTerm(bands::Int64, V::Dict{NTuple{6,Int}, T}) where {T<:AbstractFloat}
        all(k -> all(>(0), k), keys(V)) || throw(ArgumentError("V has negative indices."))
        V_hermitian, key_V = check_hermitian_dict(V)
        V_hermitian || throw(ArgumentError("V is not Hermitian. Missing or inconsistent conjugate for key $(key_V)."))
        return new{T}(bands, V)
    end
end
# Constructors
function ThreeBodyTerm(V::Vector{T}) where {T<:AbstractFloat}
    threebody = Dict{NTuple{6,Int},T}()
    for (i, val) in enumerate(V)
        threebody = dominant_threebody_term(1, i+1, val, threebody)
    end

    return ThreeBodyTerm(1, threebody)
end
function ThreeBodyTerm(V::Matrix{T}) where {T<:AbstractFloat}
    bands = size(V, 1)
    size(V, 2) % bands == 0 || throw(ArgumentError("Second dimension of V ($(size(V,2))) must be multiple of number of bands ($bands)."))
    ishermitian(V[1:bands, 1:bands]) || throw(ArgumentError("V on-site matrix is not Hermitian."))

    threebody = Dict{NTuple{6,Int},T}()
    for i in 1:bands
        for j in 1:size(V,2)
            threebody = dominant_threebody_term(i, j, V[i,j], threebody)
        end
    end

    return ThreeBodyTerm(bands, threebody)
end

"""
    MagneticField{T<:AbstractFloat} <: AbstractHamiltonianTerm

Represents a magnetic field term in the Hamiltonian `B * Sᶻ`.

# Fields
- `B::T`  
    Magnetic field strength.

# Constructors
- `MagneticField(B)` — creates the term with specified magnetic field strength.
"""
struct MagneticField{T<:AbstractFloat} <: AbstractHamiltonianTerm
    B::T            # Magnetic field strength
end

"""
    StaggeredField{T<:AbstractFloat} <: AbstractHamiltonianTerm

Represents a staggered magnetic field term in the Hamiltonian. 
For multi-band models, the staggering is applied between equivalent orbitals.

# Fields
- `J::T`  
    Inter-chain Hund's coupling.
- `Ms::T`  
    Initial staggered magnetization.

# Constructors
- `StaggeredField(J, Ms)` — creates the term with specified coupling and initial magnetization.
"""
struct StaggeredField{T<:AbstractFloat} <: AbstractHamiltonianTerm
    J::T            # Inter-chain Hund's coupling
    Ms::T           # Initial staggered magnetization
end

"""
    SpinMeanField{T<:AbstractFloat} <: AbstractHamiltonianTerm

Mean-field coupling term representing inter-chain spin interactions. 

# Fields
- `J::Matrix{T}`: Inter-chain coupling matrix of size NxN, where J[i,j] couples site i
  in the current chain to site j in neighboring chains. Multiply with the coordination number
  to account for multiple chains.
- `spins::Union{Vector{T}, Matrix{T}}`: Expected spin expectation values from neighboring chains.
    - **Collinear**: A `Vector{T}` of length N containing z-components.
    - **Noncollinear**: A `Matrix{T}` of size Nx3 containing (x, y, z) components.

# Constructors
- `SpinMeanField(J, spins)`: Creates the term with specified coupling and initial spins.
"""
struct SpinMeanField{T<:AbstractFloat} <: AbstractHamiltonianTerm
    J::Matrix{T}
    spins::Union{Vector{T}, Matrix{T}}

    function SpinMeanField(J::Matrix{T}, spins::Vector{T}) where {T<:AbstractFloat}
        @assert size(J, 1) == size(J, 2) "Coupling matrix J must be square."
        @assert size(J, 1) == length(spins) "Number of sites in J must match length of spins vector."
        return new{T}(J, spins)
    end
    function SpinMeanField(J::Matrix{T}, spins::Matrix{T}) where {T<:AbstractFloat}
        @assert size(J, 1) == size(J, 2) "Coupling matrix J must be square."
        @assert size(spins, 2) == 3 "Noncollinear spins matrix must have exactly 3 columns (x, y, z)."
        @assert size(J, 1) == size(spins, 1) "Number of sites in J must match number of rows in spins matrix."
        return new{T}(J, spins)
    end
end

abstract type AbstractInterchainMF <: AbstractHamiltonianTerm end

"""
    ChargeGapMF{T<:Real,S<:Number} <: AbstractInterchainMF

Terms used in a perturbative treatment based on the charge gap, parametrizing
effective interchain/interladder processes.

# Fields
- `t_inter::Dict{NTuple{2, Int64}, T}`
    Inter-chain hopping parameters. `t_inter[(i,j)]` is the hopping amplitude 
    from site i on chain 0 to site j on the neigboring chain. Has to be scaled with `√(z/Δ)`,
    where `z` is the coordination number and `Δ` the charge/band gap.
- `bands::Int64`
    Number of bands per unit cell. Must match the number defined in HubbardParams.
- `cell_width::Int64`
    Number of sites in the unit cell. Must match the number defined in SymmetryConfig.
- `range::Int64`
    Maximum unit-cell separation, beyond whatever spread `t_inter` itself already spans,
    between the two neighboring-chain sites tracked in the `beta_*` correlators.
- `beta_uu::Dict{NTuple{2, Int64}, S}`
    Dictionary of self-consistent parameters `⟨cₖ↑⁺cₗ↑⟩`.
- `beta_ud::Dict{NTuple{2, Int64}, S}`
    Dictionary of self-consistent parameters `⟨cₖ↑⁺cₗ↓⟩`. The corresponding
    `⟨cₖ↓⁺cₗ↑⟩` parameters are obtained from its adjoint.
- `beta_dd::Dict{NTuple{2, Int64}, S}`
    Dictionary of self-consistent parameters `⟨cₖ↓⁺cₗ↓⟩`.

    Each beta dictionary must contain exactly the keys generated from pairs of
    keys in `t_inter`. For every `(a, i)` and `(b, j)` in `keys(t_inter)`, and
    for every cell offset `cw = 0:cell_width-1` and relative displacement
    `r = -range*bands:bands:range*bands`, the required key is
    `(i + cw*bands, j + r + cw*bands)` shifted so that the first index is in
    the range `[1,bands*cell_width]`. Thus, all required combinations must be
    present and no other keys are allowed.

# Constructors
- `ChargeGapMF(t_inter, bands, cell_width, range, beta_uu, beta_ud, beta_dd)`
    Constructor accepting explicitly initialized beta dictionaries.
- `ChargeGapMF(t_inter, bands, cell_width, range)`
    Convenience constructor that creates all three beta dictionaries with the
    required keys and initializes every value to zero.

# Notes
- The beta dictionaries are not fixed couplings: they should be iterated to
    convergence together with the ground state (or other target state) to
    satisfy the chosen self-consistency condition.
"""
struct ChargeGapMF{T<:Real,S<:Number} <: AbstractInterchainMF 
    t_inter::Dict{NTuple{2, Int64}, T}
    bands::Int64
    cell_width::Int64
    range::Int64
    beta_uu::Dict{NTuple{2, Int64}, S}
    beta_ud::Dict{NTuple{2, Int64}, S}
    beta_dd::Dict{NTuple{2, Int64}, S}
    function ChargeGapMF(t_inter::Dict{NTuple{2, Int64}, T}, bands::Int64, cell_width::Int64, range::Int64,
                beta_uu::Dict{NTuple{2, Int64}, S}, beta_ud::Dict{NTuple{2, Int64}, S}, beta_dd::Dict{NTuple{2, Int64}, S}
            ) where {T<:Real,S<:Number}
        bands > 0 || throw(ArgumentError("bands must be a positive integer, got $bands."))
        cell_width > 0 || throw(ArgumentError("cell_width must be a positive integer, got $cell_width."))
        range >= 0 || throw(ArgumentError("range must be a positive integer, got $range."))
        all(k -> all(>(0), k[1]), keys(t_inter)) || throw(ArgumentError("t_inter has negative first index."))
        
        expected_keys = Set{NTuple{2, Int64}}()
        period = cell_width * bands
        for (_, i) in keys(t_inter), (_, j) in keys(t_inter), r in -range*bands:bands:range*bands, cell in 0:cell_width-1
            idx = (i + cell*bands, j + r + cell*bands)
            # Shift to have first index in central unit cell
            shift = mod1(idx[1], period) - idx[1]
            idx = (idx[1] + shift, idx[2] + shift)
            push!(expected_keys, idx)
        end

        for (dict, name) in ((beta_uu, "beta_uu"), (beta_ud, "beta_ud"), (beta_dd, "beta_dd"))
            dict_keys = keys(dict)
            missing_keys = setdiff(expected_keys, dict_keys)
            isempty(missing_keys) || throw(ArgumentError("$name is missing elements: $(join(missing_keys, ", "))"))
            extra_keys = setdiff(dict_keys, expected_keys)
            isempty(extra_keys) || throw(ArgumentError("$name contains extra elements: $(join(extra_keys, ", "))"))
        end

        return new{T,S}(t_inter, bands, cell_width, range, beta_uu, beta_ud, beta_dd)
    end
end
# Constructor
function ChargeGapMF(
            t_inter::Dict{NTuple{2, Int64}, T}, bands::Int64, cell_width::Int64, range::Int64
        ) where {T<:Real}
    beta = Dict{NTuple{2, Int64}, ComplexF64}()
    period = cell_width * bands
    for (_, i) in keys(t_inter), (_, j) in keys(t_inter), r in -range*bands:bands:range*bands, cell in 0:cell_width-1
        idx = (i + cell*bands, j + r + cell*bands)
        shift = mod1(idx[1], period) - idx[1]
        idx = (idx[1] + shift, idx[2] + shift)
        beta[idx] = 0.0 + 0.0im
    end

    return ChargeGapMF(t_inter, bands, cell_width, range, beta, copy(beta), copy(beta))
end

"""
    PairGapMF{T<:AbstractFloat} <: AbstractHamiltonianTerm

Terms used in a perturbative treatment based on the pair gap, parameterizing 
effective interchain/interladder processes. 
Ref: Bollmark et al., Phys. Rev. X 13, 011039 (2023)

# Fields
- `alpha::Vector{T}`
    Self-consistent parameters associated with pair-tunneling.
- `beta::Vector{T}`
    Self-consistent parameters associated with exchange.
# Notes
- `alpha` and `beta` are not fixed couplings: they should be iterated to convergence together
  with the ground state (or other target state) to satisfy the chosen self-consistency condition.
"""
struct PairGapMF{T<:AbstractFloat} <: AbstractInterchainMF 
    alpha::Vector{T}
    beta::Vector{T}
end

"""
    HolsteinTerm{T<:AbstractFloat} <: AbstractHamiltonianTerm

Represents Holstein-type electron–phonon coupling terms `w b⁺ᵢ bᵢ` and
`gₐ(nᵢₐ-<n>)(b⁺ⱼ + bⱼ)` in the Hamiltonian. The coupling may be local
(`i=j`) or decaying with distance `rᵢⱼ^(-ξ)`. Couplings
smaller than `threshold` are neglected.

# Fields
- `w::Vector{T}`  
    Local phonon frequencies (one per mode).
- `g::Matrix{T}`  
    Electron–phonon coupling strengths, size `(bands, nmodes)`.
- `max_b::Int64`  
    Maximum number of phonons allowed per mode (Fock-space truncation).
- `mean_ne::T`  
    Mean number of electrons per site in the bare Hubbard model, used to
    normal-order the density operator.
- `xi::T=Inf`  
    Power law decay length for non-local coupling (`Inf` means strictly local).
- `threshold::T=0`  
    Minimum coupling magnitude retained in the Hamiltonian; smaller values are
    dropped for efficiency.

# Constructor
    HolsteinTerm(w, g, max_b, mean_ne; xi=zero(T), threshold=zero(T))

All arguments are positional except `xi` and `threshold`, which are keyword
arguments with default `0`.  Pass `xi < Inf` together with a positive
`threshold` to enable decaying non-local coupling.
"""
struct HolsteinTerm{T<:AbstractFloat} <: AbstractHamiltonianTerm
    w::Vector{T}
    g::Matrix{T}
    max_b::Int64
    mean_ne::T
    xi::T
    threshold::T

    function HolsteinTerm(
                w::Vector{T},
                g::Matrix{T},
                max_b::Int64,
                mean_ne::T;
                xi::T=Inf,
                threshold::T=zero(T)
            ) where {T<:AbstractFloat}
        max_b > 0 || throw(ArgumentError("max_b must be a positive integer, got $max_b."))
        size(g, 2) == length(w) || throw(ArgumentError(
            "Number of columns of g ($(size(g,2))) must equal length of w ($(length(w)))."))
        xi > zero(T) || throw(ArgumentError("xi must be positive, got $xi."))
        threshold >= zero(T) || throw(ArgumentError("threshold must be non-negative, got $threshold."))
        !isfinite(xi) || threshold > zero(T) || throw(ArgumentError(
            "A positive threshold is required when xi < Inf to avoid retaining negligibly small long-range couplings."))
        return new{T}(w, g, max_b, mean_ne, xi, threshold)
    end
end

"""
    ImpurityTerm{T<:AbstractFloat} <: AbstractHamiltonianTerm

Represents a local impurity by modifying hopping and two-body interaction
parameters at the specified impurity sites.

# Fields
- `t_imp::Dict{NTuple{2, Int64}, T}`  
    Changes in the hopping amplitudes relative to the bare Hubbard model.
    Entries `t_imp[(i,j)]` correspond to `Δt_ij = t'_ij - t_ij`, where the
    indices `i` and `j` directly specify the sites associated with the impurity.
- `U_imp::Dict{NTuple{4, Int64}, T}`  
    Changes in the two-body interaction tensor relative to the bare Hubbard model.
    Entries `U_imp[(i,j,k,l)]` correspond to `ΔU_ijkl = U'_ijkl - U_ijkl`, where
    the indices `i`, `j`, `k`, and `l` directly specify the sites associated with
    the impurity.

# Constructors
- `ImpurityTerm(t_imp, U_imp)` — creates an impurity term with specified changes
  to the hopping and interaction parameters at the impurity sites.
"""
struct ImpurityTerm{T<:AbstractFloat} <: AbstractHamiltonianTerm
    t_imp::Dict{NTuple{2, Int64}, T}
    U_imp::Dict{NTuple{4, Int64}, T}

    function ImpurityTerm(
                t_imp::Dict{NTuple{2,Int64}, T},
                U_imp::Dict{NTuple{4,Int},T}
            ) where {T<:AbstractFloat}
        return new{T}(t_imp, U_imp)
    end
end

######################
# Calculation set up #
######################

# Check for duplicate Hamiltonian terms
function get_family(::Type{T}) where {T<:AbstractHamiltonianTerm}
    S = supertype(T)
    (S === AbstractHamiltonianTerm || S === Any) ? T : get_family(S)
end
function find_duplicate(terms::Tuple)
    for (i, t) in enumerate(terms)
        T = get_family(typeof(t))
        for s in terms[i+1:end]
            T2 = get_family(typeof(s))
            T2 === T && return T
        end
    end
    return nothing
end

"""
    CalcConfig{T<:AbstractFloat, HamiltonianTerms<:Tuple{Vararg{AbstractHamiltonianTerm}}}

Holds all configuration information for a lattice or many-body calculation,
including symmetries, the base Hubbard Hamiltonian, and optional additional Hamiltonian terms.

# Fields
- `symmetries::SymmetryConfig`  
    Contains particle-number and spin symmetries, unit-cell geometry, and optional filling.
- `hubbard::HubbardParams{T}`  
    Base electronic Hamiltonian parameters (bands, hopping, and two-body interactions).
- `terms::HamiltonianTerms`  
    Tuple of additional Hamiltonian terms (subtypes of `AbstractHamiltonianTerm`).

# Constructors
- `CalcConfig(symmetries, hubbard, terms)` — creates a configuration with specified
  symmetries, Hubbard parameters, and extra terms. Duplicate term types are checked,
  and band/shape consistency is enforced.
- `CalcConfig(symmetries, hubbard, term)` — convenience constructor with one extra term (`terms = (term,)`).
- `CalcConfig(symmetries, hubbard)` — convenience constructor with no extra terms (`terms = ()`).

# Validation notes
- `CalcConfig` throws `ArgumentError` for invalid user inputs, including duplicate term types,
  incompatible term dimensions, and invalid filling/cell-width compatibility checks.
"""
struct CalcConfig{
    T<:AbstractFloat,
    HamiltonianTerms<:Tuple{Vararg{AbstractHamiltonianTerm}}
}
    symmetries::SymmetryConfig
    hubbard::HubbardParams
    terms::HamiltonianTerms

    function CalcConfig(
                symmetries::SymmetryConfig,
                hubbard::HubbardParams{T},
                terms::HamiltonianTerms
            ) where {T<:AbstractFloat, HamiltonianTerms<:Tuple{Vararg{AbstractHamiltonianTerm}}}

        dup = find_duplicate(terms)
        dup === nothing || throw(ArgumentError("Duplicate Hamiltonian term detected: $dup."))

        bands = hubbard.bands
        cw = symmetries.cell_width
        expected_sites = bands * cw

        if symmetries.filling !== nothing
            n = numerator( symmetries.filling)
            d = denominator(symmetries.filling)
            (n > 0 && d > 0) || throw(ArgumentError("Filling numerator and denominator must be positive integers, got $n//$d."))
            necessary_width = d * (mod(n, 2) + 1)
            cw % necessary_width == 0 || throw(ArgumentError("cell_width ($(cw)) must be a multiple of $necessary_width to accommodate the specified filling ($n / $d)."))
        end

        for term in terms
            # each dict indices are only interpretable given bands; this is cross-checked against HubbardParams.bands"
            if :bands in fieldnames(typeof(term))
                term.bands == bands || throw(ArgumentError("Number of bands in HubbardParams ($bands) does not match number of bands in $(typeof(term)) ($(term.bands))."))
            end
            if :cell_width in fieldnames(typeof(term))
                term.cell_width == cw|| throw(ArgumentError("Number of bands in HubbardParams ($bands) does not match number of bands in $(typeof(term)) ($(term.bands))."))
            end
            if term isa HolsteinTerm
                size(term.g, 1) == bands || throw(ArgumentError("Number of bands in HubbardParams ($bands) does not match first dimension of HolsteinTerm.g ($(size(term.g,1)))."))
            elseif term isa SpinMeanField
                size(term.J, 1) == expected_sites || throw(ArgumentError("Number of electron sites in cell ($expected_sites) does not match first dimension of SpinMeanField.J ($(size(term.J,1)))."))
            elseif term isa ChargeGapMF
                max_index = maximum(k[1] for k in keys(term.t_inter))
                max_index <= bands || throw(ArgumentError("Index in ChargeGapMF.t_inter ($(max_index)) exceeds number of bands ($bands)."))
            elseif term isa ImpurityTerm
                all(k -> all(i -> 1 <= i <= expected_sites, k), keys(term.t_imp)) || throw(ArgumentError("All indices in ImpurityTerm.t_imp must be between 1 and $expected_sites."))
                all(k -> all(i -> 1 <= i <= expected_sites, k), keys(term.U_imp)) || throw(ArgumentError("All indices in ImpurityTerm.U_imp must be between 1 and $expected_sites."))
            end
            if symmetries.filling !== nothing && term isa HolsteinTerm
                newf = symmetries.filling * bands // (bands + length(term.w))
                
                symmetries = SymmetryConfig(
                    symmetries.particle_symmetry,
                    symmetries.spin_symmetry,
                    symmetries.cell_width,
                    newf
                )
            end
        end

        return new{T, HamiltonianTerms}(symmetries, hubbard, terms)
    end
    CalcConfig(
            symmetries::SymmetryConfig, 
            hubbard::HubbardParams{T}, 
            term::AbstractHamiltonianTerm
        ) where {T<:AbstractFloat} = CalcConfig(symmetries, hubbard, (term,))
    CalcConfig(symmetries::SymmetryConfig, hubbard::HubbardParams{T}) where {T<:AbstractFloat} = CalcConfig(symmetries, hubbard, ())
end
