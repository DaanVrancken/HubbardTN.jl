####################
# State properties #
####################

"""
    dim_state(ψ::Union{InfiniteMPS,FiniteMPS})

Determine the bond dimensions in an MPS.
"""
function dim_state(ψ::Union{InfiniteMPS,FiniteMPS})
    dimension = Int64.(zeros(length(ψ)))
    for i in 1:length(ψ)
        dimension[i] = dim(space(ψ.AL[i],1))
    end
    return dimension
end

"""
    density_e(ψ::Union{InfiniteMPS,FiniteMPS}, calc::CalcConfig)

Compute the electron density.
- Returns the density per orbital and per unit-cell position/site.
"""
function density_e(ψ::Union{InfiniteMPS,FiniteMPS}, calc::CalcConfig)
    symm = calc.symmetries
    n = number_e(symm.particle_symmetry, symm.spin_symmetry; filling=symm.filling)

    bands = calc.hubbard.bands

    idx = findfirst(t -> t isa HolsteinTerm, calc.terms)
    w = (idx === nothing ? [] : calc.terms[idx].w)
    boson_modes = (idx === nothing ? 0 : 1) * length(w)

    Ne = zeros(bands, symm.cell_width)
    for i in 1:bands
        for j in 1:symm.cell_width
            site = i + (j - 1) * (bands + boson_modes)
            Ne[i,j] = real(expectation_value(ψ, site => n))
        end
    end
    return Ne
end

"""
    density_b(ψ::Union{InfiniteMPS,FiniteMPS}, calc::CalcConfig)

Compute the number of bosons per site in the unit cell.
"""
function density_b(ψ::Union{InfiniteMPS,FiniteMPS}, calc::CalcConfig)
    symm = calc.symmetries
    idx = findfirst(t -> t isa HolsteinTerm, calc.terms)
    max_b = (idx === nothing ? error("No bosonic terms in model") : calc.terms[idx].max_b)
    w = (idx === nothing ? [] : calc.terms[idx].w)

    n = number_b(symm.particle_symmetry, symm.spin_symmetry, max_b; filling=symm.filling)
    bands = calc.hubbard.bands

    Nb = zeros(length(w), symm.cell_width)
    for i in 1:length(w)
        for j in 1:symm.cell_width
            site = i+(j-1)*(bands+length(w)) + bands
            Nb[i,j] = real(expectation_value(ψ, site => n))
        end
    end
    
    return Nb
end

"""
    density_spin(ψ::Union{InfiniteMPS,FiniteMPS}, calc::CalcConfig)

Compute the electron spin density per site in the unit cell.
"""
function density_spin(ψ::Union{InfiniteMPS,FiniteMPS}, calc::CalcConfig)
    symm = calc.symmetries

    n_up = number_up(symm.particle_symmetry, symm.spin_symmetry; filling=symm.filling)
    n_down = number_down(symm.particle_symmetry, symm.spin_symmetry; filling=symm.filling)

    bands = calc.hubbard.bands

    idx = findfirst(t -> t isa HolsteinTerm, calc.terms)
    w = (idx === nothing ? [] : calc.terms[idx].w)
    boson_modes = (idx === nothing ? 0 : 1) * length(w)
    
    Nup = zeros(bands,symm.cell_width);
    Ndown = zeros(bands,symm.cell_width);
    for i in 1:bands
        for j in 1:symm.cell_width
            site = i+(j-1)*(bands+boson_modes)
            Nup[i,j] = real(expectation_value(ψ, site => n_up))
            Ndown[i,j] = real(expectation_value(ψ, site => n_down))
        end
    end

    return Nup, Ndown
end

"""
    calc_ms(ψ::InfiniteMPS, calc::CalcConfig)

Compute the staggered magnetization in an InfiniteMPS.
"""
function calc_ms(ψ::InfiniteMPS, calc::CalcConfig)
    up, down = density_spin(ψ, calc)
    Mag = (up - down)/2
    if !all(x -> isapprox(abs(x),abs(Mag[1,1]),rtol=10^(-6)), vec(Mag))
        @warn "Staggerd magnetization varies across unit cell: returning value for first site only."
    end
    return abs(Mag[1,1])
end


"""
    get_alpha(ψ::InfiniteMPS, calc::CalcConfig, ty::T, tz::T, Ep::T) where {T<:AbstractFloat}

Compute α-coefficients from pair correlators.

Here:
- `ty` is the hopping between neighboring ladders within the same plane,
- `tz` is the hopping to the ladders below and above the ladder under consideration.

For a 1-band model, returns `[a0, a01]`, where:
- `a0`  is the onsite contribution from ⟨c↓ c↑⟩,
- `a01` is the nearest-neighbor contribution from ⟨c↑₁ c↓₂⟩.

For a 2-band model, returns `[a0, a1, a00, a01, a10, a11]`.

Only onsite and nearest-neighbor terms are included.
"""
function get_alpha(ψ::InfiniteMPS, calc::CalcConfig, ty::T, tz::T, Ep::T) where {T<:AbstractFloat}
    symm  = calc.symmetries
    ps    = symm.particle_symmetry
    ss    = symm.spin_symmetry
    bands = calc.hubbard.bands

    @assert Ep != 0 "Ep must be nonzero"

    # onsite pair annihilation Δ = c↓ c↑
    Δ   = delete_pair_onesite(ps, ss)

    if bands == 1
        c0  = real(expectation_value(ψ, 1 => Δ))
        c01 = real(expectation_value(ψ, (1,2) => HubbardOperators.u_min_d_min(ComplexF64, ps, ss)))

        a01 = 2 * 4 * ty * tz * c01 / Ep
        a0  = 2 * 4 * ty * tz * c0  / Ep

        return [a0, a01]
    elseif bands == 2
        c0 = real(expectation_value(ψ, 1 => Δ))
        c1 = real(expectation_value(ψ, 2 => Δ))
        c00 = real(expectation_value(ψ, (1,3) => HubbardOperators.u_min_d_min(ComplexF64, ps, ss)))
        c01 = real(expectation_value(ψ, (1,2) => HubbardOperators.u_min_d_min(ComplexF64, ps, ss)))
        c11 = real(expectation_value(ψ, (2,4) => HubbardOperators.u_min_d_min(ComplexF64, ps, ss)))

        a00 = 2 * (ty^2 * c11 + 2 * tz * c00) / Ep
        a11 = 2 * (ty^2 * c00 + 2 * tz * c11) / Ep
        a01 = 4 * tz^2 * c01 / Ep
        a0 = 2 * (ty^2 * c1 + 2 * tz^2 * c0) / Ep
        a1 = 2 * (ty^2 * c0 + 2 * tz^2 * c1) / Ep

        return [a0, a1, a00, a01, a01, a11]
    else
        error("get_alpha is only implemented for 1-band and 2-band models, got bands = $bands")
    end
end

"""
    get_beta(ψ::InfiniteMPS, calc::CalcConfig, ty::T, tz::T, Ep::T) where {T<:AbstractFloat}

Compute β-coefficients from density and hopping correlators.

Here:
- `ty` is the hopping between neighboring ladders within the same plane,
- `tz` is the hopping to the ladders below and above the ladder under consideration.

For a 1-band model, returns `[b0, b01]`, where:
- `b0`  is the onsite contribution from ⟨n⟩,
- `b01` is the nearest-neighbor contribution from ⟨c†₁ c₂⟩.

For a 2-band model, returns `[b00, b01, b10, b11]`.

Only onsite and nearest-neighbor terms are included.
"""
function get_beta(ψ::InfiniteMPS, calc::CalcConfig, ty::T, tz::T, E::T) where {T<:AbstractFloat}
    symm  = calc.symmetries
    ps    = symm.particle_symmetry
    ss    = symm.spin_symmetry
    bands = calc.hubbard.bands

    @assert E != 0 "E must be nonzero"

    if bands == 1
        cu  = real(expectation_value(ψ, (1,1) => c_plusmin_up(ComplexF64, ps, ss)))
        cd  = real(expectation_value(ψ, (1,1) => c_plusmin_down(ComplexF64, ps, ss)))
        c01 = real(expectation_value(ψ, (1,2) => c_plusmin_up(ComplexF64, ps, ss)))

        println("cu = ", cu)
        println("cd = ", cd)
        println("c01 = ", c01)
        
        b0_u  = 2 * 4 * ty * tz * cu  / E
        b0_d  = 2 * 4 * ty * tz * cd  / E
        b01 = 2 * 4 * ty * tz * c01 / E

        if ss == U1Irrep
            return [b0_u, b0_d, b01]
        end

        c0_ud  = real(expectation_value(ψ, (1,1) => c_plusmin_updown(ComplexF64, ps, ss)))
        c01_ud  = real(expectation_value(ψ, (1,2) => c_plusmin_updown(ComplexF64, ps, ss)))

        println("c0_ud = ", c0_ud)
        println("c01_ud = ", c01_ud)         

        b0_ud  = 2 * 4 * ty * tz * c0_ud  / E
        b01_ud = 2 * 4 * ty * tz * c01_ud / E

        return [b0_u, b0_d, b01, b0_ud, b01_ud]

    elseif bands == 2
        c0_u  = real(expectation_value(ψ, (1,1) => c_plusmin_up(ComplexF64, ps, ss)))
        c0_d  = real(expectation_value(ψ, (1,1) => c_plusmin_down(ComplexF64, ps, ss)))
        println("c0_u = ", c0_u)
        println("c0_d = ", c0_d)

        b0_u = 2 * (ty^2 * c0_d + 2 * tz^2 * c0_u) / E
        b0_d = 2 * (ty^2 * c0_u + 2 * tz^2 * c0_d) / E

        c00 = real(expectation_value(ψ, (1,3) => c_plusmin_up(ComplexF64, ps, ss)))
        c01 = real(expectation_value(ψ, (1,2) => c_plusmin_up(ComplexF64, ps, ss)))
        c10 = real(expectation_value(ψ, (2,1) => c_plusmin_up(ComplexF64, ps, ss)))
        c11 = real(expectation_value(ψ, (2,4) => c_plusmin_up(ComplexF64, ps, ss)))

        println("c00 = ", c00)
        println("c01 = ", c01)
        println("c10 = ", c10)
        println("c11 = ", c11)

        b00 = 2 * (ty^2 * c11 + 2 * tz^2 * c00) / E
        b01 = (4 * tz^2 * c10) / E
        b10 = (4 * tz^2 * c01) / E
        b11 = 2 * (ty^2 * c00 + 2 * tz^2 * c11) / E

        if ss == U1Irrep
            return [b0_u, b0_d, b00, b01, b10, b11]
        end

        c00_ud = real(expectation_value(ψ, (1,3) => c_plusmin_updown(ComplexF64, ps, ss)))
        c01_ud  = real(expectation_value(ψ, (1,2) => c_plusmin_updown(ComplexF64, ps, ss)))
        c10_ud  = real(expectation_value(ψ, (2,1) => c_plusmin_updown(ComplexF64, ps, ss)))
        c11_ud  = real(expectation_value(ψ, (2,4) => c_plusmin_updown(ComplexF64, ps, ss)))

        println("c00_ud = ", c00_ud)
        println("c01_ud = ", c01_ud)
        println("c10_ud = ", c10_ud)
        println("c11_ud = ", c11_ud)        

        b00_ud = 2 * (ty^2 * c11_ud + 2 * tz^2 * c00_ud) / E
        b01_ud = (4 * tz^2 * c10_ud) / E
        b10_ud = (4 * tz^2 * c01_ud) / E
        b11_ud = 2 * (ty^2 * c00_ud + 2 * tz^2 * c11_ud) / E

        return [b0_u, b0_d, b00, b01, b10, b11, b00_ud, b01_ud, b10_ud, b11_ud]
    else
        error("get_beta is only implemented for 1-band and 2-band models, got bands = $bands")
    end
end

"""
    density_correlations(ψ::InfiniteMPS, calc::CalcConfig; R::Int=15, thr::Float64=1e-10)

Compute connected density–density correlations

    C(r) = ⟨n₀ n_r⟩ − ⟨n⟩²

for distances `r = 1:R` in the infinite MPS `ψ`.

The function prints the correlation values and stops early if the
magnitude falls below `thr`, which helps avoid unnecessary evaluations
once correlations have effectively decayed.

The site spacing accounts for interleaved boson modes and electronic
bands defined in `calc`. The resulting correlations can later be used
to estimate the charge Luttinger parameter `Kρ` from the small-q
behavior of the density structure factor.
"""
function density_correlations(ψ::InfiniteMPS, calc::CalcConfig; R::Int=15, thr::Float64=1e-10)

    symm = calc.symmetries
    ps   = symm.particle_symmetry
    ss   = symm.spin_symmetry

    # --- Kρ estimate from density structure factor at small q ---
    idx = findfirst(t -> t isa HolsteinTerm, calc.terms)
    w = (idx === nothing ? [] : calc.terms[idx].w)
    boson_modes = (idx === nothing ? 0 : 1) * length(w)
    bands = calc.hubbard.bands

    n = number_e(ps, ss)

    # connected density correlator C(r) = <n0 n_r> - <n>^2
    C = zeros(Float64, R)
    nn = n ⊗ n
    for r in 1:R
        sr = 1 + (r-1) * (boson_modes + bands)
        C[r] = real(expectation_value(ψ, (1, sr) => nn) - (expectation_value(ψ, (1) => n) * expectation_value(ψ, (sr) => n)))
        println("r=$(r-1)  C[r]=$(C[r])")

        if abs(C[r]) < thr
            println("Correlation below threshold at r=$r → stopping.")
            C = C[1:r]   
        break
        end
    end
end


##########
# Saving #
##########

"""
    save_computation(d::Dict{String, Any}, path::String, file_name::String)

Save the output dictionary of e.g. `compute_groundstate` as a `.jld2` file at the specified path.
"""
function save_computation(d::Dict{String, Any}, path::String, file_name::String)
    ispath(path) || mkdir(path)
    @save joinpath(path, file_name*".jld2") d
end

"""
    load_computation(path_to_file::String)

Load the output dictionary of e.g. `compute_groundstate` stored as a `.jld2` file.
"""
function load_computation(path_to_file::String)
    @load joinpath(path_to_file) d
    return d
end

"""
    save_state(ψ::InfiniteMPS, path::String, name::String)

Save the tensors of an `InfiniteMPS` object to disk as individual `.jld2` files.

# Arguments
- `ψ::InfiniteMPS`: The infinite matrix product state (MPS) whose tensors will be saved.
- `path::String`: The base directory where the state folder will be created.
- `name::String`: The name of the subdirectory under `path` where the tensors will be stored.

# Description
This function creates a subdirectory `joinpath(path, name)` and saves each tensor 
`ψ.AL[i]` as a `.jld2` file named `state<i>.jld2` inside it. Each tensor is converted 
to a `Dict` before saving for serialization compatibility. The function prints 
a message after each tensor is successfully saved. This approach may be useful for
storing large MPS objects.
"""
function save_state(ψ::InfiniteMPS, path::String, name::String)
    path = joinpath(path,name)
    mkdir(path)
    for i in 1:length(ψ)
        d = convert(Dict,ψ.AL[i])
        @save joinpath(path,"state$i.jld2") d
        @info "State $i saved."
    end
end

"""
    load_state(path::String) -> InfiniteMPS

Load an `InfiniteMPS` object from a directory of saved `.jld2` tensor files.

# Arguments
- `path::String`: Path to the directory containing the saved MPS tensor files 
  (e.g., `state1.jld2`, `state2.jld2`, ...).

# Description
This function reconstructs an `InfiniteMPS` object previously saved with [`save_state`](@ref).  
It reads all `.jld2` files in the specified directory, converts each stored `Dict`
back into a `TensorMap`, and combines them into a periodic array before wrapping 
the result in an `InfiniteMPS` object.

# Returns
- `InfiniteMPS`: The reconstructed infinite matrix product state.
"""
function load_state(path::String)
    entries = readdir(path)
    file_count = count(entry -> isfile(joinpath(path, entry)), entries)

    @load joinpath(path,"state1.jld2") d
    A = [convert(TensorMap, d)]
    for i in 2:file_count
        @load joinpath(path,"state$i.jld2") d
        push!(A, convert(TensorMap, d))
    end

    return InfiniteMPS(PeriodicArray(A))
end

compact_float(x::Real) = replace(rstrip(rstrip(string(x), '0'), '.'), "-0" => "0")

"""
    dict_tag(d::Dict)

Returns a canonical string tag from a hopping/interaction dictionary for saving.
"""
function dict_tag(d::Dict)
    pairs = sort(collect(d); by = first)
    parts = String[]

    for (k, v) in pairs
        (v isa AbstractFloat && iszero(v)) && continue
        if k isa Tuple && length(k) == 2
            # hopping: t_ij, keep i<j
            i, j = k
            i > j && continue
            push!(parts, "t$(i)$(j)_" * compact_float(v))

        elseif k isa Tuple && length(k) == 4
            # interaction: U_ijkl, canonicalize symmetry
            i, j, l, m = k
            if (i, j, l, m) > (j, i, m, l)
                i, j, l, m = j, i, m, l
            end
            if (i, j) > (l, m)
                i, j, l, m = l, m, i, j
            end
            push!(parts, "U$(i)$(j)$(l)$(m)_" * compact_float(v))
        else
            error("Unsupported key in dict_tag: $k")
        end
    end

    unique!(parts)
    return isempty(parts) ? "0" : join(parts, "_")
end