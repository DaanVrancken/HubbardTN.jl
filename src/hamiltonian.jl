############################
# Hamiltonian construction #
############################

# Build the symmetry-dependent operators
function build_ops(symm::SymmetryConfig, bands::Int64, max_b::Int64, nmodes::Int64)
    ps = symm.particle_symmetry
    ss = symm.spin_symmetry
    fill = symm.filling

    electron_space = hubbard_space(ps, ss; filling=fill)

    ops = (
        c⁺c      = c_plusmin(ps, ss; filling=fill),
        n        = number_e(ps, ss; filling=fill),
        c⁺c⁺cc   = two_body(ps, ss; filling=fill),
        c⁺c⁺c⁺ccc = three_body(ps, ss; filling=fill)
    )
    if ss !== SU2Irrep
        ops = merge(ops, (Sz = Sz(ps, ss; filling=fill),))
    end
    if ss === Trivial
        ops = merge(ops, (Sx = Sx(ps, ss; filling=fill), Sy = Sy(ps, ss; filling=fill)))
        ops = merge(ops, (c⁺c_ud = c_plusmin_updown(ps, ss; filling=fill), c⁺c_du = c_plusmin_downup(ps, ss; filling=fill)))
        ops = merge(ops, (c⁺c_uu = c_plusmin_up(ps, ss; filling=fill), c⁺c_dd = c_plusmin_down(ps, ss; filling=fill)))
    end
    if ps === Trivial
        ops = merge(ops, (c⁺pair = create_pair_onesite(ps, ss; filling=fill), cpair = delete_pair_onesite(ps, ss; filling=fill)))
    end

    phonon_spaces = []
    if max_b > 0
        ops = merge(ops, (bmin = b_min(ps, ss, max_b; filling=fill),
                          bplus = b_plus(ps, ss, max_b; filling=fill),
                          nb = number_b(ps, ss, max_b; filling=fill)))

        phonon_space = boson_space(ps, ss, max_b; filling=fill)
        phonon_spaces = [phonon_space for _ in 1:nmodes]
    end

    electron_spaces = [electron_space for _ in 1:bands]
    spaces = append!(electron_spaces, phonon_spaces)

    return ops, repeat(spaces, symm.cell_width)
end

"""
    hamiltonian(calc::CalcConfig)

Constructs the many-body Hamiltonian for a system defined by configuration `calc`.

# Notes
- Lattice sites are represented using an `InfiniteChain` of length `cell_width * bands`.
- The resulting MPO can be used directly for DMRG, VUMPS, or other tensor network calculations.
"""
function hamiltonian(calc::CalcConfig{T}) where {T<:AbstractFloat}
    empty!(two_body_cache)
    empty!(three_body_cache)

    bands = calc.hubbard.bands
    t = calc.hubbard.t
    U = calc.hubbard.U

    idx = findfirst(t -> t isa HolsteinTerm, calc.terms)
    max_b = (idx === nothing ? 0 : calc.terms[idx].max_b)
    w = (idx === nothing ? [] : calc.terms[idx].w)
    boson_modes = Int(max_b>0) * length(w)
    period = bands + boson_modes

    ops, spaces = build_ops(calc.symmetries, bands, max_b, boson_modes)
    cell_width = calc.symmetries.cell_width

    h::Vector{Pair{Tuple{Vararg{Int64}}, Any}} = [(1,) => 0*ops.n]  # Initialize MPO

    # --- Hopping ---
    for cell in 0:(cell_width-1)
        site(i) = i + cell*period + div(i-1, bands)*boson_modes
        h = append!(h, [site.((i,j)) => -t_ij*ops.c⁺c for ((i,j), t_ij) in t])
    end

    # --- 2-body Interaction ---
    for cell in 0:(cell_width-1)
        site(i) = i + cell*period + div(i-1, bands)*boson_modes
        h = append!(h, [site.((i,j,k,l)) => 0.5 * U_ijkl * ops.c⁺c⁺cc for ((i,j,k,l), U_ijkl) in U])
    end

    H = InfiniteMPOHamiltonian(spaces, h...)

    # --- Extra terms ---
    for term in calc.terms
        H += hamiltonian_term(term, ops, spaces, cell_width, bands, boson_modes)
    end

    return H
end

# Three-body interaction term
function hamiltonian_term(
                    term::ThreeBodyTerm, 
                    ops, 
                    spaces,
                    cell_width::Int64,
                    bands::Int64,
                    boson_modes::Int64
                )
    V = term.V
    period = bands + boson_modes

    h = []

    for cell in 0:(cell_width-1)
        site(i) = i + cell*period + div(i-1, bands)*boson_modes
        h = append!(h, [site.((i,j,k,l,m,n)) => 1/6 * V_ijklmn * ops.c⁺c⁺c⁺ccc for ((i,j,k,l,m,n), V_ijklmn) in V])
    end

    return InfiniteMPOHamiltonian(spaces, h...)
end
# Magnetic field term
function hamiltonian_term(
                    term::MagneticField, 
                    ops, 
                    spaces,
                    cell_width::Int64,
                    bands::Int64,
                    boson_modes::Int64
                )
    B = term.B

    electron_sites = [i + div(i-1, bands)*boson_modes for i in 1:(cell_width*bands)]

    return InfiniteMPOHamiltonian(spaces, [(i,) => -B*ops.Sz for i in electron_sites]...)
end
# Staggered magnetic field term
function hamiltonian_term(
                    term::StaggeredField, 
                    ops, 
                    spaces,
                    cell_width::Int64,
                    bands::Int64,
                    boson_modes::Int64
                )
    J = term.J
    Ms = term.Ms
    period = bands + boson_modes
    phase = (-1) .^ (div.(0:(period*cell_width-1), period))

    electron_sites = [i + div(i-1, bands)*boson_modes for i in 1:(cell_width*bands)]

    return InfiniteMPOHamiltonian(spaces, [(i,) => 2*J*Ms * phase[i] * ops.Sz for i in electron_sites]...)
end
# Spin mean field term
function hamiltonian_term(
                    term::SpinMeanField, 
                    ops,
                    spaces, 
                    cell_width::Int64,
                    bands::Int64,
                    boson_modes::Int64
                )
    J = term.J
    s = term.spins

    electron_sites = [i + div(i-1, bands)*boson_modes for i in 1:(cell_width*bands)]

    if length(size(s)) == 1
        h = [(i,) => J[i,j]*s[j]*ops.Sz for i in electron_sites, j in electron_sites]
    else
        h = [(i,) => J[i,j]*(s[j,1]*ops.Sx + s[j,2]*ops.Sy + s[j,3]*ops.Sz) for i in electron_sites, j in electron_sites]
    end
    return InfiniteMPOHamiltonian(spaces, h...)
end
# Charge gap mean field term
function hamiltonian_term(
                    term::ChargeGapMF, 
                    ops,
                    spaces,
                    cell_width::Int64,
                    bands::Int64,
                    boson_modes::Int64
                )
    hasproperty(ops, :c⁺c_uu) || throw(ArgumentError("ChargeGapMF requires Trivial spin symmetry."))

    electron_site(i) = 1 + fld(i - 1, bands) * (bands + boson_modes) + mod(i - 1, bands)
    beta_index(i)    = mod1(i, bands*cell_width)

    t_inter = term.t_inter
    range   = term.range

    h = Any[]
    for cell in 0:(cell_width-1), ((i, j), t_ij) in t_inter, ((k, l), t_kl) in t_inter, r in -range*bands:bands:range*bands

        (abs(i - k - r) <= range && abs(j - l - r) <= range) || continue
        
        coefficient = 2 * t_ij * t_kl
        sites   = (electron_site(i + cell*bands), electron_site(k + r + cell*bands))
        idx = (beta_index(j + cell*bands), beta_index(l + r + cell*bands))
        append!(h, [
            sites => coefficient * term.beta_uu[idx...]  * ops.c⁺c_uu,
            sites => coefficient * term.beta_ud[idx...]  * ops.c⁺c_ud,
            sites => coefficient * term.beta_ud'[idx...] * ops.c⁺c_du,
            sites => coefficient * term.beta_uu'[idx...] * ops.c⁺c_dd
        ])
    end 

    return InfiniteMPOHamiltonian(spaces, h...)
end
# Pair gap mean field term
function hamiltonian_term(
                    term::PairGapMF, 
                    ops,
                    spaces,
                    cell_width::Int64,
                    bands::Int64,
                    boson_modes::Int64
                )

    electron_sites = [i + div(i-1, bands)*boson_modes for i in 1:(cell_width*bands)]

    if hasproperty(ops, :cpair)
        hopping_onsite = ops.c⁺pair + ops.cpair
        hopping_pair = HubbardOperators.d_plus_u_plus(ComplexF64,Trivial,U1Irrep) + HubbardOperators.u_min_d_min(ComplexF64,Trivial,U1Irrep)
    end

    if bands == 1
        a0, a01 = term.alpha
        b0, b01 = term.beta
    elseif bands == 2
        a0, a1, a00, a01, a10, a11 = term.alpha
        if hasproperty(ops, :c⁺c_ud)
            b00, b01, b10, b11, b00_ud, b01_ud, b10_ud, b11_ud = term.beta
        else
            b00, b01, b10, b11 = term.beta
        end
    else
        error("PairGapMF term: only 1-band and 2-band models are implemented, got bands = $bands.")
    end
    
    h = Any[]

    if bands ==1
        if hasproperty(ops, :cpair)
            append!(h, [
                (i,) => -a0 * hopping_onsite
                for i in electron_sites
            ])

            h = append!(h, [
                (electron_sites[n], electron_sites[n+1]) => -a01*hopping_pair
                for n in 1:(length(electron_sites)-1)
            ])
            h = append!(h, [
                (electron_sites[n+1], electron_sites[n]) => -a01*hopping_pair
                for n in 1:(length(electron_sites)-1)
            ])
        end
        h = append!(h, [
            (electron_sites[n+1], electron_sites[n]) => b01*ops.c⁺c
            for n in 1:(length(electron_sites)-1)
        ])
        h = append!(h, [
            (electron_sites[n], electron_sites[n+1]) => b01*ops.c⁺c
            for n in 1:(length(electron_sites)-1)
        ])
        return InfiniteMPOHamiltonian(spaces, h...)
    end

    if bands == 2
        if hasproperty(ops, :cpair)
            @assert a0 == a1 
            append!(h, [
                (i,) => -a0 * hopping_onsite
                for i in electron_sites
            ])
            @assert a01 == a10
            h = append!(h, [(1, 2) => -a01*hopping_pair])
            h = append!(h, [(2, 1) => -a01*hopping_pair])
            h = append!(h, [(1, 3) => -a00*hopping_pair])
            h = append!(h, [(3, 1) => -a00*hopping_pair])
            h = append!(h, [(2, 4) => -a11*hopping_pair])
            h = append!(h, [(4, 2) => -a11*hopping_pair])
        end
        @assert b01 == b10
        h = append!(h, [(1, 2) => b01*ops.c⁺c])
        h = append!(h, [(2, 1) => b01*ops.c⁺c])
        h = append!(h, [(1, 3) => b00*ops.c⁺c])
        h = append!(h, [(3, 1) => b00*ops.c⁺c])
        h = append!(h, [(2, 4) => b11*ops.c⁺c])
        h = append!(h, [(4, 2) => b11*ops.c⁺c])

        if hasproperty(ops, :c⁺c_ud)
            @assert b01_ud == b10_ud
            h = append!(h, [(1, 2) => b01_ud*ops.c⁺c_ud + b01_ud*ops.c⁺c_du])
            h = append!(h, [(2, 1) => b01_ud*ops.c⁺c_ud + b01_ud*ops.c⁺c_du])
            h = append!(h, [(1, 3) => b00_ud*ops.c⁺c_ud + b00_ud*ops.c⁺c_du])
            h = append!(h, [(3, 1) => b00_ud*ops.c⁺c_ud + b00_ud*ops.c⁺c_du])
            h = append!(h, [(2, 4) => b11_ud*ops.c⁺c_ud + b11_ud*ops.c⁺c_du])
            h = append!(h, [(4, 2) => b11_ud*ops.c⁺c_ud + b11_ud*ops.c⁺c_du])
        end

        return InfiniteMPOHamiltonian(spaces, h...)
    end
end
# Holstein coupling term
function hamiltonian_term(
                    term::HolsteinTerm, 
                    ops,
                    spaces, 
                    cell_width::Int64,
                    bands::Int64,
                    boson_modes::Int64
                )
    w = term.w
    g = term.g
    mean_ne = term.mean_ne
    xi = term.xi

    period = bands + boson_modes

    electron_sites = [i + div(i-1, bands)*boson_modes for i in 1:(cell_width*bands)]
    electron_ind(i) = mod1(i, period)
    phonon_sites = [i + bands + div(i-1, boson_modes)*bands for i in 1:(cell_width*boson_modes)]
    phonon_ind(i) = mod1(i, period) - bands
    cell(i) = div(i-1, period)

    H_ph = InfiniteMPOHamiltonian(spaces, [(i,) => w[phonon_ind(i)] * ops.nb for i in phonon_sites]...)

    H_ep = 0 * H_ph

    # Precompute non-local exponential fit for a power-law
    if xi != Inf
        K = 1
        cs, λs, err = inv_power_expsum(xi, K)

        while err ≥ term.threshold
            K += 1
            cs, λs, err = inv_power_expsum(xi, K)
        end

        cs = real.(cs)
        cs ./= sum(cs)
        λs = real.(λs)

        @info "Created exponential fit for non-local Holstein coupling: K=$K err=$err"
    end

    for e in electron_sites
        ce = cell(e)
        be = electron_ind(e)
        for p in phonon_sites
            cp = cell(p)
            m = phonon_ind(p)
            O_e = g[be, m] * (ops.n - mean_ne * id(domain(ops.n)))
            O_p = ops.bmin + ops.bplus
            O_ep = O_e ⊗ O_p

            if xi == Inf # Pure local Holstein coupling
                if ce == cp
                    H_ep += InfiniteMPOHamiltonian(spaces, (e, p) => O_ep)
                end
            else # Nonlocal Holstein coupling in terms of exponentials
                if ce == cp
                    println(e,p,g[be,m])
                    for (c, λ) in zip(cs, λs)
                        H_ep += exponential_mpo(spaces, (e, p), c * O_ep, λ^2)
                    end

                elseif abs(ce - cp) == 1
                    println(e,p,g[be,m])
                    for (c, λ) in zip(cs, λs)
                        H_ep += exponential_mpo(spaces, (e, p), c * λ * O_ep, λ^2)
                    end
                end
            end
        end
    end

    return H_ph + H_ep
end
