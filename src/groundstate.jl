##################
# Initialize MPS #
##################

function build_virtualspace(I, charge_ranges::Tuple, maxdim::Int)
    base_key = ntuple(_ -> 0, length(charge_ranges))
    base_key = length(base_key) == 1 ? base_key[1] : base_key

    V = Vect[(I)](base_key => 1)
    for key in Iterators.product(charge_ranges...)
        k = length(key) == 1 ? key[1] : key
        V = V ⊕ Vect[(I)](k => maxdim)
    end
    return V
end

function maximal_virtualspace(::Type{Trivial}, ::Type{Trivial}, total_width::Int, maxdim::Int, P::Int)
    return build_virtualspace(FermionParity, ((0:1),), maxdim)
end
function maximal_virtualspace(::Type{Trivial}, ::Type{U1Irrep}, total_width::Int, maxdim::Int, P::Int)
    return build_virtualspace(FermionParity ⊠ U1Irrep, (0:1, -total_width:1//2:total_width), maxdim)
end
function maximal_virtualspace(::Type{Trivial}, ::Type{SU2Irrep}, total_width::Int, maxdim::Int, P::Int)
    return build_virtualspace(FermionParity ⊠ SU2Irrep, (0:1, 0:1//2:3), maxdim)
end
function maximal_virtualspace(::Type{U1Irrep}, ::Type{Trivial}, total_width::Int, maxdim::Int, P::Int)
    return build_virtualspace(FermionParity ⊠ U1Irrep, (0:1, -total_width*P:total_width*P), maxdim)
end
function maximal_virtualspace(::Type{U1Irrep}, ::Type{U1Irrep}, total_width::Int, maxdim::Int, P::Int)
    return build_virtualspace(FermionParity ⊠ U1Irrep ⊠ U1Irrep,
                        (0:1, -total_width*P:total_width*P, -total_width:1//2:total_width), maxdim)
end
function maximal_virtualspace(::Type{U1Irrep}, ::Type{SU2Irrep}, total_width::Int, maxdim::Int, P::Int)
    return build_virtualspace(FermionParity ⊠ U1Irrep ⊠ SU2Irrep,
                        (0:1, -total_width*P:total_width*P, 0:1//2:3), maxdim)
end
function maximal_virtualspace(::Type{SU2Irrep}, ::Type{Trivial}, total_width::Int, maxdim::Int, P::Int)
    return build_virtualspace(FermionParity ⊠ SU2Irrep, (0:1, 0:1//2:3), maxdim)
end
function maximal_virtualspace(::Type{SU2Irrep}, ::Type{U1Irrep}, total_width::Int, maxdim::Int, P::Int)
    return build_virtualspace(FermionParity ⊠ SU2Irrep ⊠ U1Irrep,
                        (0:1, 0:1//2:3, -total_width:1//2:total_width), maxdim)
end
function maximal_virtualspace(::Type{SU2Irrep}, ::Type{SU2Irrep}, total_width::Int, maxdim::Int, P::Int)
    return build_virtualspace(FermionParity ⊠ SU2Irrep ⊠ SU2Irrep, (0:1, 0:1//2:3, 0:1//2:3), maxdim)
end

function initialize_mps(H::InfiniteMPOHamiltonian, calc::CalcConfig; max_dimension::Int=50)
    Ps = physicalspace.(parent(H))

    # Compute left and right fusion spaces
    V_right = accumulate(fuse, Ps)
    V_left = accumulate(fuse, dual.(Ps); init=one(first(Ps))) |> reverse
    len = length(V_left)
    step = len - 1
    V_left = [view(V_left, len - step + 1:len); view(V_left, 1:len - step)]

    # Intersect left and right spaces
    V = TensorKit.infimum.(V_left, V_right)

    # Construct maximal symmetry-allowed virtual space
    P = calc.symmetries.filling === nothing ? 1 : numerator(calc.symmetries.filling)
    Vmax = maximal_virtualspace(calc.symmetries.particle_symmetry, calc.symmetries.spin_symmetry, length(Ps), max_dimension, P)

    V_trunc = TensorKit.infimum.(V, fill(Vmax, length(V)))

    return InfiniteMPS(Ps, V_trunc)
end

function initialize_mps(H::FiniteMPOHamiltonian, calc::CalcConfig; max_dimension::Int=50, Sz_target=0)
    sym = fℤ₂ ⊠ calc.symmetries.particle_symmetry ⊠ calc.symmetries.spin_symmetry

    use_ps = calc.symmetries.particle_symmetry!=Trivial
    use_ss = calc.symmetries.spin_symmetry!=Trivial
    charges_left = Any[0]
    charges_right = Any[0]

    P = use_ps ? numerator(calc.symmetries.filling) : 1

    if use_ps
        Nphys = Int(calc.symmetries.filling * calc.hubbard.bands * calc.symmetries.cell_width)
        charges_right[1] = mod(Nphys, 2)
        push!(charges_left, 0)
        push!(charges_right, 0)
    end
    if use_ss
        push!(charges_left, 0)
        push!(charges_right, Sz_target)
    end

    left  = Vect[sym](Tuple(charges_left) => 1) 
    right = Vect[sym](Tuple(charges_right) => 1)

    Ps = physicalspace(H)
    Vmax = maximal_virtualspace(calc.symmetries.particle_symmetry, calc.symmetries.spin_symmetry, length(Ps), max_dimension, P)
    return FiniteMPS(Ps, Vmax; left=left, right=right)
end 


####################
# Find groundstate #
####################

"""
    compute_groundstate(calc::CalcConfig;
                        svalue::Float64=2.0,
                        tol::Float64=1e-8,
                        init_state::Union{Nothing, InfiniteMPS, FiniteMPS}=nothing,
                        maxiter::Int=1000,
                        max_init_dim::Int=50,
                        verbosity::Int=0,
                        finite_mps::Bool=false,
                        Sz_target=0,)

Compute the ground state of the Hamiltonian defined by the CalcConfig `calc`.

# Keyword Arguments
- `svalue::Float64=2.0`:
    Exponent controlling the Schmidt value truncation cutoff: tensors with Schmidt
    values below `10^(-svalue)` are discarded.
- `tol::Float64=1e-8`:
    Requested convergence tolerance for the iterative solver (VUMPS or IDMRG2).
    **Important:** `tol` is automatically floored to `schmidtcut / 10`
    (i.e. `10^(-(svalue+1))`), so values tighter than that threshold have no
    effect.  In practice, set `svalue` first and let `tol` follow.
- `init_state::Union{Nothing, InfiniteMPS, FiniteMPS}=nothing`:
    Optional initial MPS. If `finite_mps=false`, expects an `InfiniteMPS`; if `true`, a `FiniteMPS`.
- `maxiter::Int=1000`:
    Maximum number of iterations for ground state optimization.
- `max_init_dim::Int=50`:
    Maximum bond dimension for the initial MPS construction.
- `verbosity::Int=0`:
    Controls the level of printed output from the solver.
- `finite_mps::Bool=false`:
    If `true`, performs a finite-size DMRG calculation with periodic boundary conditions;
    otherwise uses infinite algorithms (VUMPS/IDMRG).
- Sz_target=0,
    Target total spin projection for finite MPS calculations. Only relevant if `finite_mps=true`

# Returns
A `Dict` with the following entries:
- `"groundstate"` → optimized infinite MPS representing the ground state.
- `"environments"` → left and right environment tensors.
- `"ham"` → Hamiltonian MPO used in the optimization.
- `"error"` → final convergence error.
"""
function compute_groundstate(
                calc::CalcConfig;
                svalue::Float64=2.0,
                tol::Float64=1e-8, 
                init_state::Union{Nothing, InfiniteMPS, FiniteMPS}=nothing,
                maxiter::Int64=1000,
                max_init_dim::Int=50,
                verbosity::Int64=0,
                finite_mps::Bool=false,
                Sz_target=0,
            )
    schmidtcut = 10.0^(-svalue)
    tol = max(tol, schmidtcut/10)
    H = hamiltonian(calc)
    total_width = length(H)

    if finite_mps
        H = periodic_boundary_conditions(H, total_width)
        ψ₀ = isnothing(init_state) ? initialize_mps(H, calc; max_dimension=max_init_dim, Sz_target=Sz_target) : init_state
        ψ, envs, δ = find_groundstate(ψ₀, H, DMRG2(; maxiter=maxiter, trscheme=trunctol(; atol=schmidtcut), tol=tol, verbosity=verbosity))
    else
        ψ₀ = isnothing(init_state) ? initialize_mps(H, calc; max_dimension=max_init_dim) : init_state
        if total_width > 1
            ψ₀, envs, = find_groundstate(ψ₀, H, IDMRG2(; maxiter=maxiter, trscheme=trunctol(; atol=schmidtcut), tol=tol, verbosity=verbosity))
        else
            ψ₀, envs, = find_groundstate(ψ₀, H, VUMPS(; maxiter=maxiter, tol=tol, verbosity=verbosity))
            ψ₀ = changebonds(ψ₀, SvdCut(; trscheme=trunctol(; atol=schmidtcut)))
            χ = sum(i -> dim(left_virtualspace(ψ₀, i)), 1:total_width)
            for i in 1:maxiter
                ψ₀, envs = changebonds(ψ₀, H, VUMPSSvdCut(; trscheme=trunctol(; atol=schmidtcut)))
                ψ₀, = find_groundstate(ψ₀, H, VUMPS(; tol=max(tol, schmidtcut / 10), verbosity=verbosity), envs)
                ψ₀ = changebonds(ψ₀, SvdCut(; trscheme=trunctol(; atol=schmidtcut)))
                χ′ = sum(i -> dim(left_virtualspace(ψ₀, i)), 1:total_width)
                isapprox(χ, χ′; rtol=0.05) && break
                χ = χ′
            end
        end
        
        alg = VUMPS(; maxiter=maxiter, tol=tol, verbosity=verbosity) &
            GradientGrassmann(; maxiter=maxiter, tol=tol, verbosity=verbosity)
        ψ, envs, δ = find_groundstate(ψ₀, H, alg)        
    end

    return Dict("groundstate" => ψ, "environments" => envs, "ham" => H, "error" => δ)
end


###########################
# Find chemical potential #
###########################

function change_chemical_potential(model::HubbardParams, μ::Float64)
    bands = model.bands
    t = model.t
    U = model.U

    @assert length(μ) == bands "Length of μ vector must match number of bands."

    t_new = copy(t)
    μ_old = Inf
    for i in 1:bands
        if !haskey(t_new, (i,i))
            t_new[(i,i)] = 0.0
        end
        μ_old = min(t_new[(i,i)], μ_old)
    end

    for i in 1:bands
        # necessary to keep Δ energy constant between bands
        t_new[(i,i)] += μ - μ_old
    end
    
    return HubbardParams(bands, t_new, U)
end

function calculate_filling(calc::CalcConfig, 
                            μ::Float64; 
                            svalue::Float64=2.0, 
                            init_state::Union{Nothing, InfiniteMPS}=nothing
                        )
    simul = CalcConfig(calc.symmetries, change_chemical_potential(calc.hubbard, μ), calc.terms)
    gs = compute_groundstate(simul; svalue=svalue, init_state=init_state)
    ψ = gs["groundstate"]
    filling = sum(density_e(ψ, calc)) / length(ψ)

    return filling, ψ
end

"""
    find_chemical_potential(calc::CalcConfig, 
                            filling::Float64;
                            mu_lower_init::Float64=0.0,
                            mu_upper_init::Float64=1.0,
                            svalue::Float64=2.0,
                            tol_filling::Float64=1e-8, 
                            maxiter::Int64=25,
                            verbosity::Int64=0)
                            
Find the chemical potential μ that yields the desired filling in the ground state.

!!! note
    This function is only implemented for `Trivial` particle symmetry, where the
    filling is not fixed by the symmetry sector.  For `U1Irrep` particle symmetry,
    set the target filling directly via `SymmetryConfig`.
"""
function find_chemical_potential(
                calc::CalcConfig, 
                filling::Float64;
                mu_lower_init::Float64=0.0,
                mu_upper_init::Float64=1.0,
                svalue::Float64=2.0,
                tol_filling::Float64=1e-8, 
                maxiter::Int64=25,
                verbosity::Int64=0
            )

    @assert 0.0 < filling < 2.0 "Filling must be between 0 and 2 electrons per site."
    if calc.symmetries.particle_symmetry != Trivial
        error("Chemical potential search only implemented for trivial particle symmetry.")
    end

    # --- Perform bisection method to find μ iteratively ---

    # Initial bounds
    mu_lower = mu_lower_init
    mu_upper = mu_upper_init
    N_lower, ψ_lower = calculate_filling(calc, mu_lower; svalue=svalue)
    N_upper, ψ_upper = calculate_filling(calc, mu_upper; svalue=svalue)

    if abs(N_upper - filling) < tol_filling
        verbosity > 0 && @info "Converged in 0 iterations. μ = $mu_upper, Filling = $N_upper."
        return mu_upper
    end
    if abs(N_lower - filling) < tol_filling
        verbosity > 0 && @info "Converged in 0 iterations. μ = $mu_lower, Filling = $N_lower."
        return mu_lower
    end

    # Bracketing filling interval
    while (N_lower - filling) * (N_upper - filling) > 0.0
        # If the target filling is not bracketed, expand the interval.
        if N_lower > filling # Both fillings are too high, increase mu_lower
            mu_diff = abs(mu_upper - mu_lower)
            mu_upper = mu_lower
            mu_lower -= mu_diff * 2.0 # Expand downwards
            verbosity > 0 && @info "Target filling $filling not bracketed. Decreasing mu_lower to $mu_lower."
            N_upper = N_lower
            N_lower, ψ_lower = calculate_filling(calc, mu_lower; svalue=svalue, init_state=ψ_lower)
        else # Both fillings are too low, increase mu_upper
            mu_diff = abs(mu_upper - mu_lower)
            mu_lower = mu_upper
            mu_upper += mu_diff * 2.0 # Expand upwards
            verbosity > 0 && @info "Target filling $filling not bracketed. Increasing mu_upper to $mu_upper."
            N_lower = N_upper
            N_upper, ψ_upper = calculate_filling(calc, mu_upper; svalue=svalue, init_state=ψ_upper)
        end
        if abs(N_upper - filling) < tol_filling
            verbosity > 0 && @info "Converged in 0 iterations. μ = $mu_upper, Filling = $N_upper."
            return mu_upper
        end
        if abs(N_lower - filling) < tol_filling
            verbosity > 0 && @info "Converged in 0 iterations. μ = $mu_lower, Filling = $N_lower."
            return mu_lower
        end
        if abs(mu_upper - mu_lower) > 1e3 # Safety break for too large range
             error("Bisection range expanded too much (range [$mu_lower, $mu_upper]). Check initial guess or model parameters.")
        end
    end
    
    verbosity > 0 && @info "Initial bracket: μ ∈ [$mu_lower, $mu_upper], Filling ∈ [$N_lower, $N_upper]"

    # Bisection Iteration
    mu_mid = (mu_lower + mu_upper) / 2.0
    ψ_mid = ψ_lower

    for i in 1:maxiter
        mu_mid = (mu_lower + mu_upper) / 2.0
        N_mid, ψ_mid = calculate_filling(calc, mu_mid; svalue=svalue, init_state=ψ_mid)
        
        # Check for convergence
        fill_diff = abs(N_mid - filling)
        if fill_diff < tol_filling
            verbosity > 0 && @info "Converged in $i iterations. μ = $mu_mid, Filling = $N_mid."
            return mu_mid
        end

        # Bisection logic: Narrow the interval
        if (N_mid - filling) * (N_lower - filling) < 0
            mu_upper = mu_mid
            N_upper = N_mid
        else
            mu_lower = mu_mid
            N_lower = N_mid
        end

        verbosity > 0 && @info "Iter $i: μ ∈ [$mu_lower, $mu_upper], N ∈ [$N_lower, $N_upper]. Diff: $fill_diff"

        # Secondary convergence check based on mu interval width
        if abs(mu_upper - mu_lower) < tol_filling * 10 # Heuristic μ tolerance
            @warn "μ interval width converged (width < $(tol_filling * 10)). μ = $mu_mid, Filling = $N_mid."
            return mu_mid
        end
    end

    @warn "Bisection did not converge in $maxiter iterations. μ ∈ [$mu_lower, $mu_upper], N ∈ [$N_lower, $N_upper]. Final μ = $mu_mid."

    return mu_mid
end