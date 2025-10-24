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

function initialize_mps(symm::SymmetryConfig, total_width::Int; max_dimension::Int=50)
    ps = hubbard_space(symm.particle_symmetry, symm.spin_symmetry; filling=symm.filling)
    Ps = [ps for _ in 1:total_width]

    # Compute left and right fusion spaces
    V_right = accumulate(fuse, Ps)
    V_left = accumulate(fuse, dual.(Ps); init=one(first(Ps))) |> reverse
    len = length(V_left)
    step = len - 1
    V_left = [view(V_left, len - step + 1:len); view(V_left, 1:len - step)]

    # Intersect left and right spaces
    V = TensorKit.infimum.(V_left, V_right)

    # Construct maximal symmetry-allowed virtual space
    P = symm.filling === nothing ? 1 : symm.filling[1]
    Vmax = maximal_virtualspace(symm.particle_symmetry, symm.spin_symmetry, total_width, max_dimension, P)

    V_trunc = TensorKit.infimum.(V, fill(Vmax, length(V)))

    return InfiniteMPS(Ps, V_trunc)
end


####################
# Find groundstate #
####################

"""
    compute_groundstate(calc::CalcConfig;
                        svalue::Float64=2.0,
                        tol::Float64=1e-8,
                        init_state::Union{Nothing, InfiniteMPS}=nothing,
                        maxiter::Int=1000,
                        max_init_dim::Int=50,
                        verbosity::Int=0)

Compute the ground state of the Hamiltonian defined by the CalcConfig `calc`.

# Keyword Arguments
- `svalue::Float64=2.0`: 
    Exponent used to define the truncation cutoff as `10^(-svalue)` for Schmidt value truncation.
- `tol::Float64=1e-8`: 
    Convergence tolerance for iterative solvers (used by VUMPS or IDMRG2).
- `init_state::Union{Nothing, InfiniteMPS}=nothing`: 
    Optional initial infinite MPS. If not provided, a random symmetry-consistent MPS.
- `maxiter::Int=1000`: 
    Maximum number of iterations for ground state optimization.
- `max_init_dim::Int=50`: 
    Maximum bond dimension for the initial MPS construction.
- `verbosity::Int=0`: 
    Controls the level of printed output from the solver.

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
                init_state::Union{Nothing, InfiniteMPS}=nothing,
                maxiter::Int64=1000,
                max_init_dim::Int=50,
                verbosity::Int64=0
            )
    H = hamiltonian(calc)

    symm = calc.symmetries
    total_width = calc.model.bands * symm.cell_width
    ψ₀ = isnothing(init_state) ? initialize_mps(symm, total_width; max_dimension=max_init_dim) : init_state

    schmidtcut = 10.0^(-svalue)
    tol = max(tol, schmidtcut/10)
    
    if total_width > 1
        ψ₀, envs, = find_groundstate(ψ₀, H, IDMRG2(; maxiter=maxiter, trscheme=truncbelow(schmidtcut), tol=tol, verbosity=verbosity))
    else
        ψ₀, envs, = find_groundstate(ψ₀, H, VUMPS(; maxiter=maxiter, tol=tol, verbosity=verbosity))
        ψ₀ = changebonds(ψ₀, SvdCut(; trscheme=truncbelow(schmidtcut)))
        χ = sum(i -> dim(left_virtualspace(ψ₀, i)), 1:total_width)
        for i in 1:maxiter
            ψ₀, envs = changebonds(ψ₀, H, VUMPSSvdCut(;trscheme=truncbelow(schmidtcut)))
            ψ₀, = find_groundstate(ψ₀, H, VUMPS(; tol=max(tol, schmidtcut / 10), verbosity=verbosity), envs)
            ψ₀ = changebonds(ψ₀, SvdCut(; trscheme=truncbelow(schmidtcut)))
            χ′ = sum(i -> dim(left_virtualspace(ψ₀, i)), 1:total_width)
            isapprox(χ, χ′; rtol=0.05) && break
            χ = χ′
        end
    end
    
    alg = VUMPS(; maxiter=maxiter, tol=tol, verbosity=verbosity) &
        GradientGrassmann(; maxiter=maxiter, tol=tol, verbosity=verbosity)
    ψ, envs, δ = find_groundstate(ψ₀, H, alg)
    
    return Dict("groundstate" => ψ, "environments" => envs, "ham" => H, "error" => δ)
end

"""
    find_chemical_potential(calc::CalcConfig;
                            svalue::Float64=2.0,
                            tol::Float64=1e-8,
                            init_state::Union{Nothing, InfiniteMPS}=nothing,
                            maxiter::Int=1000,
                            max_init_dim::Int=50,
                            verbosity::Int=0)
                            
Find the chemical potential μ that yields the desired filling in the ground state.
"""
function find_chemical_potential(
                calc::CalcConfig;
                svalue::Float64=2.0,
                tol::Float64=1e-8, 
                init_state::Union{Nothing, InfiniteMPS}=nothing,
                maxiter::Int64=1000,
                max_init_dim::Int=50,
                verbosity::Int64=0
            )
    error("Not yet implemnted.")
    verbosity_mu = get(simul.kwargs, :verbosity_mu, 0)
    t = simul.t
    u = simul.u
    s = simul.svalue
    bond_dim=simul.bond_dim 
    period = simul.period
    kwargs = simul.kwargs

    if simul.μ !== nothing
        simul2 = OBC_Sim2(t,u,simul.μ,s,bond_dim,period;kwargs)
        dictionary = compute_groundstate(simul2; tol=tol, verbosity=verbosity, maxiter=maxiter);
        dictionary["μ"] = simul.μ
    else 
        f = simul.f
        tol_mu = get(kwargs, :tol_mu, 1e-8)
        maxiter_mu = get(kwargs, :maxiter_mu, 20)
        step_size = get(kwargs, :step_size, 1.0)
        flag = false

        lower_bound = get(simul.kwargs, :lower_mu, 0.0)
        upper_bound = get(simul.kwargs, :upper_mu, 0.0)
        mid_point = (lower_bound + upper_bound)/2
        i = 1

        simul2 = OBC_Sim2(t,u,lower_bound,s,bond_dim,period;kwargs)
        dictionary_l = compute_groundstate(simul2; tol=tol, verbosity=verbosity, maxiter=maxiter);
        dictionary_u = deepcopy(dictionary_l)
        dictionary_sp = deepcopy(dictionary_l)
        while i<=maxiter_mu
            if abs(density_state(dictionary_u["groundstate"]) - f) < tol_mu
                flag=true
                dictionary_sp = deepcopy(dictionary_u)
                mid_point = upper_bound
                break
            elseif abs(density_state(dictionary_l["groundstate"]) - f) < tol_mu
                flag=true
                dictionary_sp = deepcopy(dictionary_l)
                mid_point = lower_bound
                break
            elseif density_state(dictionary_u["groundstate"]) < f
                lower_bound = copy(upper_bound)
                upper_bound += step_size
                simul2 = OBC_Sim2(t,u,upper_bound,s,bond_dim,period;kwargs)
                dictionary_u = compute_groundstate(simul2; tol=tol, verbosity=verbosity, maxiter=maxiter)
            elseif density_state(dictionary_l["groundstate"]) > f
                upper_bound = copy(lower_bound)
                lower_bound -= step_size
                simul2 = OBC_Sim2(t,u,lower_bound,s,bond_dim,period;kwargs)
                dictionary_l = compute_groundstate(simul2; tol=tol, verbosity=verbosity, maxiter=maxiter)
            else
                break
            end
            verbosity_mu>0 && @info "Iteration μ: $i => Lower bound: $lower_bound; Upper bound: $upper_bound"
            i+=1
        end
        if upper_bound>0.0
            value = "larger"
            dictionary = dictionary_u
        else
            value = "smaller"
            dictionary = dictionary_l
        end
        if i>maxiter_mu
            max_value = (i-1)*step_size
            @warn "The chemical potential is $value than: $max_value. Increase the stepsize."
        end

        while abs(density_state(dictionary["groundstate"]) - f)>tol_mu && i<=maxiter_mu && !flag
            mid_point = (lower_bound + upper_bound)/2
            simul2 = OBC_Sim2(t,u,mid_point,s,bond_dim,period;kwargs)
            dictionary = compute_groundstate(simul2)
            if density_state(dictionary["groundstate"]) < f
                lower_bound = copy(mid_point)
            else
                upper_bound = copy(mid_point)
            end
            verbosity_mu>0 && @info "Iteration μ: $i => Lower bound: $lower_bound; Upper bound: $upper_bound"
            i+=1
        end
        if i>maxiter_mu
            @warn "The chemical potential lies between $lower_bound and $upper_bound, but did not converge within the tolerance. Increase maxiter_mu."
        else
            verbosity_mu>0 && @info "Final chemical potential = $mid_point"
        end

        if flag
            dictionary = dictionary_sp
        end

        dictionary["μ"] = mid_point
    end

    return dictionary
end
