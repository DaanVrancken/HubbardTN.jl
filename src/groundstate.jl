###############
# Groundstate #
###############

function initialize_mps(operator, P::Int64, max_dimension::Int64, spin::Bool)
    Ps = physicalspace.(parent(operator))
    L = length(Ps)
    V_right = accumulate(fuse, Ps)
    
    V_l = accumulate(fuse, dual.(Ps); init=one(first(Ps)))
    V_left = reverse(V_l)
    len = length(V_left)
    step = length(V_left)-1
    V_left = [view(V_left,len-step+1:len); view(V_left,1:len-step)]   # same as circshift(V_left,1)

    V = TensorKit.infimum.(V_left, V_right)

    if !spin
        Vmax = Vect[(FermionParity ⊠ Irrep[SU₂] ⊠ Irrep[U₁])]((0,0,0)=>1)     # find maximal virtual space
        for i in 0:1
            for j in 0:1//2:3
                for k in -(L*P):1:(L*P)
                    Vmax = Vect[(FermionParity ⊠ Irrep[SU₂] ⊠ Irrep[U₁])]((i,j,k)=>max_dimension) ⊕ Vmax
                end
            end
        end
    else
        Vmax = Vect[(FermionParity ⊠ Irrep[U₁] ⊠ Irrep[U₁])]((0,0,0)=>1)
        for i in 0:1
            for j in -L:1:L
                for k in -(L*P):1:(L*P)
                    Vmax = Vect[(FermionParity ⊠ Irrep[U₁] ⊠ Irrep[U₁])]((i,j,k)=>max_dimension) ⊕ Vmax
                end
            end
        end
    end

    V_max = copy(V)

    for i in 1:length(V_right)
        V_max[i] = Vmax
    end

    V_trunc = TensorKit.infimum.(V,V_max)

    return InfiniteMPS(Ps, V_trunc)
end

function initialize_mps(operator, max_dimension::Int64)
    Ps = physicalspace.(parent(operator))

    V_right = accumulate(fuse, Ps)
    
    V_l = accumulate(fuse, dual.(Ps); init=one(first(Ps)))
    V_left = reverse(V_l)
    len = length(V_left)
    step = length(V_left)-1
    V_left = [view(V_left,len-step+1:len); view(V_left,1:len-step)]   # same as circshift(V_left,1)

    V = TensorKit.infimum.(V_left, V_right)

    Vmax = Vect[(FermionParity ⊠ Irrep[SU₂])]((0,0)=>1)     # find maximal virtual space

    for i in 0:1
        for j in 0:1//2:3
            Vmax = Vect[(FermionParity ⊠ Irrep[SU₂])]((i,j)=>max_dimension) ⊕ Vmax
        end
    end

    V_max = copy(V)      # if no copy(), V will change along when V_max is changed

    for i in 1:length(V_right)
        V_max[i] = Vmax
    end

    V_trunc = TensorKit.infimum.(V,V_max)

    return InfiniteMPS(Ps, V_trunc)
end

function compute_groundstate(simul::Union{OB_Sim, MB_Sim, OBC_Sim2, MBC_Sim}; tol::Float64=1e-6, verbosity::Int64=0, maxiter::Int64=1000)
    H = hamiltonian(simul)
    spin::Bool = get(simul.kwargs, :spin, false)
    init_state = get(simul.kwargs, :init_state, nothing)

    if isnothing(init_state)
        if hasproperty(simul, :P)
            ψ₀ = initialize_mps(H,simul.P,simul.bond_dim,spin)
        else
            ψ₀ = initialize_mps(H,simul.bond_dim)
        end
    else
        ψ₀ = init_state
    end
    
    schmidtcut = 10.0^(-simul.svalue)
    
    if length(H) > 1
        ψ₀, envs, = find_groundstate(ψ₀, H, IDMRG2(; trscheme=truncbelow(schmidtcut), tol=tol, verbosity=verbosity))
    else
        ψ₀, envs, = find_groundstate(ψ₀, H, VUMPS(; tol=max(tol, schmidtcut/10), verbosity=verbosity))
        ψ₀ = changebonds(ψ₀, SvdCut(; trscheme=truncbelow(schmidtcut)))
        χ = sum(i -> dim(left_virtualspace(ψ₀, i)), 1:length(H))
        for i in 1:maxiter
            ψ₀, envs = changebonds(ψ₀, H, VUMPSSvdCut(;trscheme=truncbelow(schmidtcut)))
            ψ₀, = find_groundstate(ψ₀, H, VUMPS(; tol=max(tol, schmidtcut / 10), verbosity=verbosity), envs)
            ψ₀ = changebonds(ψ₀, SvdCut(; trscheme=truncbelow(schmidtcut)))
            χ′ = sum(i -> dim(left_virtualspace(ψ₀, i)), 1:length(H))
            isapprox(χ, χ′; rtol=0.05) && break
            χ = χ′
        end
    end
    
    alg = VUMPS(; maxiter=maxiter, tol=tol, verbosity=verbosity) &
        GradientGrassmann(; maxiter=maxiter, tol=tol, verbosity=verbosity)
    ψ, envs, δ = find_groundstate(ψ₀, H, alg)
    
    return Dict("groundstate" => ψ, "environments" => envs, "ham" => H, "delta" => δ, "config" => simul)
end

function compute_groundstate(simul::OBC_Sim; tol::Float64=1e-6, verbosity::Int64=0, maxiter::Int64=1000)
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

"""
    produce_groundstate(model::Simulation; force::Bool=false, path::String="")

Compute or load groundstate of the `model`. It can be stored at or loaded from `path`. If `force=true`, overwrite existing calculation.
"""
function produce_groundstate(simul::Union{MB_Sim, MBC_Sim}; force::Bool=false, path::String="")
    code = get(simul.kwargs, :code, "")
    S = "nospin_"
    spin::Bool = get(simul.kwargs, :spin, false)
    if spin
        S = "spin_"
    end

    data, _ = produce_or_load(compute_groundstate, simul, path; prefix="groundstate_"*S*code, force=force)
    return data
end

function produce_groundstate(simul::Union{OB_Sim, OBC_Sim}; force::Bool=false, path::String="")
    t = simul.t 
    u = simul.u
    if hasproperty(simul, :J)
        J = simul.J
    else
        J = 0
    end
    S_spin = "nospin_"
    spin::Bool = get(simul.kwargs, :spin, false)
    if spin
        S_spin = "spin_"
    end
    U13::Vector{Float64} = get(simul.kwargs, :U13, [0.0])
    JMs::Tuple{Float64, Float64} = get(simul.kwargs, :JMs, (0.0,0.0))
    J_inter = JMs[1]
    Ms = JMs[2]
    S = "groundstate_"*S_spin*"t$(t)_u$(u)_J$(J)_U13$(U13)_JMs$(J_inter)_$(Ms)"
    S = replace(S, ", " => "_")
    data, _ = produce_or_load(compute_groundstate, simul, path; prefix=S, force=force)
    return data
end