###############
# Excitations #
###############

function compute_excitations(simul::Simulation, momenta, nums::Int64; 
                                    charges::Vector{Float64}=[0,0.0,0], path_gs::String="",
                                    trunc_dim::Int64=0, trunc_scheme::Int64=0, DW = false, shift=1,
                                    solver=Arnoldi(;krylovdim=30,tol=1e-6,eager=true))
    if trunc_dim<0
        return error("Trunc_dim should be a positive integer.")
    end
    spin::Bool = get(simul.kwargs, :spin, false)

    if hasproperty(simul, :Q)
        Q = simul.Q
        if !spin
            sector = fℤ₂(charges[1]) ⊠ SU2Irrep(charges[2]) ⊠ U1Irrep(charges[3]*Q)
        else
            sector = fℤ₂(charges[1]) ⊠ U1Irrep(charges[2]) ⊠ U1Irrep(charges[3]*Q)
        end
    else
        sector = fℤ₂(charges[1]) ⊠ SU2Irrep(charges[2])
    end

    dictionary = produce_groundstate(simul; path=path_gs)
    ψ = dictionary["groundstate"]
    H = dictionary["ham"]
    if trunc_dim==0
        envs = dictionary["environments"]
    else
        dict_trunc = produce_TruncState(simul, trunc_dim; path=path_gs, trunc_scheme=trunc_scheme)
        ψ = dict_trunc["ψ_trunc"]
        envs = dict_trunc["envs_trunc"]
    end
    if DW
        ψ_s = circshift(ψ, shift)
        envs_s = environments(ψ_s, H);
        Es, qps = excitations(H, QuasiparticleAnsatz(solver, MPSKit.Defaults.alg_environments(;dynamic_tols=false)), momenta./length(H), ψ, envs, ψ_s, envs_s; num=nums, sector=sector)
    else
        Es, qps = excitations(H, QuasiparticleAnsatz(solver, MPSKit.Defaults.alg_environments(;dynamic_tols=false)), momenta./length(H), ψ, envs; num=nums, sector=sector)
    end
    
    return Dict("Es" => Es, "qps" => qps, "momenta" => momenta)
end

"""
    produce_excitations(model::Simulation, momenta, nums::Int64; path::String="", path_gs::String=path, force::Bool=false, charges::Vector{Float64}=[0,0.0,0], kwargs...)

Compute or load quasiparticle excitations of the desired `model`.

# Arguments
- `model`: Model for which excitations are sought.
- `momenta`: Momenta of the quasiparticle excitations.
- `nums`: Number of excitations.
- `path`: Path to save/load the calculation.
- `path_gs`: Path to load the groundstate from, if different from `path`.
- `force`: If true, overwrite existing calculation.
- `charges`: Charges of the symmetry sector of the excitations.
"""
function produce_excitations(simul::Simulation, momenta, nums::Int64; path::String="", path_gs::String=path,
                                    force::Bool=false, charges::Vector{Float64}=[0,0.0,0], 
                                    trunc_dim::Int64=0, trunc_scheme::Int64=0, 
                                    solver=Arnoldi(;krylovdim=30,tol=1e-6,eager=true))
    spin::Bool = get(simul.kwargs, :spin, false)
    S = ""
    param_string = ""
    if size(simul.u)[1] == 1
        if hasproperty(simul, :J)
            J = simul.J
        else
            J = 0
        end
        U13::Vector{Float64} = get(simul.kwargs, :U13, [0.0])
        JMs::Tuple{Float64, Float64} = get(simul.kwargs, :JMs, (0.0,0.0))
        J_inter = JMs[1]
        Ms = JMs[2]
        param_string = "t$(simul.t)u$(simul.u)J$(J)U$(U13)m$(J_inter)_$(Ms)_"
    end
    if typeof(momenta)==Float64
        momenta_string = "_k=$momenta"
    else
        momenta_string = "_k=$(first(momenta))to$(last(momenta))div$(length(momenta))"
    end
    if hasproperty(simul, :Q)
        if !spin
            charge_string = "f$(Int(charges[1]))su$(charges[2])u$(Int(charges[3]))"
        else
            charge_string = "f$(Int(charges[1]))u$(charges[2])u$(Int(charges[3]))"
            S = "spin_"
        end
    else
        charge_string = "f$(Int(charges[1]))su$(charges[2])"
    end
    code = get(simul.kwargs, :code, "")
    Prefix = "exc_"*S*param_string*code*"_N=$nums"*"c="*charge_string*momenta_string*"_tr=$trunc_dim"
    Prefix = replace(Prefix, "__" => "_")
    Prefix = replace(Prefix, "3.141592653589793" => "pi")

    data, _ = produce_or_load(simul, path; prefix=Prefix, force=force) do cfg
        return compute_excitations(cfg, momenta, nums; path_gs=path_gs, charges=charges, trunc_dim=trunc_dim, trunc_scheme=trunc_scheme, solver=solver)
    end
    return data
end

"""
    produce_bandgap(model::Union{OB_Sim, MB_Sim}; resolution::Int64=5, force::Bool=false)

Compute or load the band gap of the desired model.
"""
function produce_bandgap(simul::Union{OB_Sim, MB_Sim}; path::String="", path_gs::String=path, resolution::Int64=5, force::Bool=false)
    momenta = range(0, π, resolution)
    spin::Bool = get(simul.kwargs, :spin, false)

    if spin
        error("Band gap for spin systems not implemented.")
    end

    Exc_hole = produce_excitations(simul, momenta, 1; path=path, path_gs=path_gs, force=force, charges=[1,1/2,-1])
    Exc_elec = produce_excitations(simul, momenta, 1; path=path, path_gs=path_gs, force=force, charges=[1,1/2,1])

    E_hole = real(Exc_hole["Es"])
    E_elec = real(Exc_elec["Es"])

    E_tot = E_hole + E_elec

    gap, k = findmin(E_tot[:,1])

    if k != 1
        @warn "Indirect band gap! Higher resolution might be required."
    end

    return gap, momenta[k]
end

function produce_domainwalls(simul::Simulation, momenta, nums::Int64; path::String="", path_gs::String=path,
                                    force::Bool=false, charges::Vector{Float64}=[0,0.0,1], 
                                    trunc_dim::Int64=0, trunc_scheme::Int64=0, shift=1,
                                    solver=Arnoldi(;krylovdim=30,tol=1e-6,eager=true))
    spin::Bool = get(simul.kwargs, :spin, false)
    S = ""
    param_string = ""
    if size(simul.u)[1] == 1
        if hasproperty(simul, :J)
            J = simul.J
        else
            J = 0
        end
        U13::Vector{Float64} = get(simul.kwargs, :U13, [0.0])
        JMs::Tuple{Float64, Float64} = get(simul.kwargs, :JMs, (0.0,0.0))
        J_inter = JMs[1]
        Ms = JMs[2]
        param_string = "t$(simul.t)u$(simul.u)J$(J)U$(U13)m$(J_inter)_$(Ms)_"
    end
    if typeof(momenta)==Float64
        momenta_string = "_k=$momenta"
    else
        momenta_string = "_k=$(first(momenta))to$(last(momenta))div$(length(momenta))"
    end
    if hasproperty(simul, :Q)
        if !spin
            charge_string = "f$(Int(charges[1]))su$(charges[2])u$(Int(charges[3]))"
        else
            charge_string = "f$(Int(charges[1]))u$(charges[2])u$(Int(charges[3]))"
            S = "spin_"
        end
    else
        charge_string = "f$(Int(charges[1]))su$(charges[2])"
    end
    code = get(simul.kwargs, :code, "")
    Prefix = "exc_"*S*param_string*code*"_N=$nums"*"c="*charge_string*momenta_string*"_tr=$trunc_dim"
    Prefix = replace(Prefix, "__" => "_")
    Prefix = replace(Prefix, "3.141592653589793" => "pi")

    data, _ = produce_or_load(simul, path; prefix=Prefix, force=force) do cfg
        return compute_excitations(cfg, momenta, nums; path_gs=path_gs, charges=charges, trunc_dim=trunc_dim, trunc_scheme=trunc_scheme, DW=true, shift=shift, solver=solver)
    end
    return data
end