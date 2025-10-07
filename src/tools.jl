##############
# Truncation #
##############

function TruncState(simul::Simulation, trunc_dim::Int64; trunc_scheme::Int64=0, path_gs=path_gs)
    if trunc_dim<=0
        return error("trunc_dim should be a positive integer.")
    end
    if trunc_scheme!=0 && trunc_scheme!=1
        return error("trunc_scheme should be either 0 (VUMPSSvdCut) or 1 (SvdCut).")
    end

    dictionary = produce_groundstate(simul; path=path_gs)
    ψ = dictionary["groundstate"]
    H = dictionary["ham"]
    if trunc_scheme==0
        ψ, envs = changebonds(ψ,H,VUMPSSvdCut(; trscheme=truncdim(trunc_dim)))
    else
        ψ, envs = changebonds(ψ,H,SvdCut(; trscheme=truncdim(trunc_dim)))
    end
    return  Dict("ψ_trunc" => ψ, "envs_trunc" => envs)
end

"""
    produce_truncstate(model::Simulation, trunc_dim::Int64; trunc_scheme::Int64=0, force::Bool=false, path::String="", path_gs::String=path)

Compute or load a truncated approximation of the groundstate.

# Arguments
- `model`: Model for which the groundstate is to be truncated.
- `trunc_dim`: Maximal bond dimension of the truncated state.
- `trunc_scheme`: Scheme to perform the truncation. 0 = VUMPSSvdCut. 1 = SvdCut.
- `force`: If true, overwrite existing calculation.
- `path`: Path to save/load the calculation.
- `path_gs`: Path to load the groundstate from, if different from `path`.
"""
function produce_TruncState(simul::Simulation, trunc_dim::Int64; trunc_scheme::Int64=0, force::Bool=false, path::String="", path_gs::String=path)
    code = get(simul.kwargs, :code, "")
    data, _ = produce_or_load(simul, path; prefix="Trunc_GS_"*code*"_dim=$trunc_dim"*"_scheme=$trunc_scheme", force=force) do cfg
        return TruncState(cfg, trunc_dim; trunc_scheme=trunc_scheme, path_gs=path_gs)
    end
    return data
end


####################
# State properties #
####################

"""
    dim_state(ψ::InfiniteMPS)

Determine the bond dimensions in an infinite MPS.
"""
function dim_state(ψ::InfiniteMPS)
    dimension = Int64.(zeros(length(ψ)))
    for i in 1:length(ψ)
        dimension[i] = dim(space(ψ.AL[i],1))
    end
    return dimension
end

"""
    density_spin(model::Union{OB_Sim,MB_Sim}; path::String="")

Compute the density of spin up and spin down per site in the unit cell for the ground state stored at `path`.
"""
function density_spin(simul::Union{OB_Sim,MB_Sim}; path::String="")
    P = simul.P;
    Q = simul.Q

    dictionary = produce_groundstate(simul; path=path);
    ψ₀ = dictionary["groundstate"];
    
    spin::Bool = get(simul.kwargs, :spin, false)

    if !spin
        error("This system is spin independent.")
    end

    return density_spin(ψ₀, P, Q)
end

function density_spin(ψ₀::InfiniteMPS, P::Int64, Q::Int64)
    I, Ps = SymSpace(P,Q,true)
    if iseven(P)
        T = Q
    else 
        T = 2*Q
    end
    Bands = Int(length(ψ₀)/T)

    nup = zeros(ComplexF64, Ps ← Ps)
    blocks(nup)[I((0, 0, 2*Q-P))] .= 1
    blocks(nup)[I((1, 1, Q-P))] .= 1
    ndown = zeros(ComplexF64, Ps ← Ps)
    blocks(ndown)[I((0, 0, 2*Q-P))] .= 1
    blocks(ndown)[I((1, -1, Q-P))] .= 1

    Nup = zeros(Bands,T);
    Ndown = zeros(Bands,T);
    for i in 1:Bands
        for j in 1:T
            Nup[i,j] = real(expectation_value(ψ₀, (i+(j-1)*Bands) => nup))
            Ndown[i,j] = real(expectation_value(ψ₀, (i+(j-1)*Bands) => ndown))
        end
    end

    return Nup, Ndown
end

"""
    calc_ms(model::Union{OB_Sim,MB_Sim})

Compute the staggered magnetization of the ground state.
"""
function calc_ms(model::Union{OB_Sim,MB_Sim})
    up, down = hf.density_spin(model)
    Mag = up - down
    if !all(x -> isapprox(abs(x),abs(Mag[1,1]),rtol=10^(-6)), vec(Mag))
        @warn "Spin-density wave?"
    end
    return abs(Mag[1,1])
end

"""
    density_state(model::Simulation; path::String="")

Compute the number of electrons per site in the unit cell for the ground state stored at `path`.
"""
function density_state(simul::Union{OB_Sim,MB_Sim}; path::String="")
    P = simul.P;
    Q = simul.Q

    dictionary = produce_groundstate(simul; path=path);
    ψ₀ = dictionary["groundstate"];
    
    spin::Bool = get(simul.kwargs, :spin, false)

    return density_state(ψ₀, P, Q, spin)
end

function density_state(simul::Union{OBC_Sim, MBC_Sim}; path::String="")
    dictionary = produce_groundstate(simul; path=path);
    ψ = dictionary["groundstate"];

    return density_state(ψ)
end

# For Hubbard models without chemical potential
function density_state(ψ₀::InfiniteMPS,P::Int64,Q::Int64,spin::Bool)
    if iseven(P)
        T = Q
    else 
        T = 2*Q
    end
    Bands = Int(length(ψ₀)/T)

    n = Number(P,Q,spin)

    Nₑ = zeros(Bands*T,1);
    for i in 1:(Bands*T)
        Nₑ[i] = real(expectation_value(ψ₀, i => n))
    end
    
    N_av = zeros(Bands,1)
    for i in 1:Bands
        av = 0
        for j in 0:(T-1)
            av = Nₑ[i+Bands*j] + av
        end
        N_av[i,1] = av/T
    end

    check = (sum(Nₑ)/(T*Bands) ≈ P/Q)
    println("Filling is conserved: $check")

    return Nₑ
end

# For Hubbard models involving a chemical potential
function density_state(ψ::InfiniteMPS)
    Bands = length(ψ)

    n = Number()

    Nₑ = zeros(Bands);
    for i in 1:Bands
        Nₑ[i] = real(expectation_value(ψ, i => n))
    end

    if Bands==1
        # convert 1x1 matrix into scalar
        Nₑ = sum(Nₑ)
    end

    return Nₑ
end


####################
# Tools & Plotting #
####################

"""
    plot_excitations(momenta, energies; title="Excitation_energies", l_margin=[15mm 0mm])

Plot the obtained energy levels in functions of the momentum.
"""
function plot_excitations(momenta, Es; title="Excitation energies", l_margin=[15mm 0mm])
    _, nums = size(Es)
    plot(momenta,real(Es[:,1]), label="", linecolor=:blue, title=title, left_margin=l_margin)
    for i in 2:nums
        plot!(momenta,real(Es[:,i]), label="", linecolor=:blue)
    end
    xlabel!("k")
    ylabel!("Energy density")
end

"""
    plot_spin(model::Simulation; title="Spin Density", l_margin=[15mm 0mm])

Plot the spin density of the model throughout the unit cell as a heatmap.
"""
function plot_spin(model::Simulation; title="Spin Density", l_margin=[15mm 0mm])
    up, down = hf.density_spin(model)
    Sz = up - down
    heatmap(Sz, color=:grays, c=:grays, label="", xlabel="Site", ylabel="Band", title=title, clims=(-1, 1))
end

"""
    extract_params(path::String; range_u::Int64= 1, range_t::Int64=2, range_J::Int64=1, 
                        range_U13::Int64=1, r_1111::Int64 = 1, r_112::Int64 = 1)

Extract the parameters from a params.jl file located at `path` in PyFoldHub format.
"""
function extract_params(path::String; range_u::Int64= 1, range_t::Int64=2, range_J::Int64=1, 
                        range_U13::Int64=1, r_1111::Int64 = 1, r_112::Int64 = 1)
    # Wmn should be rank 8 tensor (only one frequency point)
    include(path)

    B = size(Wmn)[5]
    site_0 = ceil(Int,size(Wmn)[1]/2)

    t = zeros(B,B*range_t)
    U = zeros(B,B*range_u)
    J = zeros(B,B*range_J)
    U13_OS = zeros(B,B)
    if range_U13 == 1
        U13_IS = zeros(B,B*range_U13,4)
    else
        U13_IS = zeros(B,B*(range_U13-1),4)
    end
    for i in 1:B
        for j in 1:B
            for r in 0:(range_t-1)
                t[i,j+r*B] = tmn[site_0+r,i,j] + corr_H[site_0+r,i,j] #+ corr_G_HW[site_0+r,i,j] + corr_v_xc[site_0+r,i,j], check minus sign...
            end
            for r in 0:(range_u-1)
                U[i,j+r*B] = Wmn[site_0,site_0,site_0+r,site_0+r,i,i,j,j]
            end
            for r in 0:(range_J-1)
                if r!=0 || i!=j
                    J[i,j+r*B] = Wmn[site_0,site_0+r,site_0+r,site_0,i,j,j,i]
                    if !(J[i,j+r*B] ≈ Wmn[site_0,site_0+r,site_0,site_0+r,i,j,i,j])
                        error("J1 is not equal to J2 at (r,i,j)=($r,$i,$j).")
                    end
                end
            end
            for r in 1:(range_U13-1)
                U13_IS[i,j+(r-1)*B,1] = Wmn[site_0,site_0+r,site_0+r,site_0+r,i,j,j,j]
                U13_IS[i,j+(r-1)*B,2] = Wmn[site_0+r,site_0+r,site_0,site_0+r,j,j,i,j]
                U13_IS[i,j+(r-1)*B,3] = Wmn[site_0+r,site_0,site_0,site_0,j,i,i,i]
                U13_IS[i,j+(r-1)*B,4] = Wmn[site_0,site_0,site_0+r,site_0,i,i,j,i]
                if !(U13_IS[i,j+(r-1)*B,1] ≈ Wmn[site_0+r,site_0,site_0+r,site_0+r,j,i,j,j]) || !(U13_IS[i,j+(r-1)*B,2] ≈ Wmn[site_0+r,site_0+r,site_0+r,site_0,j,j,j,i]) ||
                    !(U13_IS[i,j+(r-1)*B,3] ≈ Wmn[site_0,site_0+r,site_0,site_0,i,j,i,i]) || !(U13_IS[i,j+(r-1)*B,4] ≈ Wmn[site_0,site_0,site_0,site_0+r,i,i,i,j])
                    error("U13_IS not consistent.")
                end
            end
            if i != j
                U13_OS[i,j] = Wmn[site_0,site_0,site_0,site_0,i,j,j,j]
                if !isapprox(U13_OS[i,j], Wmn[site_0,site_0,site_0,site_0,j,i,j,j], rtol=1e-3) || !isapprox(U13_OS[i,j], Wmn[site_0,site_0,site_0,site_0,j,j,i,j], rtol=1e-3) || 
                    !isapprox(U13_OS[i,j], Wmn[site_0,site_0,site_0,site_0,j,j,j,i], rtol=1e-3)
                    @warn "U13_OS not consistent at i=$i, j=$j, for rtol=1e-3."
                    if !isapprox(U13_OS[i,j], Wmn[site_0,site_0,site_0,site_0,j,i,j,j], atol=1e-3) || !isapprox(U13_OS[i,j], Wmn[site_0,site_0,site_0,site_0,j,j,i,j], atol=1e-3) || 
                        !isapprox(U13_OS[i,j], Wmn[site_0,site_0,site_0,site_0,j,j,j,i], atol=1e-3)
                        error("U13_OS not consistent at i=$i, j=$j.")
                    end
                end
            end
        end
    end

    #shift chemical potential
    mu = minimum(diag(t[:,1:B]))
    t[:,1:B] -= mu.*I

    U112 = Dict{Tuple{Int, Int, Int, Int}, Float64}()
    for i in 1:r_112*B, j in 1:r_112*B, k in 1:r_112*B, l in 1:r_112*B
        if (i == j || i == k || i == l || j == k || j == l || k == l) && length(unique((i, j, k, l))) == 3 && minimum((i,j,k,l)) <= B
            mod_i = mod(i-1,B) + 1; r_i = (i - 1) ÷ B;
            mod_j = mod(j-1,B) + 1; r_j = (j - 1) ÷ B;
            mod_k = mod(k-1,B) + 1; r_k = (k - 1) ÷ B;
            mod_l = mod(l-1,B) + 1; r_l = (l - 1) ÷ B;
            # change index order to those of operators
            U112[(i,k,l,j)] = Wmn[site_0+r_i,site_0+r_j,site_0+r_k,site_0+r_l,mod_i,mod_j,mod_k,mod_l]
        end
    end

    U1111 = Dict{Tuple{Int, Int, Int, Int}, Float64}()
    for i in 1:r_1111*B, j in 1:r_1111*B, k in 1:r_1111*B, l in 1:r_1111*B
        if length(unique((i, j, k, l))) == 4 && minimum((i,j,k,l)) <= B
            mod_i = mod(i-1,B) + 1; r_i = (i - 1) ÷ B;
            mod_j = mod(j-1,B) + 1; r_j = (j - 1) ÷ B;
            mod_k = mod(k-1,B) + 1; r_k = (k - 1) ÷ B;
            mod_l = mod(l-1,B) + 1; r_l = (l - 1) ÷ B;
            # change index order to those of operators
            U1111[(i,k,l,j)] = Wmn[site_0+r_i,site_0+r_j,site_0+r_k,site_0+r_l,mod_i,mod_j,mod_k,mod_l]
        end
    end

    return t, U, J, U13_OS, U13_IS, U112, U1111
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
a message after each tensor is successfully saved.
"""
function save_state(ψ::InfiniteMPS, path::String, name::String)
    path = joinpath(path,name)
    mkdir(path)
    for i in 1:length(ψ)
        d = convert(Dict,ψ.AL[i])
        @save joinpath(path,"state$i.jld2") d
        println("State $i saved.")
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