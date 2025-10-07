abstract type Simulation end
name(s::Simulation) = string(typeof(s))

function Base.string(s::TensorKit.ProductSector{Tuple{FermionParity,SU2Irrep,U1Irrep}})
    parts = map(x -> sprint(show, x; context=:typeinfo => typeof(x)), s.sectors)
    return "[fℤ₂×SU₂×U₁]$(parts)"
end

function Base.string(s::TensorKit.ProductSector{Tuple{FermionParity,U1Irrep,U1Irrep}})
    parts = map(x -> sprint(show, x; context=:typeinfo => typeof(x)), s.sectors)
    return "[fℤ₂×U₁×U₁]$(parts)"
end

function Base.string(s::TensorKit.ProductSector{Tuple{FermionParity,SU2Irrep}})
    parts = map(x -> sprint(show, x; context=:typeinfo => typeof(x)), s.sectors)
    return "[fℤ₂×SU₂]$(parts)"
end

"""
    OB_Sim(t::Vector{Float64}, u::Vector{Float64}, μ=0.0, J::Vector{Float64}, P=1, Q=1, svalue=2.0, bond_dim=50, period=0; kwargs...)

Construct a parameter set for a 1D one-band Hubbard model with a fixed number of particles.

# Arguments
- `t`: Vector in which element ``n`` is the value of the hopping parameter of distance ``n``. The first element is the nearest-neighbour hopping.
- `u`: Vector in which element ``n`` is the value of the Coulomb interaction with site at distance ``n-1``. The first element is the on-site interaction.
- `J`: Vector in which element ``n`` is the value of the exchange interaction with site at distance ``n``. The first element is the nearest-neighbour exchange.
- `µ`: The chemical potential.
- `P`,`Q`: The ratio `P`/`Q` defines the number of electrons per site, which should be larger than 0 and smaller than 2.
- `svalue`: The Schmidt truncation value, used to truncate in the iDMRG2 algorithm for the computation of the groundstate.
- `bond_dim`: The maximal bond dimension used to initialize the state.
- `Period`: Perform simulations on a helix with circumference `Period`. Value 0 corresponds to an infinite chain.

Put the optional argument `spin=true` to perform spin-dependent calculations.
"""
struct OB_Sim <: Simulation
    t::Vector{Float64}
    u::Vector{Float64}
    μ::Float64
    J::Vector{Float64}
    P::Int64
    Q::Int64
    svalue::Float64
    bond_dim::Int64
    period::Int64
    kwargs
    function OB_Sim(t::Vector{Float64}, u::Vector{Float64}, μ::Float64=0.0, P::Int64=1, Q::Int64=1, svalue=2.0, bond_dim = 50, period = 0; kwargs...)
        return new(t, u, μ, [0.0], P, Q, svalue, bond_dim, period, kwargs)
    end
    function OB_Sim(t::Vector{Float64}, u::Vector{Float64}, μ::Float64=0.0, J::Vector{Float64}=[0.0], P::Int64=1, Q::Int64=1, svalue=2.0, bond_dim = 50, period = 0; kwargs...)
        return new(t, u, μ, J, P, Q, svalue, bond_dim, period, kwargs)
    end
end
name(::OB_Sim) = "OB"

"""
    MB_Sim(t::Matrix{Float64}, u::Matrix{Float64}, J::Matrix{Float64}, U13::Matrix{Float64}, P=1, Q=1, svalue=2.0, bond_dim=50; kwargs...)

Construct a parameter set for a 1D B-band Hubbard model with a fixed number of particles.

# Arguments
- `t`: Bx(nB) matrix in which element ``(i,j)`` is the hopping parameter from band ``i`` to band ``j``. The on-site, nearest neighbour, next-to-nearest neighbour... hopping matrices are concatenated horizontally.
- `u`: Bx(nB) matrix in which element ``(i,j)`` is the Coulomb repulsion ``U_{ij}=U_{iijj}`` between band ``i`` and band ``j``. The on-site, nearest neighbour, next-to-nearest neighbour... matrices are concatenated horizontally.
- `J`: Bx(nB) matrix in which element ``(i,j)`` is the exchange ``J_{ij}=U_{ijji}=U_{ijij}`` between band ``i`` and band ``j``. The on-site, nearest neighbour, next-to-nearest neighbour... matrices are concatenated horizontally. The diagonal terms of the on-site matrix are ignored.
- `U13`: BxB matrix in which element ``(i,j)`` is the parameter ``U_{ijjj}=U_{jijj}=U_{jjij}=U_{jjji}`` between band ``i`` and band ``j``. Only on-site. The diagonal terms of the on-site matrix are ignored. This argument is optional.
- `P`,`Q`: The ratio `P`/`Q` defines the number of electrons per site, which should be larger than 0 and smaller than 2.
- `svalue`: The Schmidt truncation value, used to truncate in the iDMRG2 algorithm for the computation of the groundstate.
- `bond_dim`: The maximal bond dimension used to initialize the state.

Put the optional argument 'spin=true' to perform spin-dependent calculations. 

U13 inter-site, Uijkk, and Uijkl can be inserted using kwargs.

Use the optional argument `name` to assign a name to the model. 
This is used to destinguish between different parameter sets: Wrong results could be loaded or overwritten if not used consistently!!!
"""
struct MB_Sim <: Simulation
    t::Matrix{Float64}                        #convention: number of bands = number of rows, BxB for on-site + Bx(B*range) matrix for IS
    u::Matrix{Float64}                        #convention: BxB matrix for OS (with OB on diagonal) + Bx(B*range) matrix for IS
    J::Matrix{Float64}                        #convention: BxB matrix for OS (with OB zeros) + Bx(B*range) matrix for IS
    U13::Matrix{Float64}                      #Matrix with iiij, iiji... parameters. Same convention.
    P::Int64
    Q::Int64
    svalue::Float64
    bond_dim::Int64
    kwargs
    function MB_Sim(t::Matrix{Float64}, u::Matrix{Float64}, J::Matrix{Float64}, P=1, Q=1, svalue=2.0, bond_dim = 50; kwargs...)
        Bands,_ = size(t)
        return new(t, u, J, zeros(Bands,Bands), P, Q, svalue, bond_dim, kwargs)
    end
    function MB_Sim(t::Matrix{Float64}, u::Matrix{Float64}, J::Matrix{Float64}, U13::Matrix{Float64}, P=1, Q=1, svalue=2.0, bond_dim = 50; kwargs...)
        return new(t, u, J, U13, P, Q, svalue, bond_dim, kwargs)
    end
end
name(::MB_Sim) = "MB"

"""
    OBC_Sim(t::Vector{Float64}, u::Vector{Float64}, μf::Float64, svalue=2.0, bond_dim=50, period=0; mu=true, kwargs...)

Construct a parameter set for a 1D one-band Hubbard model with the number of particles determined by a chemical potential.

# Arguments
- `t`: Vector in which element ``n`` is the value of the hopping parameter of distance ``n``. The first element is the nearest-neighbour hopping.
- `u`: Vector in which element ``n`` is the value of the Coulomb interaction with site at distance ``n-1``. The first element is the on-site interaction.
- `µf`: The chemical potential, if `mu=true`. Otherwise, the filling of the system. The chemical potential corresponding to the given filling is determined automatically.
- `svalue`: The Schmidt truncation value, used to truncate in the iDMRG2 algorithm for the computation of the groundstate.
- `bond_dim`: The maximal bond dimension used to initialize the state.
- `Period`: Perform simulations on a helix with circumference `Period`. Value 0 corresponds to an infinite chain.

Spin-dependent calculations are not yet implemented.
"""
struct OBC_Sim <: Simulation
    t::Vector{Float64}
    u::Vector{Float64}
    μ::Union{Float64, Nothing}    # Imposed chemical potential
    f::Union{Float64, Nothing}    # Fraction indicating the filling
    svalue::Float64
    bond_dim::Int64
    period::Int64
    kwargs
    function OBC_Sim(t, u, μf::Float64, svalue=2.0, bond_dim = 50, period = 0; mu=true, kwargs...)
        spin::Bool = get(kwargs, :spin, false)
        if spin
            error("Spin not implemented.")
        end
        if mu
            return new(t, u, μf, nothing, svalue, bond_dim, period, kwargs)
        else
            if 0 < μf < 2
                return new(t, u, nothing, μf, svalue, bond_dim, period, kwargs)
            else
                return error("Filling should be between 0 and 2.")
            end
        end
    end
end
name(::OBC_Sim) = "OBC"

# used to compute groundstates in µ iterations
struct OBC_Sim2 <: Simulation
    t::Vector{Float64}
    u::Vector{Float64}
    μ::Union{Float64, Nothing}    # Imposed chemical potential
    svalue::Float64
    bond_dim::Int64
    period::Int64
    kwargs
    function OBC_Sim2(t, u, μ::Float64, svalue=2.0, bond_dim = 50, period = 0; kwargs...)
        return new(t, u, μ, svalue, bond_dim, period, kwargs)
    end
end
name(::OBC_Sim2) = "OBC2"

"""
    MBC_Sim(t::Matrix{Float64}, u::Matrix{Float64}, J::Matrix{Float64}, U13::Matrix{Float64}, svalue=2.0, bond_dim=50; kwargs...)

Construct a parameter set for a 1D ``B``-band Hubbard model with the number of particles determined by a chemical potential.

# Arguments
- `t`: ``B\\times nB`` matrix in which element ``(i,j)`` is the hopping parameter from band ``i`` to band ``j``. The on-site, nearest neighbour, next-to-nearest neighbour... hopping matrices are concatenated horizontally. The diagonal terms of the on-site matrix determine the filling.
- `u`: ``B\\times nB`` matrix in which element ``(i,j)`` is the Coulomb repulsion ``U_{ij}=U_{iijj}`` between band ``i`` and band ``j``. The on-site, nearest neighbour, next-to-nearest neighbour... matrices are concatenated horizontally.
- `J`: ``B\\times nB`` matrix in which element ``(i,j)`` is the exchange ``J_{ij}=U_{ijji}=U_{ijij}`` between band ``i`` and band ``j``. The on-site, nearest neighbour, next-to-nearest neighbour... matrices are concatenated horizontally. The diagonal terms of the on-site matrix are ignored.
- `U13`: ``B\\times B`` matrix in which element ``(i,j)`` is the parameter ``U_{ijjj}=U_{jijj}=U_{jjij}=U_{jjji}`` between band ``i`` and band ``j``. Only on-site. The diagonal terms of the on-site matrix are ignored. This argument is optional.
- `svalue`: The Schmidt truncation value, used to truncate in the iDMRG2 algorithm for the computation of the groundstate.
- `bond_dim`: The maximal bond dimension used to initialize the state.

Spin-dependent calculations are not yet implemented.

U13 inter-site, Uijkk, and Uijkl can be inserted using kwargs.

Use the optional argument `name` to assign a name to the model. 
This is used to destinguish between different parameter sets: Wrong results could be loaded or overwritten if not used consistently!!!
"""
struct MBC_Sim <: Simulation
    t::Matrix{Float64}                        #convention: number of bands = number of rows, BxB for on-site + Bx(B*range) matrix for IS
    u::Matrix{Float64}                        #convention: BxB matrix for OS (with OB on diagonal) + Bx(B*range) matrix for IS
    J::Matrix{Float64}                        #convention: BxB matrix for OS (with OB zeros) + Bx(B*range) matrix for IS
    U13::Matrix{Float64}                      #Matrix with iiij, iiji... parameters. Same convention.
    svalue::Float64
    bond_dim::Int64
    kwargs
    function MBC_Sim(t::Matrix{Float64}, u::Matrix{Float64}, J::Matrix{Float64}, svalue=2.0, bond_dim = 50; kwargs...)
        spin::Bool = get(kwargs, :spin, false)
        if spin
            error("Spin not implemented.")
        end
        Bands,_ = size(t)
        return new(t, u, J, zeros(Bands,Bands), svalue, bond_dim, kwargs)
    end
    function MBC_Sim(t::Matrix{Float64}, u::Matrix{Float64}, J::Matrix{Float64}, U13::Matrix{Float64}, svalue=2.0, bond_dim = 50; kwargs...)
        spin::Bool = get(kwargs, :spin, false)
        if spin
            error("Spin not implemented.")
        end
        return new(t, u, J, U13, svalue, bond_dim, kwargs)
    end
end
name(::MBC_Sim) = "MBC"
