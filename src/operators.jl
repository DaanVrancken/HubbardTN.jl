##########
# Spaces #
##########

"""
    hubbard_space(particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; (P,Q)::Tuple{Int64, Int64}=(1,1))

Return the local hilbert space for a Hubbard-type model with the given particle and spin symmetries.
The possible symmetries are `Trivial`, `U1Irrep`, and `SU2Irrep`, for both particle number and spin.
When using `U1Irrep` particle symmetry, the filling can be adjusted to `P` particles per `Q` sites by
passing the keyword `filling=(P,Q)`. The default is `(1,1)`, i.e., half-filling.
"""
function hubbard_space(::Type{Trivial} = Trivial, ::Type{Trivial} = Trivial; kwargs...)
    return Vect[FermionParity](0 => 2, 1 => 2)
end
function hubbard_space(::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    return Vect[FermionParity ⊠ U1Irrep]((0, 0) => 2, (1, 1 // 2) => 1, (1, -1 // 2) => 1)
end
function hubbard_space(::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    return Vect[FermionParity ⊠ SU2Irrep]((0, 0) => 2, (1, 1 // 2) => 1)
end
function hubbard_space(::Type{U1Irrep}, ::Type{Trivial}; filling::Tuple{Int64,Int64}=(1,1))
    P, Q = filling
    return Vect[FermionParity ⊠ U1Irrep]((0, -P) => 1, (1, Q-P) => 2, (0, 2Q-P) => 1)
end
function hubbard_space(::Type{U1Irrep}, ::Type{U1Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    P, Q = filling
    return Vect[FermionParity ⊠ U1Irrep ⊠ U1Irrep](
        (0, -P, 0) => 1, (1, Q-P, 1 // 2) => 1,
        (1, Q-P, -1 // 2) => 1, (0, 2Q-P, 0) => 1
    )
end
function hubbard_space(::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    P, Q = filling
    return Vect[FermionParity ⊠ U1Irrep ⊠ SU2Irrep](
        (0, -P, 0) => 1, (1, Q-P, 1 // 2) => 1, (0, 2Q-P, 0) => 1
    )
end
# How to pass (P,Q) to the SU2 case?
function hubbard_space(::Type{SU2Irrep}, ::Type{Trivial}; kwargs...)
    return Vect[FermionParity ⊠ SU2Irrep]((0, 0) => 2, (1, 1 // 2) => 1)
end
function hubbard_space(::Type{SU2Irrep}, ::Type{U1Irrep}; kwargs...)
    return Vect[FermionParity ⊠ SU2Irrep ⊠ U1Irrep]((0, 0, 0) => 1, (1, 1 // 2, 1) => 1)
end
function hubbard_space(::Type{SU2Irrep}, ::Type{SU2Irrep}; kwargs...)
    return Vect[FermionParity ⊠ SU2Irrep ⊠ SU2Irrep]((1, 1 // 2, 1 // 2) => 1)
end


#############
# Operators #
#############

function single_site_operator(
        T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; filling::Tuple{Int64,Int64}=(1,1)
    )
    V = hubbard_space(particle_symmetry, spin_symmetry; filling=filling)
    return zeros(T, V ← V)
end

function two_site_operator(
        T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; filling::Tuple{Int64,Int64}=(1,1)
    )
    V = hubbard_space(particle_symmetry, spin_symmetry; filling=filling)
    return zeros(T, V ⊗ V ← V ⊗ V)
end

"""
    c_plusmin_up(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the two-body operator ``c†_{1,↑}, c_{2,↑}`` that creates a spin-up electron at the first site and annihilates a spin-up electron at the second.
"""
c_plusmin_up(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_plusmin_up(ComplexF64, P, S; kwargs...)
function c_plusmin_up(T, ::Type{Trivial}, ::Type{Trivial}; kwargs...)
    t = two_site_operator(T, Trivial, Trivial)
    I = sectortype(t)
    t[(I(1), I(0), dual(I(0)), dual(I(1)))][1, 1, 1, 1] = 1
    t[(I(1), I(1), dual(I(0)), dual(I(0)))][1, 2, 1, 2] = 1
    t[(I(0), I(0), dual(I(1)), dual(I(1)))][2, 1, 2, 1] = -1
    t[(I(0), I(1), dual(I(1)), dual(I(0)))][2, 2, 2, 2] = -1
    return t
end
function c_plusmin_up(T, ::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    t = two_site_operator(T, Trivial, U1Irrep)
    I = sectortype(t)
    t[(I(1, 1 // 2), I(0, 0), dual(I(0, 0)), dual(I(1, 1 // 2)))][1, 1, 1, 1] = 1
    t[(I(1, 1 // 2), I(1, -1 // 2), dual(I(0, 0)), dual(I(0, 0)))][1, 1, 1, 2] = 1
    t[(I(0, 0), I(0, 0), dual(I(1, -1 // 2)), dual(I(1, 1 // 2)))][2, 1, 1, 1] = -1
    t[(I(0, 0), I(1, -1 // 2), dual(I(1, -1 // 2)), dual(I(0, 0)))][2, 1, 1, 2] = -1
    return t
end
function c_plusmin_up(T, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    return error("Not implemented")
end
function c_plusmin_up(T, ::Type{U1Irrep}, ::Type{Trivial}; filling::Tuple{Int64,Int64}=(1,1))
    t = two_site_operator(T, U1Irrep, Trivial; filling=filling)
    P, Q = filling
    I = sectortype(t)
    t[(I(1, Q-P), I(0, -P), dual(I(0, -P)), dual(I(1, Q-P)))][1, 1, 1, 1] = 1
    t[(I(1, Q-P), I(1, Q-P), dual(I(0, -P)), dual(I(0, 2Q-P)))][1, 2, 1, 1] = 1
    t[(I(0, 2Q-P), I(0, -P), dual(I(1, Q-P)), dual(I(1, Q-P)))][1, 1, 2, 1] = -1
    t[(I(0, 2Q-P), I(1, Q-P), dual(I(1, Q-P)), dual(I(0, 2Q-P)))][1, 2, 2, 1] = -1
    return t
end
function c_plusmin_up(T, ::Type{U1Irrep}, ::Type{U1Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    t = two_site_operator(T, U1Irrep, U1Irrep; filling=filling)
    P, Q = filling
    I = sectortype(t)
    t[(I(1, Q-P, 1 // 2), I(0, -P, 0), dual(I(0, -P, 0)), dual(I(1, Q-P, 1 // 2)))] .= 1
    t[(I(1, Q-P, 1 // 2), I(1, Q-P, -1 // 2), dual(I(0, -P, 0)), dual(I(0, 2Q-P, 0)))] .= 1
    t[(I(0, 2Q-P, 0), I(0, -P, 0), dual(I(1, Q-P, -1 // 2)), dual(I(1, Q-P, 1 // 2)))] .= -1
    t[(I(0, 2Q-P, 0), I(1, Q-P, -1 // 2), dual(I(1, Q-P, -1 // 2)), dual(I(0, 2Q-P, 0)))] .= -1
    return t
end
function c_plusmin_up(T, ::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    return error("Not implemented")
end
function c_plusmin_up(T, ::Type{SU2Irrep}, ::Type{Trivial}; kwargs...)
    return error("Not implemented")
end
function c_plusmin_up(T, ::Type{SU2Irrep}, ::Type{U1Irrep}; kwargs...)
    return error("Not implemented")
end
function c_plusmin_up(T, ::Type{SU2Irrep}, ::Type{SU2Irrep}; kwargs...)
    return error("Not implemented")
end

"""
    c_plusmin_down(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the two-body operator ``c†_{1,↓}, c_{2,↓}`` that creates a spin-down electron at the first site and annihilates a spin-down electron at the second.
"""
c_plusmin_down(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_plusmin_down(ComplexF64, P, S; kwars...)
function c_plusmin_down(T, ::Type{Trivial}, ::Type{Trivial}; kwargs...)
    t = two_site_operator(T, Trivial, Trivial)
    I = sectortype(t)
    t[(I(1), I(0), dual(I(0)), dual(I(1)))][2, 1, 1, 2] = 1
    t[(I(1), I(1), dual(I(0)), dual(I(0)))][2, 1, 1, 2] = -1
    t[(I(0), I(0), dual(I(1)), dual(I(1)))][2, 1, 1, 2] = 1
    t[(I(0), I(1), dual(I(1)), dual(I(0)))][2, 1, 1, 2] = -1
    return t
end
function c_plusmin_down(T, ::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    t = two_site_operator(T, Trivial, U1Irrep)
    I = sectortype(t)
    t[(I(1, -1 // 2), I(0, 0), dual(I(0, 0)), dual(I(1, -1 // 2)))][1, 1, 1, 1] = 1
    t[(I(1, -1 // 2), I(1, 1 // 2), dual(I(0, 0)), dual(I(0, 0)))][1, 1, 1, 2] = -1
    t[(I(0, 0), I(0, 0), dual(I(1, 1 // 2)), dual(I(1, -1 // 2)))][2, 1, 1, 1] = 1
    t[(I(0, 0), I(1, 1 // 2), dual(I(1, 1 // 2)), dual(I(0, 0)))][2, 1, 1, 2] = -1
    return t
end
function c_plusmin_down(T, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    return error("Not implemented")
end
function c_plusmin_down(T, ::Type{U1Irrep}, ::Type{Trivial}; filling::Tuple{Int64,Int64}=(1,1))
    t = two_site_operator(T, U1Irrep, Trivial; filling=filling)
    P, Q = filling
    I = sectortype(t)
    t[(I(1, Q-P), I(0, -P), dual(I(0, -P)), dual(I(1, Q-P)))][2, 1, 1, 2] = 1
    t[(I(1, Q-P), I(1, Q-P), dual(I(0, -P)), dual(I(0, 2Q-P)))][2, 1, 1, 1] = -1
    t[(I(0, 2Q-P), I(0, -P), dual(I(1, Q-P)), dual(I(1, Q-P)))][1, 1, 1, 2] = 1
    t[(I(0, 2Q-P), I(1, Q-P), dual(I(1, Q-P)), dual(I(0, 2Q-P)))][1, 1, 1, 1] = -1
    return t
end
function c_plusmin_down(T, ::Type{U1Irrep}, ::Type{U1Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    t = two_site_operator(T, U1Irrep, U1Irrep; filling=filling)
    P, Q = filling
    I = sectortype(t)
    t[(I(1, Q-P, -1 // 2), I(0, -P, 0), dual(I(0, -P, 0)), dual(I(1, Q-P, -1 // 2)))] .= 1
    t[(I(1, Q-P, -1 // 2), I(1, Q-P, 1 // 2), dual(I(0, -P, 0)), dual(I(0, 2Q-P, 0)))] .= -1
    t[(I(0, 2Q-P, 0), I(0, -P, 0), dual(I(1, Q-P, 1 // 2)), dual(I(1, Q-P, -1 // 2)))] .= 1
    t[(I(0, 2Q-P, 0), I(1, Q-P, 1 // 2), dual(I(1, Q-P, 1 // 2)), dual(I(0, 2Q-P, 0)))] .= -1
    return t
end
function c_plusmin_down(T, ::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    return error("Not implemented")
end
function c_plusmin_down(T, ::Type{SU2Irrep}, ::Type{Trivial}; kwargs...)
    return error("Not implemented")
end
function c_plusmin_down(T, ::Type{SU2Irrep}, ::Type{U1Irrep}; kwargs...)
    return error("Not implemented")
end
function c_plusmin_down(T, ::Type{SU2Irrep}, ::Type{SU2Irrep}; kwargs...)
    return error("Not implemented")
end

"""
    c_minplus_up(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the Hermitian conjugate of `c_plusmin_up`, i.e.
``(c†_{1,↑}, c_{2,↑})† = -c_{1,↑}, c†_{2,↑}`` (note the extra minus sign). 
It annihilates a spin-up electron at the first site and creates a spin-up electron at the second.
"""
c_minplus_up(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_minplus_up(ComplexF64, P, S; kwargs...)
function c_minplus_up(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; filling::Tuple{Int64,Int64}=(1,1))
    return copy(adjoint(c_plusmin_up(T, particle_symmetry, spin_symmetry; filling=filling)))
end

"""
    c_minplus_down(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the Hermitian conjugate of `c_plusmin_down`, i.e.
``(c†_{1,↓}, c_{2,↓})† = -c_{1,↓}, c†_{2,↓}`` (note the extra minus sign). 
It annihilates a spin-down electron at the first site and creates a spin-down electron at the second.
"""
c_minplus_down(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_minplus_down(ComplexF64, P, S; kwargs...)
function c_minplus_down(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; filling::Tuple{Int64,Int64}=(1,1))
    return copy(adjoint(c_plusmin_down(T, particle_symmetry, spin_symmetry; filling=filling)))
end

"""
    c_plusmin(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the two-body operator that creates a particle at the first site and annihilates a particle at the second.
This is the sum of `c_plusmin_up` and `c_plusmin_down`.
"""
c_plusmin(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_plusmin(ComplexF64, P, S; kwargs...)
function c_plusmin(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; filling::Tuple{Int64,Int64}=(1,1))
    return c_plusmin_up(T, particle_symmetry, spin_symmetry; filling=filling) +
        c_plusmin_down(T, particle_symmetry, spin_symmetry; filling=filling)
end
function c_plusmin(T, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    t = two_site_operator(T, Trivial, SU2Irrep)
    I = sectortype(t)
    f1 = only(fusiontrees((I(0, 0), I(1, 1 // 2)), I(1, 1 // 2)))
    f2 = only(fusiontrees((I(1, 1 // 2), I(0, 0)), I(1, 1 // 2)))
    t[f1, f2][1, 1, 1, 1] = 1
    f3 = only(fusiontrees((I(1, 1 // 2), I(0, 0)), I(1, 1 // 2)))
    f4 = only(fusiontrees((I(0, 0), I(1, 1 // 2)), I(1, 1 // 2)))
    t[f3, f4][1, 2, 2, 1] = -1
    f5 = only(fusiontrees((I(0, 0), I(0, 0)), I(0, 0)))
    f6 = only(fusiontrees((I(1, 1 // 2), I(1, 1 // 2)), I(0, 0)))
    t[f5, f6][1, 2, 1, 1] = sqrt(2)
    f7 = only(fusiontrees((I(1, 1 // 2), I(1, 1 // 2)), I(0, 0)))
    f8 = only(fusiontrees((I(0, 0), I(0, 0)), I(0, 0)))
    t[f7, f8][1, 1, 2, 1] = sqrt(2)
    return t
end
function c_plusmin(T, ::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    P, Q = filling
    t = two_site_operator(T, U1Irrep, SU2Irrep; filling=filling)
    I = sectortype(t)

    # t = two_site_operator(T, U1Irrep, SU2Irrep; filling=filling)
    # t[(I(1, Q-P, 1 // 2), I(0, -P, 0), dual(I(0, -P, 0)), dual(I(1, Q-P, 1 // 2)))] .= 1
    # t[(I(1, Q-P, 1 // 2), I(1, Q-P, 1 // 2), dual(I(0, -P, 0)), dual(I(0, 2Q-P, 0)))] .= 1
    # t[(I(0, 2Q-P, 0), I(0, -P, 0), dual(I(1, Q-P, 1 // 2)), dual(I(1, Q-P, 1 // 2)))] .= 1
    # t[(I(0, 2Q-P, 0), I(1, Q-P, 1 // 2), dual(I(1, Q-P, 1 // 2)), dual(I(0, 2Q-P, 0)))] .= 1

    Ps = hubbard_space(U1Irrep, SU2Irrep; filling=filling)
    Vs = Vect[I]((1, Q, 1 // 2) => 1)

    c_plus = zeros(T, Ps ← Ps ⊗ Vs)
    blocks(c_plus)[I((1, Q-P, 1 // 2))] .= 1
    blocks(c_plus)[I((0, 2Q-P, 0))] .= sqrt(2)

    c_min = zeros(T, Vs ⊗ Ps ← Ps)
    blocks(c_min)[I((1, Q-P, 1 // 2))] .= 1
    blocks(c_min)[I((0, 2Q-P, 0))] .= sqrt(2)

    @planar twosite[-1 -2; -3 -4] := c_plus[-1; -3 1] * c_min[1 -2; -4]
    return twosite
end

"""
    c_minplus(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the two-body operator that annihilates a particle at the first site and creates a particle at the second.
This is the sum of `c_minplus_up` and `c_minplus_down`.
"""
c_minplus(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_minplus(ComplexF64, P, S; kwargs...)
function c_minplus(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; filling::Tuple{Int64,Int64}=(1,1))
    return copy(adjoint(c_plusmin(T, particle_symmetry, spin_symmetry; filling=filling)))
end

"""
    number_up(particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the one-body operator that counts the number of spin-up electrons.
"""
number_up(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = number_up(ComplexF64, P, S; kwargs...)
function number_up(T::Type{<:Number}, ::Type{Trivial} = Trivial, ::Type{Trivial} = Trivial; kwargs...)
    t = single_site_operator(T, Trivial, Trivial)
    I = sectortype(t)
    t[(I(1), I(1))][1, 1] = 1
    t[(I(0), I(0))][2, 2] = 1
    return t
end
function number_up(T, ::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, U1Irrep)
    I = sectortype(t)
    t[(I(1, 1 // 2), dual(I(1, 1 // 2)))][1, 1] = 1
    t[(I(0, 0), dual(I(0, 0)))][2, 2] = 1
    return t
end
function number_up(T, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    throw(ArgumentError("`number_up` is not symmetric under `SU2Irrep` spin symmetry"))
end
function number_up(T, ::Type{U1Irrep}, ::Type{Trivial}; filling::Tuple{Int64,Int64}=(1,1))
    t = single_site_operator(T, U1Irrep, Trivial; filling=filling)
    P, Q = filling
    I = sectortype(t)
    block(t, I(1, Q-P))[1, 1] = 1
    block(t, I(0, 2Q-P))[1, 1] = 1
    return t
end
function number_up(T, ::Type{U1Irrep}, ::Type{U1Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    t = single_site_operator(T, U1Irrep, U1Irrep; filling=filling)
    P, Q = filling
    I = sectortype(t)
    block(t, I(1, Q-P, 1 // 2)) .= 1
    block(t, I(0, 2Q-P, 0)) .= 1
    return t
end
function number_up(T, ::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    throw(ArgumentError("`number_up` is not symmetric under `SU2Irrep` spin symmetry"))
end
function number_up(T, ::Type{SU2Irrep}, ::Type{Trivial}; kwargs...)
    return error("Not implemented")
end
function number_up(T, ::Type{SU2Irrep}, ::Type{U1Irrep}; kwargs...)
    return error("Not implemented")
end
function number_up(T, ::Type{SU2Irrep}, ::Type{SU2Irrep}; kwargs...)
    throw(ArgumentError("`number_up` is not symmetric under `SU2Irrep` spin symmetry"))
end

"""
    number_down(particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the one-body operator that counts the number of spin-down electrons.
"""
number_down(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = number_down(ComplexF64, P, S; kwargs...)
function number_down(T::Type{<:Number}, ::Type{Trivial} = Trivial, ::Type{Trivial} = Trivial; kwargs...)
    t = single_site_operator(T, Trivial, Trivial)
    I = sectortype(t)
    t[(I(1), I(1))][2, 2] = 1
    t[(I(0), I(0))][2, 2] = 1
    return t
end
function number_down(T, ::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, U1Irrep)
    I = sectortype(t)
    t[(I(1, -1 // 2), dual(I(1, -1 // 2)))][1, 1] = 1
    t[(I(0, 0), I(0, 0))][2, 2] = 1
    return t
end
function number_down(T, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    throw(ArgumentError("`number_down` is not symmetric under `SU2Irrep` spin symmetry"))
end
function number_down(T, ::Type{U1Irrep}, ::Type{Trivial}; filling::Tuple{Int64,Int64}=(1,1))
    t = single_site_operator(T, U1Irrep, Trivial; filling=filling)
    P, Q = filling
    I = sectortype(t)
    block(t, I(1, Q-P))[2, 2] = 1 # expected to be [1,2]
    block(t, I(0, 2Q-P))[1, 1] = 1
    return t
end
function number_down(T, ::Type{U1Irrep}, ::Type{U1Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    t = single_site_operator(T, U1Irrep, U1Irrep; filling=filling)
    P, Q = filling
    I = sectortype(t)
    block(t, I(1, Q-P, -1 // 2)) .= 1
    block(t, I(0, 2Q-P, 0)) .= 1
    return t
end
function number_down(T, ::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    throw(ArgumentError("`number_down` is not symmetric under `SU2Irrep` spin symmetry"))
end
function number_down(T, ::Type{SU2Irrep}, ::Type{Trivial}; kwargs...)
    return error("Not implemented")
end
function number_down(T, ::Type{SU2Irrep}, ::Type{U1Irrep}; kwargs...)
    return error("Not implemented")
end
function number_down(T, ::Type{SU2Irrep}, ::Type{SU2Irrep}; kwargs)
    throw(ArgumentError("`number_down` is not symmetric under `SU2Irrep` spin symmetry"))
end

"""
    number_e(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the one-body operator that counts the number of particles.
"""
number_e(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = number_e(ComplexF64, P, S; kwargs...)
function number_e(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return number_up(T, particle_symmetry, spin_symmetry; kwargs...) +
        number_down(T, particle_symmetry, spin_symmetry; kwargs...)
end
function number_e(T, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, SU2Irrep; kwargs...)
    I = sectortype(t)
    block(t, I(1, 1 // 2))[1, 1] = 1
    block(t, I(0, 0))[2, 2] = 2
    return t
end
function number_e(T, ::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    t = single_site_operator(T, U1Irrep, SU2Irrep; filling=filling)
    P, Q = filling
    I = sectortype(t)
    block(t, I(1, Q-P, 1 // 2)) .= 1
    block(t, I(0, 2Q-P, 0)) .= 2
    return t
end

"""
    number_pair(T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the one-body operator that counts the number of doubly occupied sites.
"""
number_pair(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = number_pair(ComplexF64, P, S; kwargs...)
function number_pair(
        T, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...
    )
    return number_up(T, particle_symmetry, spin_symmetry; kwargs...) *
        number_down(T, particle_symmetry, spin_symmetry; kwargs...)
end
function number_pair(T, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, SU2Irrep)
    I = sectortype(t)
    block(t, I(0, 0))[2, 2] = 1
    return t
end
function number_pair(T, ::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Tuple{Int64,Int64}=(1,1))
    t = single_site_operator(T, U1Irrep, SU2Irrep; filling=filling)
    P, Q = filling
    I = sectortype(t)
    block(t, I(0, 2Q-P, 0)) .= 1
    return t
end
