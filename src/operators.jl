##########
# Spaces #
##########

"""
    hubbard_space(particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; filling::Rational{Int}=1//1)

Construct the local Hilbert space for a Hubbard-type model with the specified particle-number
and spin symmetries.

Supported symmetries are `Trivial` and `U1Irrep` for both particle number and spin, with
`SU2Irrep` additionally supported for spin. When `particle_symmetry` is `U1Irrep`, the
filling can be specified as a rational number `P//Q`, representing the number of particles
per site. The default is `1//1`, corresponding to half-filling.
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
function hubbard_space(::Type{U1Irrep}, ::Type{Trivial}; filling::Rational{Int}=1//1)
    P = numerator(filling); Q = denominator(filling)
    return Vect[FermionParity ⊠ U1Irrep]((0, -P) => 1, (1, Q-P) => 2, (0, 2Q-P) => 1)
end
function hubbard_space(::Type{U1Irrep}, ::Type{U1Irrep}; filling::Rational{Int}=1//1)
    P = numerator(filling); Q = denominator(filling)
    return Vect[FermionParity ⊠ U1Irrep ⊠ U1Irrep](
        (0, -P, 0) => 1, (1, Q-P, 1 // 2) => 1,
        (1, Q-P, -1 // 2) => 1, (0, 2Q-P, 0) => 1
    )
end
function hubbard_space(::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Rational{Int}=1//1)
    P = numerator(filling); Q = denominator(filling)
    return Vect[FermionParity ⊠ U1Irrep ⊠ SU2Irrep](
        (0, -P, 0) => 1, (1, Q-P, 1 // 2) => 1, (0, 2Q-P, 0) => 1
    )
end


#############
# Operators #
#############

function single_site_operator(
        T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; filling::Rational{Int}=1//1
    )
    V = hubbard_space(particle_symmetry, spin_symmetry; filling=filling)
    return zeros(T, V ← V)
end

function two_site_operator(
        T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; filling::Rational{Int}=1//1
    )
    V = hubbard_space(particle_symmetry, spin_symmetry; filling=filling)
    return zeros(T, V ⊗ V ← V ⊗ V)
end

function boson_single_site_operator(
        T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}, max_b::Int64
    )
    V = holstein_space(particle_symmetry, spin_symmetry, max_b)
    return zeros(T, V ← V)
end

"""
    c_plusmin_up(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the two-site operator ``c†_{1,↑} c_{2,↑}`` that creates a spin-up electron at the first site and annihilates one at the second.
"""
c_plusmin_up(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_plusmin_up(ComplexF64, P, S; kwargs...)
function c_plusmin_up(T::Type{<:Number}, ::Type{Trivial}, ::Type{Trivial}; kwargs...)
    t = two_site_operator(T, Trivial, Trivial)
    I = sectortype(t)
    t[(I(1), I(0), dual(I(0)), dual(I(1)))][1, 1, 1, 1] = 1
    t[(I(1), I(1), dual(I(0)), dual(I(0)))][1, 2, 1, 2] = 1
    t[(I(0), I(0), dual(I(1)), dual(I(1)))][2, 1, 2, 1] = -1
    t[(I(0), I(1), dual(I(1)), dual(I(0)))][2, 2, 2, 2] = -1
    return t
end
function c_plusmin_up(T::Type{<:Number}, ::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    t = two_site_operator(T, Trivial, U1Irrep)
    I = sectortype(t)
    t[(I(1, 1 // 2), I(0, 0), dual(I(0, 0)), dual(I(1, 1 // 2)))][1, 1, 1, 1] = 1
    t[(I(1, 1 // 2), I(1, -1 // 2), dual(I(0, 0)), dual(I(0, 0)))][1, 1, 1, 2] = 1
    t[(I(0, 0), I(0, 0), dual(I(1, -1 // 2)), dual(I(1, 1 // 2)))][2, 1, 1, 1] = -1
    t[(I(0, 0), I(1, -1 // 2), dual(I(1, -1 // 2)), dual(I(0, 0)))][2, 1, 1, 2] = -1
    return t
end
function c_plusmin_up(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{Trivial}; filling::Rational{Int}=1//1)
    t = two_site_operator(T, U1Irrep, Trivial; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    t[(I(1, Q-P), I(0, -P), dual(I(0, -P)), dual(I(1, Q-P)))][1, 1, 1, 1] = 1
    t[(I(1, Q-P), I(1, Q-P), dual(I(0, -P)), dual(I(0, 2Q-P)))][1, 2, 1, 1] = 1
    t[(I(0, 2Q-P), I(0, -P), dual(I(1, Q-P)), dual(I(1, Q-P)))][1, 1, 2, 1] = -1
    t[(I(0, 2Q-P), I(1, Q-P), dual(I(1, Q-P)), dual(I(0, 2Q-P)))][1, 2, 2, 1] = -1
    return t
end
function c_plusmin_up(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{U1Irrep}; filling::Rational{Int}=1//1)
    t = two_site_operator(T, U1Irrep, U1Irrep; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    t[(I(1, Q-P, 1 // 2), I(0, -P, 0), dual(I(0, -P, 0)), dual(I(1, Q-P, 1 // 2)))] .= 1
    t[(I(1, Q-P, 1 // 2), I(1, Q-P, -1 // 2), dual(I(0, -P, 0)), dual(I(0, 2Q-P, 0)))] .= 1
    t[(I(0, 2Q-P, 0), I(0, -P, 0), dual(I(1, Q-P, -1 // 2)), dual(I(1, Q-P, 1 // 2)))] .= -1
    t[(I(0, 2Q-P, 0), I(1, Q-P, -1 // 2), dual(I(1, Q-P, -1 // 2)), dual(I(0, 2Q-P, 0)))] .= -1
    return t
end
function c_plusmin_up(T::Type{<:Number}, ::Type{<:Sector}, ::Type{SU2Irrep}; filling::Rational{Int}=1//1)
    throw(ArgumentError("`c_plusmin_up` is not symmetric under `SU2Irrep` spin symmetry"))
end

"""
    c_plusmin_down(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the two-site operator ``c†_{1,↓} c_{2,↓}`` that creates a spin-down electron at the first site and annihilates one at the second.
"""
c_plusmin_down(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_plusmin_down(ComplexF64, P, S; kwargs...)
function c_plusmin_down(T::Type{<:Number}, ::Type{Trivial}, ::Type{Trivial}; kwargs...)
    t = two_site_operator(T, Trivial, Trivial)
    I = sectortype(t)
    t[(I(1), I(0), dual(I(0)), dual(I(1)))][2, 1, 1, 2] = 1
    t[(I(1), I(1), dual(I(0)), dual(I(0)))][2, 1, 1, 2] = -1
    t[(I(0), I(0), dual(I(1)), dual(I(1)))][2, 1, 1, 2] = 1
    t[(I(0), I(1), dual(I(1)), dual(I(0)))][2, 1, 1, 2] = -1
    return t
end
function c_plusmin_down(T::Type{<:Number}, ::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    t = two_site_operator(T, Trivial, U1Irrep)
    I = sectortype(t)
    t[(I(1, -1 // 2), I(0, 0), dual(I(0, 0)), dual(I(1, -1 // 2)))][1, 1, 1, 1] = 1
    t[(I(1, -1 // 2), I(1, 1 // 2), dual(I(0, 0)), dual(I(0, 0)))][1, 1, 1, 2] = -1
    t[(I(0, 0), I(0, 0), dual(I(1, 1 // 2)), dual(I(1, -1 // 2)))][2, 1, 1, 1] = 1
    t[(I(0, 0), I(1, 1 // 2), dual(I(1, 1 // 2)), dual(I(0, 0)))][2, 1, 1, 2] = -1
    return t
end
function c_plusmin_down(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{Trivial}; filling::Rational{Int}=1//1)
    t = two_site_operator(T, U1Irrep, Trivial; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    t[(I(1, Q-P), I(0, -P), dual(I(0, -P)), dual(I(1, Q-P)))][2, 1, 1, 2] = 1
    t[(I(1, Q-P), I(1, Q-P), dual(I(0, -P)), dual(I(0, 2Q-P)))][2, 1, 1, 1] = -1
    t[(I(0, 2Q-P), I(0, -P), dual(I(1, Q-P)), dual(I(1, Q-P)))][1, 1, 1, 2] = 1
    t[(I(0, 2Q-P), I(1, Q-P), dual(I(1, Q-P)), dual(I(0, 2Q-P)))][1, 1, 1, 1] = -1
    return t
end
function c_plusmin_down(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{U1Irrep}; filling::Rational{Int}=1//1)
    t = two_site_operator(T, U1Irrep, U1Irrep; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    t[(I(1, Q-P, -1 // 2), I(0, -P, 0), dual(I(0, -P, 0)), dual(I(1, Q-P, -1 // 2)))] .= 1
    t[(I(1, Q-P, -1 // 2), I(1, Q-P, 1 // 2), dual(I(0, -P, 0)), dual(I(0, 2Q-P, 0)))] .= -1
    t[(I(0, 2Q-P, 0), I(0, -P, 0), dual(I(1, Q-P, 1 // 2)), dual(I(1, Q-P, -1 // 2)))] .= 1
    t[(I(0, 2Q-P, 0), I(1, Q-P, 1 // 2), dual(I(1, Q-P, 1 // 2)), dual(I(0, 2Q-P, 0)))] .= -1
    return t
end
function c_plusmin_down(T::Type{<:Number}, ::Type{<:Sector}, ::Type{SU2Irrep}; filling::Rational{Int}=1//1)
    throw(ArgumentError("`c_plusmin_down` is not symmetric under `SU2Irrep` spin symmetry"))
end

"""
    c_minplus_up(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the Hermitian conjugate of `c_plusmin_up`, i.e.
``(c†_{1,↑} c_{2,↑})† = -c_{1,↑} c†_{2,↑}`` (note the extra minus sign).
It annihilates a spin-up electron at the first site and creates one at the second.
"""
c_minplus_up(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_minplus_up(ComplexF64, P, S; kwargs...)
function c_minplus_up(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return copy(adjoint(c_plusmin_up(T, particle_symmetry, spin_symmetry; kwargs...)))
end

"""
    c_minplus_down(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the Hermitian conjugate of `c_plusmin_down`, i.e.
``(c†_{1,↓} c_{2,↓})† = -c_{1,↓} c†_{2,↓}`` (note the extra minus sign).
It annihilates a spin-down electron at the first site and creates one at the second.
"""
c_minplus_down(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_minplus_down(ComplexF64, P, S; kwargs...)
function c_minplus_down(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return copy(adjoint(c_plusmin_down(T, particle_symmetry, spin_symmetry; kwargs...)))
end

"""
    c_plusmin_updown(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the two-site spin-flip operator ``c†_{1,↑} c_{2,↓}``.
It is only defined when the spin symmetry is `Trivial`.
"""
c_plusmin_updown(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_plusmin_updown(ComplexF64, P, S; kwargs...)
function c_plusmin_updown(T::Type{<:Number}, ::Type{Trivial}, ::Type{Trivial}; kwargs...)
    t = two_site_operator(T, Trivial, Trivial)
    I = sectortype(t)
    #  I(0) (even): 1 = |0⟩,  2 = |↑↓⟩
    #  I(1) (odd) : 1 = |↑⟩,  2 = |↓⟩
    # |0,↓⟩ -> |↑,0⟩
    t[(I(1), I(0), dual(I(0)), dual(I(1)))][1, 1, 1, 2] = 1
    # |0,↑↓⟩ -> -|↑,↑⟩
    t[(I(1), I(1), dual(I(0)), dual(I(0)))][1, 1, 1, 2] = -1
    # |↓,↓⟩ -> -|↑↓,0⟩
    t[(I(0), I(0), dual(I(1)), dual(I(1)))][2, 1, 2, 2] = -1
    # |↓,↑↓⟩ -> |↑↓,↑⟩
    t[(I(0), I(1), dual(I(1)), dual(I(0)))][2, 1, 2, 2] = 1
    
    return t
end
function c_plusmin_updown(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{Trivial}; filling::Rational{Int}=1//1)
    t = two_site_operator(T, U1Irrep, Trivial; filling=filling)
    P = numerator(filling)
    Q = denominator(filling)
    I = sectortype(t)
    # 1 = ↑, 2 = ↓
    # |0,↓> -> |↑,0>
    t[(I(1, Q-P), I(0, -P), dual(I(0, -P)), dual(I(1, Q-P)))][1, 1, 1, 2] = 1
    # |0,↑↓> -> - |↑,↑>
    t[(I(1, Q-P), I(1, Q-P), dual(I(0, -P)), dual(I(0, 2Q-P)))][1, 1, 1, 1] = -1
    # |↓,↓> -> - |↑↓,0>
    t[(I(0, 2Q-P), I(0, -P), dual(I(1, Q-P)), dual(I(1, Q-P)))][1, 1, 2, 2] = -1
    # |↓,↑↓> -> |↑↓,↑>
    t[(I(0, 2Q-P), I(1, Q-P), dual(I(1, Q-P)), dual(I(0, 2Q-P)))][1, 1, 2, 1] = 1
    return t
end
function c_plusmin_updown(T::Type{<:Number}, ::Type{<:Sector}, ::Type{U1Irrep}; kwargs...)
    throw(ArgumentError("`c_plusmin_updown` is not symmetric under `U1Irrep` spin symmetry"))
end
function c_plusmin_updown(T::Type{<:Number}, ::Type{<:Sector}, ::Type{SU2Irrep}; kwargs...)
    throw(ArgumentError("`c_plusmin_updown` is not symmetric under `SU2Irrep` spin symmetry"))
end

"""
    c_minplus_updown(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the adjoint of `c_plusmin_updown`, namely the reverse spin-flip operator
``c†_{2,↓} c_{1,↑}``.
"""
c_minplus_updown(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_minplus_updown(ComplexF64, P, S; kwargs...)
function c_minplus_updown(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return copy(adjoint(c_plusmin_updown(T, particle_symmetry, spin_symmetry; kwargs...)))
end

"""
    c_plusmin_downup(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the two-site spin-flip operator ``c†_{1,↓} c_{2,↑}``.
It is only defined when the spin symmetry is `Trivial`.
"""
c_plusmin_downup(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_plusmin_downup(ComplexF64, P, S; kwargs...)
function c_plusmin_downup(T::Type{<:Number}, ::Type{Trivial}, ::Type{Trivial}; kwargs...)
    t = two_site_operator(T, Trivial, Trivial)
    I = sectortype(t)
    #  I(0) (even): 1 = |0⟩,  2 = |↑↓⟩
    #  I(1) (odd) : 1 = |↑⟩,  2 = |↓⟩
    # |0,↑> -> |↓,0>
    t[(I(1), I(0), dual(I(0)), dual(I(1)))][2, 1, 1, 1] = 1
    # |0,↑↓> -> |↓,↓>
    t[(I(1), I(1), dual(I(0)), dual(I(0)))][2, 2, 1, 2] = 1
    # |↑,↑> -> |↑↓,0>
    t[(I(0), I(0), dual(I(1)), dual(I(1)))][2, 1, 1, 1] = 1
    # |↑,↑↓> -> |↑↓,↓>
    t[(I(0), I(1), dual(I(1)), dual(I(0)))][2, 2, 1, 2] = 1
    
    return t
end
function c_plusmin_downup(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{Trivial}; filling::Rational{Int}=1//1)
    t = two_site_operator(T, U1Irrep, Trivial; filling=filling)
    P = numerator(filling)
    Q = denominator(filling)
    I = sectortype(t)
    # 1 = ↑, 2 = ↓
    # |0,↑> -> |↓,0>
    t[(I(1, Q-P), I(0, -P), dual(I(0, -P)), dual(I(1, Q-P)))][2, 1, 1, 1] = 1
    # |0,↑↓> -> |↓,↓>
    t[(I(1, Q-P), I(1, Q-P), dual(I(0, -P)), dual(I(0, 2Q-P)))][2, 2, 1, 1] = 1
    # |↑,↑> -> |↑↓,0>
    t[(I(0, 2Q-P), I(0, -P), dual(I(1, Q-P)), dual(I(1, Q-P)))][1, 1, 1, 1] = 1
    # |↑,↑↓> -> |↑↓,↓>
    t[(I(0, 2Q-P), I(1, Q-P), dual(I(1, Q-P)), dual(I(0, 2Q-P)))][1, 2, 1, 1] = 1
    return t
end
function c_plusmin_downup(T::Type{<:Number}, ::Type{<:Sector}, ::Type{U1Irrep}; kwargs...)
    throw(ArgumentError("`c_plusmin_downup` is not symmetric under `U1Irrep` spin symmetry"))
end
function c_plusmin_downup(T::Type{<:Number}, ::Type{<:Sector}, ::Type{SU2Irrep}; kwargs...)
    throw(ArgumentError("`c_plusmin_downup` is not symmetric under `SU2Irrep` spin symmetry"))
end

"""
    c_minplus_downup(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the adjoint of `c_plusmin_downup`, namely the reverse spin-flip operator
``c†_{2,↑} c_{1,↓}``.
"""
c_minplus_downup(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_minplus_downup(ComplexF64, P, S; kwargs...)
function c_minplus_downup(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return copy(adjoint(c_plusmin_downup(T, particle_symmetry, spin_symmetry; kwargs...)))
end

"""
    c_plusmin(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the two-site operator that creates a particle at the first site and annihilates one at the second.
This is the sum of `c_plusmin_up` and `c_plusmin_down`.
"""
c_plusmin(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_plusmin(ComplexF64, P, S; kwargs...)
function c_plusmin(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return c_plusmin_up(T, particle_symmetry, spin_symmetry; kwargs...) +
        c_plusmin_down(T, particle_symmetry, spin_symmetry; kwargs...)
end
function c_plusmin(T::Type{<:Number}, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    t = two_site_operator(T, Trivial, SU2Irrep)
    I = sectortype(t)
    f1 = only(fusiontrees((I(0, 0), I(1, 1 // 2)), I(1, 1 // 2)))
    f2 = only(fusiontrees((I(1, 1 // 2), I(0, 0)), I(1, 1 // 2)))
    t[f2, f1][1, 1, 1, 1] = 1
    f3 = only(fusiontrees((I(1, 1 // 2), I(0, 0)), I(1, 1 // 2)))
    f4 = only(fusiontrees((I(0, 0), I(1, 1 // 2)), I(1, 1 // 2)))
    t[f4, f3][2, 1, 1, 2] = -1
    f5 = only(fusiontrees((I(0, 0), I(0, 0)), I(0, 0)))
    f6 = only(fusiontrees((I(1, 1 // 2), I(1, 1 // 2)), I(0, 0)))
    t[f6, f5][1, 1, 1, 2] = sqrt(2)
    f7 = only(fusiontrees((I(1, 1 // 2), I(1, 1 // 2)), I(0, 0)))
    f8 = only(fusiontrees((I(0, 0), I(0, 0)), I(0, 0)))
    t[f8, f7][2, 1, 1, 1] = sqrt(2)
    return t
end
function c_plusmin(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Rational{Int}=1//1)
    P = numerator(filling); Q = denominator(filling)
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
    c_minplus(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the two-site operator that annihilates a particle at the first site and creates one at the second.
This is the sum of `c_minplus_up` and `c_minplus_down`.
"""
c_minplus(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = c_minplus(ComplexF64, P, S; kwargs...)
function c_minplus(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; filling::Rational{Int}=1//1)
    return copy(adjoint(c_plusmin(T, particle_symmetry, spin_symmetry; filling=filling)))
end

"""
    create_pair_onesite(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the one-body onsite pair creation operator Δ† = c†_↑ c†_↓.
It maps the empty state |0⟩ to the doubly occupied state |↑↓⟩.
"""
create_pair_onesite(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = create_pair_onesite(ComplexF64, P, S; kwargs...)
function create_pair_onesite(T::Type{<:Number}, ::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, U1Irrep)
    I = sectortype(t)
    t[(I(0, 0), dual(I(0, 0)))][2, 1] = 1
    return t
end
function create_pair_onesite(T::Type{<:Number}, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, SU2Irrep)
    I = sectortype(t)
    block(t, I(0, 0))[2, 1] = 1
    return t
end
function create_pair_onesite(T::Type{<:Number}, ::Type{Trivial}, ::Type{Trivial}; kwargs...)
    t = single_site_operator(T, Trivial, Trivial)
    I = sectortype(t)
    t[(I(0), dual(I(0)))][2, 1] = 1
    return t
end
function create_pair_onesite(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{<:Sector}; kwargs...)
    throw(ArgumentError("`create_pair_onesite` is not symmetric under `U1Irrep` particle symmetry"))
end

"""
    delete_pair_onesite(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the one-body onsite pair annihilation operator Δ = c_↓ c_↑.
It maps the doubly occupied state |↑↓⟩ to the empty state |0⟩.
"""
delete_pair_onesite(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = delete_pair_onesite(ComplexF64, P, S; kwargs...)
function delete_pair_onesite(T::Type{<:Number}, ::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, U1Irrep)
    I = sectortype(t)
    t[(I(0, 0), dual(I(0, 0)))][1, 2] = 1
    return t
end
function delete_pair_onesite(T::Type{<:Number}, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, SU2Irrep)
    I = sectortype(t)
    block(t, I(0, 0))[1,2] = 1
    return t
end
function delete_pair_onesite(T::Type{<:Number}, ::Type{Trivial}, ::Type{Trivial}; kwargs...)
    t = single_site_operator(T, Trivial, Trivial)
    I = sectortype(t)
    t[(I(0), dual(I(0)))][1, 2] = 1
    return t
end
function delete_pair_onesite(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{<:Sector}; kwargs...)
    throw(ArgumentError("`delete_pair_onesite` is not symmetric under `U1Irrep` particle symmetry"))
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
function number_up(T::Type{<:Number}, ::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, U1Irrep)
    I = sectortype(t)
    t[(I(1, 1 // 2), dual(I(1, 1 // 2)))][1, 1] = 1
    t[(I(0, 0), dual(I(0, 0)))][2, 2] = 1
    return t
end
function number_up(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{Trivial}; filling::Rational{Int}=1//1)
    t = single_site_operator(T, U1Irrep, Trivial; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    block(t, I(1, Q-P))[1, 1] = 1
    block(t, I(0, 2Q-P))[1, 1] = 1
    return t
end
function number_up(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{U1Irrep}; filling::Rational{Int}=1//1)
    t = single_site_operator(T, U1Irrep, U1Irrep; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    block(t, I(1, Q-P, 1 // 2)) .= 1
    block(t, I(0, 2Q-P, 0)) .= 1
    return t
end
function number_up(T::Type{<:Number}, ::Type{<:Sector}, ::Type{SU2Irrep}; filling::Rational{Int}=1//1)
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
function number_down(T::Type{<:Number}, ::Type{Trivial}, ::Type{U1Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, U1Irrep)
    I = sectortype(t)
    t[(I(1, -1 // 2), dual(I(1, -1 // 2)))][1, 1] = 1
    t[(I(0, 0), I(0, 0))][2, 2] = 1
    return t
end
function number_down(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{Trivial}; filling::Rational{Int}=1//1)
    t = single_site_operator(T, U1Irrep, Trivial; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    block(t, I(1, Q-P))[2, 2] = 1
    block(t, I(0, 2Q-P))[1, 1] = 1
    return t
end
function number_down(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{U1Irrep}; filling::Rational{Int}=1//1)
    t = single_site_operator(T, U1Irrep, U1Irrep; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    block(t, I(1, Q-P, -1 // 2)) .= 1
    block(t, I(0, 2Q-P, 0)) .= 1
    return t
end
function number_down(T::Type{<:Number}, ::Type{<:Sector}, ::Type{SU2Irrep}; filling::Rational{Int}=1//1)
    throw(ArgumentError("`number_down` is not symmetric under `SU2Irrep` spin symmetry"))
end

"""
    number_e(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the one-body operator that counts the number of particles.
"""
number_e(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = number_e(ComplexF64, P, S; kwargs...)
function number_e(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return number_up(T, particle_symmetry, spin_symmetry; kwargs...) +
        number_down(T, particle_symmetry, spin_symmetry; kwargs...)
end
function number_e(T::Type{<:Number}, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, SU2Irrep)
    I = sectortype(t)
    block(t, I(1, 1 // 2))[1, 1] = 1
    block(t, I(0, 0))[2, 2] = 2
    return t
end
function number_e(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Rational{Int}=1//1)
    t = single_site_operator(T, U1Irrep, SU2Irrep; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    block(t, I(1, Q-P, 1 // 2)) .= 1
    block(t, I(0, 2Q-P, 0)) .= 2
    return t
end

"""
    number_pair(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the one-body operator that counts the number of doubly occupied sites.
"""
number_pair(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = number_pair(ComplexF64, P, S; kwargs...)
function number_pair(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return number_up(T, particle_symmetry, spin_symmetry; kwargs...) *
        number_down(T, particle_symmetry, spin_symmetry; kwargs...)
end
function number_pair(T::Type{<:Number}, ::Type{Trivial}, ::Type{SU2Irrep}; kwargs...)
    t = single_site_operator(T, Trivial, SU2Irrep)
    I = sectortype(t)
    block(t, I(0, 0))[2, 2] = 1
    return t
end
function number_pair(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{SU2Irrep}; filling::Rational{Int}=1//1)
    t = single_site_operator(T, U1Irrep, SU2Irrep; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    block(t, I(0, 2Q-P, 0)) .= 1
    return t
end

"""
    S_plus(T::Type{<:Number}, P::Type{<:Sector}, S::Type{<:Sector})

Return the spin-raising operator ``S^+``.
"""
S_plus(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = S_plus(ComplexF64, P, S; kwargs...)
function S_plus(elt::Type{<:Number}, ::Type{Trivial}, ::Type{Trivial}; kwargs...)
    t = single_site_operator(elt, Trivial, Trivial)
    I = sectortype(t)
    t[(I(1), dual(I(1)))][1, 2] = 1.0
    return t
end
function S_plus(T::Type{<:Number}, ::Type{U1Irrep}, ::Type{Trivial}; filling::Rational{Int}=1//1)
    t = single_site_operator(T, U1Irrep, Trivial; filling=filling)
    P = numerator(filling); Q = denominator(filling)
    I = sectortype(t)
    t[(I(1, Q-P), dual(I(1, Q-P)))][1, 2] = 1.0
    return t
end
function S_plus(T::Type{<:Number}, ::Type{<:Sector}, ::Type{U1Irrep}; kwargs...)
    throw(ArgumentError("`S_plus`, `S_min` are not symmetric under `U1Irrep` spin symmetry"))
end
function S_plus(T::Type{<:Number}, ::Type{<:Sector}, ::Type{SU2Irrep}; kwargs...)
    throw(ArgumentError("`S_plus`, `S_min` are not symmetric under `SU2Irrep` spin symmetry"))
end

"""
    S_min(T::Type{<:Number}, P::Type{<:Sector}, S::Type{<:Sector})

Return the spin-lowering operator ``S^-``, the Hermitian conjugate of `S_plus`.
"""
S_min(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = S_min(ComplexF64, P, S; kwargs...)
function S_min(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return copy(adjoint(S_plus(T, particle_symmetry, spin_symmetry; kwargs...)))
end

"""
    Sx(T::Type{<:Number}, P::Type{<:Sector}, S::Type{<:Sector})

Return the spin operator ``S^x = (S^+ + S^-)/2``.
"""
Sx(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = Sx(ComplexF64, P, S; kwargs...)
function Sx(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return (S_plus(T, particle_symmetry, spin_symmetry; kwargs...) + S_min(T, particle_symmetry, spin_symmetry; kwargs...)) / 2
end

"""
    Sy(T::Type{<:Number}, P::Type{<:Sector}, S::Type{<:Sector})

Return the spin operator ``S^y = (S^+ - S^-)/(2i)``.
"""
Sy(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = Sy(ComplexF64, P, S; kwargs...)
function Sy(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return (S_plus(T, particle_symmetry, spin_symmetry; kwargs...) - S_min(T, particle_symmetry, spin_symmetry; kwargs...)) / (2im)
end

"""
    Sz(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the one-site spin operator ``S^z = (n_↑ - n_↓)/2``.
"""
Sz(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = Sz(ComplexF64, P, S; kwargs...)
function Sz(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    return 0.5 * (number_up(T, particle_symmetry, spin_symmetry; kwargs...) - number_down(T, particle_symmetry, spin_symmetry; kwargs...))
end

"""
    two_body(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the general two-body operator ``c†_{i} c†_{j} c_{k} c_{l}``.
"""
two_body(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = two_body(ComplexF64, P, S; kwargs...)
function two_body(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    c⁺c = c_plusmin(T, particle_symmetry, spin_symmetry; kwargs...)
    @tensor foursite[-1 -2 -3 -4; -5 -6 -7 -8] := c⁺c[-1 -4; -5 -8] * c⁺c[-2 -3; -6 -7]
    return foursite
end

"""
    three_body(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector})

Return the general three-body operator ``c†_{i} c†_{j} c†_{k} c_{l} c_{m} c_{n}``.
"""
three_body(P::Type{<:Sector}, S::Type{<:Sector}; kwargs...) = three_body(ComplexF64, P, S; kwargs...)
function three_body(T::Type{<:Number}, particle_symmetry::Type{<:Sector}, spin_symmetry::Type{<:Sector}; kwargs...)
    c⁺c = c_plusmin(T, particle_symmetry, spin_symmetry; kwargs...)
    @tensor sixsite[-1 -2 -3 -4 -5 -6; -7 -8 -9 -10 -11 -12] := c⁺c[-1 -6; -7 -12] * c⁺c[-2 -5; -8 -11] * c⁺c[-3 -4; -9 -10]
    return sixsite
end