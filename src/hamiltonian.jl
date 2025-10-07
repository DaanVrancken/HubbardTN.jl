###############
# Hamiltonian #
###############

function SymSpace(P,Q,spin)
    if spin
        I = fℤ₂ ⊠ U1Irrep ⊠ U1Irrep
        Ps = Vect[I]((0, 0, -P) => 1, (0, 0, 2*Q-P) => 1, (1, 1, Q-P) => 1, (1, -1, Q-P) => 1)
    else
        I = fℤ₂ ⊠ SU2Irrep ⊠ U1Irrep
        Ps = Vect[I]((0, 0, -P) => 1, (0, 0, 2*Q-P) => 1, (1, 1 // 2, Q-P) => 1)
    end

    return I, Ps
end

function Hopping(P,Q,spin)
    I, Ps = SymSpace(P,Q,spin)

    if spin
        Vup = Vect[I]((1, 1, Q) => 1)
        Vdown = Vect[I]((1, -1, Q) => 1)
    
        c⁺u = zeros(ComplexF64, Ps ← Ps ⊗ Vup)
        blocks(c⁺u)[I((1, 1, Q-P))] .= 1
        blocks(c⁺u)[I((0, 0, 2*Q-P))] .= -1
        cu = zeros(ComplexF64, Vup ⊗ Ps ← Ps)
        blocks(cu)[I((1, 1, Q-P))] .= 1
        blocks(cu)[I((0, 0, 2*Q-P))] .= 1
        
        c⁺d = zeros(ComplexF64, Ps ← Ps ⊗ Vdown)
        blocks(c⁺d)[I((1, -1, Q-P))] .= 1
        blocks(c⁺d)[I((0, 0, 2*Q-P))] .= 1
        cd = zeros(ComplexF64, Vdown ⊗ Ps ← Ps)
        blocks(cd)[I((1, -1, Q-P))] .= 1
        blocks(cd)[I((0, 0, 2*Q-P))] .= -1
    
        @planar twosite_up[-1 -2; -3 -4] := c⁺u[-1; -3 1] * cu[1 -2; -4]
        @planar twosite_down[-1 -2; -3 -4] := c⁺d[-1; -3 1] * cd[1 -2; -4]
        twosite = twosite_up + twosite_down
    else
        Vs = Vect[I]((1, 1 / 2, Q) => 1)

        c⁺ = zeros(ComplexF64, Ps ← Ps ⊗ Vs)
        blocks(c⁺)[I((1, 1 // 2, Q-P))] .= 1
        blocks(c⁺)[I((0, 0, 2*Q-P))] .= sqrt(2)

        c = zeros(ComplexF64, Vs ⊗ Ps ← Ps)
        blocks(c)[I((1, 1 / 2, Q-P))] .= 1
        blocks(c)[I((0, 0, 2*Q-P))] .= sqrt(2)

        @planar twosite[-1 -2; -3 -4] := c⁺[-1; -3 1] * c[1 -2; -4]
    end

    return twosite
end

function OSInteraction(P,Q,spin)
    I, Ps = SymSpace(P,Q,spin)

    if spin
        onesite = zeros(ComplexF64, Ps ← Ps)
        blocks(onesite)[I((0, 0, 2*Q-P))] .= 1
    else
        onesite = zeros(ComplexF64, Ps ← Ps)
        blocks(onesite)[I((0, 0, 2*Q-P))] .= 1
    end

    return onesite
end

function Number(P,Q,spin)
    I, Ps = SymSpace(P,Q,spin)

    if spin
        n = zeros(ComplexF64, Ps ← Ps)
        blocks(n)[I((0, 0, 2*Q-P))] .= 2
        blocks(n)[I((1, 1, Q-P))] .= 1
        blocks(n)[I((1, -1, Q-P))] .= 1
    else
        n = zeros(ComplexF64, Ps ← Ps)
        blocks(n)[I((0, 0, 2*Q-P))] .= 2
        blocks(n)[I((1, 1 // 2, Q-P))] .= 1
    end

    return n
end

function Sz(P, Q)
    I, Ps = SymSpace(P, Q, true)
 
    sz = zeros(ComplexF64, Ps ← Ps)
 
    blocks(sz)[I((0, 0, 2*Q-P))] .= 0.0
    blocks(sz)[I((1, 1, Q-P))] .= 0.5
    blocks(sz)[I((1, -1, Q-P))] .= -0.5
 
    return sz
end

function SymSpace()
    I = fℤ₂ ⊠ SU2Irrep
    Ps = Vect[I]((0, 0) => 2, (1, 1 // 2) => 1)

    return I, Ps
end

function Hopping()
    I, Ps = SymSpace()
    Vs = Vect[I]((1, 1 / 2) => 1)

    c⁺ = zeros(ComplexF64, Ps ← Ps ⊗ Vs)
    blocks(c⁺)[I((1, 1 // 2))][1] = 1.0+0.0im
    blocks(c⁺)[I((0, 0))][2] = sqrt(2)+0.0im

    c = zeros(ComplexF64, Vs ⊗ Ps ← Ps)
    blocks(c)[I((1, 1 // 2))][1] = 1.0+0.0im
    blocks(c)[I((0, 0))][2] = sqrt(2)+0.0im

    @planar twosite[-1 -2; -3 -4] := c⁺[-1; -3 1] * c[1 -2; -4]
    
    return twosite
end

function OSInteraction()
    I, Ps = SymSpace()

    onesite = zeros(ComplexF64, Ps ← Ps)
    blocks(onesite)[I((0, 0))][2,2] = 1.0

    return onesite
end

function Number()
    I, Ps = SymSpace()

    n = zeros(ComplexF64, Ps ← Ps)
    blocks(n)[I((0, 0))][2,2] = 2.0 
    blocks(n)[I((1, 1 // 2))] .= 1.0

    return n
end

# ONEBAND #

function hamiltonian(simul::Union{OB_Sim,OBC_Sim2})
    t = simul.t
    u = simul.u
    μ = simul.μ
    if hasproperty(simul, :J)
        J = simul.J
        D_exc = length(J)
    end
    L = simul.period
    spin::Bool = get(simul.kwargs, :spin, false)
    U13::Vector{Float64} = get(simul.kwargs, :U13, [0.0])
    JMs::Tuple{Float64, Float64} = get(simul.kwargs, :JMs, (0.0,0.0))
    J_inter = JMs[1]
    Ms = JMs[2]

    D_hop = length(t)
    D_int = length(u)
    D_U13 = length(U13)
    
    if hasproperty(simul, :P)
        P = simul.P
        Q = simul.Q
        if iseven(P)
            T = Q
        else 
            T = 2*Q
        end
        cdc = Hopping(P,Q,spin)    
        n = Number(P,Q,spin)
        OSI = OSInteraction(P,Q,spin)
    else
        T = 1
        cdc = Hopping()
        n = Number()
        OSI = OSInteraction()
    end

    twosite = cdc + cdc'
    onesite = u[1]*OSI - μ*n

    @planar nn[-1 -2; -3 -4] := n[-1; -3] * n[-2; -4]
    @tensor J1[-1 -2; -3 -4] := cdc[-1 2; 3 -4] * cdc[-2 3; 2 -3]
    @tensor J2[-1 -2; -3 -4] := cdc[-1 2; 3 -4] * cdc[3 -2; -3 2]
    @tensor C1[-1 -2; -3 -4] := cdc[-1 2; -3 -4] * cdc[-2 3; 3 2]
    @tensor C2[-1 -2; -3 -4] := cdc[-1 2; -3 4] * cdc[-2 4; 2 -4]

    C1 = C1 + C1'
    C2 = C2 + C2'
    
    H = @mpoham sum(onesite{i} for i in vertices(InfiniteChain(T)))
    if L == 0
        for range_hop in 1:D_hop
            h = @mpoham sum(-t[range_hop]*twosite{i,i+range_hop} for i in vertices(InfiniteChain(T)))
            H += h
        end
        for range_int in 2:D_int  # first element is on-site interaction
            h = @mpoham sum(u[range_int]*nn{i,i+(range_int-1)} for i in vertices(InfiniteChain(T)))
            H += h
        end
        if hasproperty(simul, :J)
            for range_exc in 1:D_exc
                h1 = @mpoham sum(J[range_exc]*J1{i,i+range_exc} for i in vertices(InfiniteChain(T)))
                h2 = @mpoham sum(0.5*J[range_exc]*J2{i,i+range_exc} + 0.5*J[range_exc]*J2{i+range_exc,i} for i in vertices(InfiniteChain(T)))
                H += h1 + h2
            end
        end
        if U13 != [0.0]
            for range_U13 in 1:D_U13
                h3 = @mpoham sum(0.5*U13[range_U13]*C1{i,i+range_U13} + 0.5*U13[range_U13]*C2{i,i+range_U13} for i in vertices(InfiniteChain(T)))
                h4 = @mpoham sum(0.5*U13[range_U13]*C1{i+range_U13,i} + 0.5*U13[range_U13]*C2{i+range_U13,i} for i in vertices(InfiniteChain(T)))
                H += h3 + h4
            end
        end
        if Ms!=0.0 && spin
            sz = Sz(P,Q)
            h = @mpoham sum(J_inter*Ms*sz{vertex}*(-1)^i for (i, vertex) in enumerate(vertices(InfiniteChain(T))))
            H += h
        end
    elseif D_hop==1 && D_int==1
        h = @mpoham sum(-t[1]*twosite{i,i+1} -t[1]*twosite{i,i+L} for i in vertices(InfiniteChain(T)))
        H += h
    else
        return error("Extended models in 2D not implemented.")
    end

    return H
end

# MULTIBAND #

# t[i,j] gives the hopping of band i on one site to band j on the same site (i≠j)
function OS_Hopping(t,T,cdc)
    Bands,Bands2 = size(t)
    
    if Bands ≠ Bands2 || typeof(t) ≠ Matrix{Float64}
        @warn "t_OS is not a float square matrix."
    end
    for i in 1:Bands
        for j in (i+1):Bands
            if !(t[i,j] ≈ t'[i,j])
                @warn "t_OS is not Hermitian."
            end
        end
    end
    
    Lattice = InfiniteStrip(Bands,T*Bands)
        
    # Define necessary different indices of sites/orbitals in the lattice
    # Diagonal terms are taken care of in chem_pot
    Indices = [(div(l-1,Bands^2)+1, div((l-1)%(Bands^2),Bands)+1, mod(l-1,Bands)+1) 
               for l in 1:(T*Bands^2) if div((l-1)%(Bands^2),Bands)+1 ≠ mod(l-1,Bands)+1]
    
    return @mpoham sum(-t[bi,bf]*cdc{Lattice[bf,site],Lattice[bi,site]} for (site, bi, bf) in Indices)
end

# t[i,j] gives the hopping of band i on one site to band j on the range^th next site
# parameter must be equal in both directions (1i->2j=2j->1i) to guarantee hermiticity
function IS_Hopping(t,range,T,cdc)
    Bands,Bands2 = size(t)
    if Bands ≠ Bands2 || typeof(t) ≠ Matrix{Float64}
        @warn "t_IS is not a float square matrix"
    end
    
    twosite = cdc + cdc'
    Lattice = InfiniteStrip(Bands,T*Bands)
        
    # Define necessary different indices of sites/orbitals in the lattice
    Indices = [(div(l-1,Bands^2)+1, div((l-1)%(Bands^2),Bands)+1, mod(l-1,Bands)+1) for l in 1:(T*Bands^2)]
    
    return @mpoham sum(-t[bi,bf]*twosite{Lattice[bf,site+range],Lattice[bi,site]} for (site, bi, bf) in Indices)
end

# μ[i] gives the hopping of band i on one site to band i on the same site.
function Chem_pot(μ,T,n)
    Bands = length(μ)
    
    Lattice = InfiniteStrip(Bands,T*Bands)
    
    Indices = [(div(l-1,Bands)+1, mod(l-1,Bands)+1) for l in 1:(T*Bands)]
    
    return @mpoham sum(-μ[i]*n{Lattice[i,j]} for (j,i) in Indices)
end

# u[i] gives the interaction on band i
function OB_interaction(u,T,OSI)
    Bands = length(u)
    
    Lattice = InfiniteStrip(Bands,T*Bands)
    
    Indices = [(div(l-1,Bands)+1, mod(l-1,Bands)+1) for l in 1:(T*Bands)]
    
    return @mpoham sum(u[i]*OSI{Lattice[i,j]} for (j,i) in Indices)
end

# U[i,j] gives the direct interaction between band i on one site to band j on the same site. Averaged over U[i,j] and U[j,i]
function Direct_OS(U,T,n)
    Bands,Bands2 = size(U)
    
    if Bands ≠ Bands2 || typeof(U) ≠ Matrix{Float64}
        @warn "U_OS is not a float square matrix"
    end
    
    U_av = zeros(Bands,Bands2)
    for i in 2:Bands    
        for j in 1:(i-1)
            U_av[i,j] = 0.5*(U[i,j]+U[j,i])
        end
    end
    
    @planar nn[-1 -2; -3 -4] := n[-1; -3] * n[-2; -4]
    Lattice = InfiniteStrip(Bands,T*Bands)
    
    Indices = [(div(l-1,Bands^2)+1, div((l-1)%(Bands^2),Bands)+1, mod(l-1,Bands)+1) 
               for l in 1:(T*Bands^2) if div((l-1)%(Bands^2),Bands)+1 ≠ mod(l-1,Bands)+1]
    
    return @mpoham sum(U_av[bi,bf]*nn{Lattice[bi,site],Lattice[bf,site]} for (site,bi,bf) in Indices if U_av[bi,bf]≠0.0)
end

# J[i,j] gives the exchange interaction between band i on one site to band j on the same site.
function Exchange1_OS(J,T,cdc)
    Bands,Bands2 = size(J)
    
    if Bands ≠ Bands2 || typeof(J) ≠ Matrix{Float64}
        @warn "J_OS is not a float square matrix"
    end
    diagonal = zeros(Bands,1)
    diagonal_zeros = zeros(Bands,1)
    for i in 1:Bands
        diagonal[i] = J[i,i]
    end
    if diagonal≠diagonal_zeros
        @warn "On-band interaction is not taken into account in Exchange_OS."
    end
    
    @tensor C4[-1 -2; -3 -4] := cdc[-1 2; 3 -4] * cdc[-2 3; 2 -3]
    Lattice = InfiniteStrip(Bands,T*Bands)
    
    Indices = [(div(l-1,Bands^2)+1, div((l-1)%(Bands^2),Bands)+1, mod(l-1,Bands)+1) 
               for l in 1:(T*Bands^2) if div((l-1)%(Bands^2),Bands)+1 ≠ mod(l-1,Bands)+1]
    
    return @mpoham sum(0.5*J[bi,bf]*C4{Lattice[bi,site],Lattice[bf,site]} for (site,bi,bf) in Indices)
end;

function Exchange2_OS(J,T,cdc)
    Bands,Bands2 = size(J)
    
    if Bands ≠ Bands2 || typeof(J) ≠ Matrix{Float64}
        @warn "J_OS is not a float square matrix"
    end
    diagonal = zeros(Bands,1)
    diagonal_zeros = zeros(Bands,1)
    for i in 1:Bands
        diagonal[i] = J[i,i]
    end
    if diagonal≠diagonal_zeros
        @warn "On-band interaction is not taken into account in Exchange_OS."
    end
    
    @tensor C4[-1 -2; -3 -4] := cdc[-1 2; 3 -4] * cdc[3 -2; -3 2]
    Lattice = InfiniteStrip(Bands,T*Bands)
    
    Indices = [(div(l-1,Bands^2)+1, div((l-1)%(Bands^2),Bands)+1, mod(l-1,Bands)+1) 
               for l in 1:(T*Bands^2) if div((l-1)%(Bands^2),Bands)+1 ≠ mod(l-1,Bands)+1]
    
    return @mpoham sum(0.5*J[bi,bf]*C4{Lattice[bi,site],Lattice[bf,site]} for (site,bi,bf) in Indices)
end;

function Exchange_OS(J,T,cdc)
    return Exchange1_OS(J,T,cdc) + Exchange2_OS(J,T,cdc)
end;

function Uijjj_OS(U,T,cdc)
    Bands,Bands2 = size(U)
    
    if Bands ≠ Bands2 || typeof(U) ≠ Matrix{Float64}
        @warn "U13_OS is not a float square matrix"
    end
    diagonal = zeros(Bands,1)
    diagonal_zeros = zeros(Bands,1)
    for i in 1:Bands
        diagonal[i] = U[i,i]
    end
    if diagonal≠diagonal_zeros
        @warn "On-band interaction is not taken into account in Exchange_OS."
    end
    
    @tensor C1[-1 -2; -3 -4] := cdc[-1 2; -3 -4] * cdc[-2 3; 3 2]
    #@tensor C2[-1 -2; -3 -4] := cdc[-2 -1; 3 -3] * cdc[3 2; 2 -4]
    @tensor C2[-1 -2; -3 -4] := cdc[-1 2; -3 4] * cdc[-2 4; 2 -4]
    #@tensor C4[-1 -2; -3 -4] := cdc[1 -1; 3 -3] * cdc[-2 3; 1 -4]

    C1 = C1 + C1'
    C2 = C2 + C2'

    Lattice = InfiniteStrip(Bands,T*Bands)
    
    Indices = [(div(l-1,Bands^2)+1, div((l-1)%(Bands^2),Bands)+1, mod(l-1,Bands)+1) 
               for l in 1:(T*Bands^2) if div((l-1)%(Bands^2),Bands)+1 ≠ mod(l-1,Bands)+1]
    
    H = @mpoham sum(0.5*U[bi,bf]*C1{Lattice[bi,site],Lattice[bf,site]} for (site,bi,bf) in Indices)
    H += @mpoham sum(0.5*U[bi,bf]*C2{Lattice[bi,site],Lattice[bf,site]} for (site,bi,bf) in Indices)

    return H
end;

# V[i,j] gives the direct interaction between band i on one site to band j on the range^th next site.
function Direct_IS(V,range,T,n)
    Bands,Bands2 = size(V)
    
    if Bands ≠ Bands2 || typeof(V) ≠ Matrix{Float64}
        @warn "V is not a float square matrix"
    end
    
    @planar nn[-1 -2; -3 -4] := n[-1; -3] * n[-2; -4]
    Lattice = InfiniteStrip(Bands,T*Bands)
    
    Indices = [(div(l-1,Bands^2)+1, div((l-1)%(Bands^2),Bands)+1, mod(l-1,Bands)+1) for l in 1:(T*Bands^2)]
    
    return @mpoham sum(V[bi,bf]*nn{Lattice[bi,site],Lattice[bf,site+range]} for (site,bi,bf) in Indices)
end

# J[i,j] gives the exchange interaction between band i on one site to band j on the range^th next site.
function Exchange1_IS(J,range,T,cdc)
    Bands,Bands2 = size(J)
    
    if Bands ≠ Bands2 || typeof(J) ≠ Matrix{Float64}
        @warn "J_IS is not a float square matrix"
    end
    
    @tensor C4[-1 -2; -3 -4] := cdc[-1 2; 3 -4] * cdc[-2 3; 2 -3]
    Lattice = InfiniteStrip(Bands,T*Bands)
    
    Indices = [(div(l-1,Bands^2)+1, div((l-1)%(Bands^2),Bands)+1, mod(l-1,Bands)+1) for l in 1:(T*Bands^2)]
    
    return @mpoham sum(J[bi,bf]*C4{Lattice[bi,site],Lattice[bf,site+range]} for (site,bi,bf) in Indices)    # operator has no direction
end;

function Exchange2_IS(J,range,T,cdc)
    Bands,Bands2 = size(J)
    
    if Bands ≠ Bands2 || typeof(J) ≠ Matrix{Float64}
        @warn "J_IS is not a float square matrix"
    end
    
    @tensor C4[-1 -2; -3 -4] := cdc[-1 2; 3 -4] * cdc[3 -2; -3 2]
    Lattice = InfiniteStrip(Bands,T*Bands)
    
    Indices = [(div(l-1,Bands^2)+1, div((l-1)%(Bands^2),Bands)+1, mod(l-1,Bands)+1) for l in 1:(T*Bands^2)]
    
    return @mpoham sum(0.5*J[bi,bf]*C4{Lattice[bi,site],Lattice[bf,site+range]} + 0.5*J[bi,bf]*C4{Lattice[bf,site+range],Lattice[bi,site]} for (site,bi,bf) in Indices) #operator has direction
end;

function Exchange_IS(J,range,T,cdc)
    return Exchange1_IS(J,range,T,cdc) + Exchange2_IS(J,range,T,cdc)
end;

# Four different matrices required: two for U13 and two for U31
function Uijjj_IS(U,range,T,cdc)
    Bands,Bands2,num = size(U)
    
    if Bands ≠ Bands2
        @warn "U13_IS is not a float square matrix"
    elseif num != 4
        # i = orbital 1 on site 0, j = orbital 2 on site "range"
        # index 1: Uijjj=Ujjji, index 2: Ujiii=Uiiij, index 3: Ujijj=Ujjij, index 4: Uijii=Uiiji
        error("U13_IS shoud be a BxBx4 array.")
    end
    
    @tensor C1[-1 -2; -3 -4] := cdc[-1 2; -3 -4] * cdc[-2 3; 3 2]
    #@tensor C2[-1 -2; -3 -4] := cdc[-2 -1; 3 -3] * cdc[3 2; 2 -4]
    @tensor C2[-1 -2; -3 -4] := cdc[-1 2; -3 4] * cdc[-2 4; 2 -4]
    #@tensor C4[-1 -2; -3 -4] := cdc[1 -1; 3 -3] * cdc[-2 3; 1 -4]

    C1 = C1 + C1'
    C2 = C2 + C2'

    Lattice = InfiniteStrip(Bands,T*Bands)
    
    Indices = [(div(l-1,Bands^2)+1, div((l-1)%(Bands^2),Bands)+1, mod(l-1,Bands)+1) for l in 1:(T*Bands^2)]
    
    H = @mpoham sum(0.5*U[bi,bf,1]*C1{Lattice[bi,site],Lattice[bf,site+range]} + 0.5*U[bi,bf,3]*C1{Lattice[bf,site+range],Lattice[bi,site]} for (site,bi,bf) in Indices) #operator has direction
    H += @mpoham sum(0.5*U[bi,bf,2]*C2{Lattice[bi,site],Lattice[bf,site+range]} + 0.5*U[bi,bf,4]*C2{Lattice[bf,site+range],Lattice[bi,site]} for (site,bi,bf) in Indices)

    return H
end;

function Uijkk(U::Dict{NTuple{4, Int64}, Float64},B,T,cdc)
    # input is dict with permutations i,j,k,l (in order Cdi Cdj Ck Cl, NOT Uijkl). Indices range over i,j,k,l = 1,...,r*B
    # At least one index in every tuple (i,j,k,l) has to be at site 0

    Ind1 = []
    Ind2 = []
    Ind3 = []
    for (i,j,k,l) in keys(U)
        if minimum((i,j,k,l)) > B
            error("At least one index in every tuple (i,j,k,l) has to be at site 0.")
        elseif length(unique((i, j, k, l))) != 3
            error("Two indices should be the same. Not more, not less.")
        end
        for site in 1:T
            if k==l
                push!(Ind1,(site,i,j,k,l))
            elseif j==k
                push!(Ind2,(site,i,j,k,l))
            elseif j==l
                push!(Ind3,(site,i,j,k,l))
            end
        end
    end

    @tensor C1[-1 -2 -3; -4 -5 -6] := cdc[-1 2; -4 -6] * cdc[-2 -3; -5 2]
    @tensor C2[-1 -2 -3; -4 -5 -6] := cdc[-1 -3; -4 -6] * cdc[-2 2; 2 -5]
    @tensor C3[-1 -2 -3; -4 -5 -6] := cdc[-1 2; -4 -5] * cdc[-2 -3; 2 -6]
    C1 = C1 + C1'
    C2 = C2 + C2'
    C3 = C3 + C3'

    Lattice = InfiniteStrip(B,T*B)
    
    @tensor init_operator[-1; -2] := cdc[-1 2; 2 -2]
    H = @mpoham sum(0.0*init_operator{i} for i in vertices(Lattice))
    if !isempty(Ind1)
        H += @mpoham sum(0.5*U[(i,j,k,l)]*C1{Lattice[mod(i-1,B)+1,site+(i-1)÷B],Lattice[mod(j-1,B)+1,site+(j-1)÷B],Lattice[mod(k-1,B)+1,site+(k-1)÷B]} for (site,i,j,k,l) in Ind1)
    end
    if !isempty(Ind2)
        H += @mpoham sum(U[(i,j,k,l)]*C2{Lattice[mod(i-1,B)+1,site+(i-1)÷B],Lattice[mod(j-1,B)+1,site+(j-1)÷B],Lattice[mod(l-1,B)+1,site+(l-1)÷B]} for (site,i,j,k,l) in Ind2)
    end
    if !isempty(Ind3)
        H += @mpoham sum(0.5*U[(i,j,k,l)]*C3{Lattice[mod(i-1,B)+1,site+(i-1)÷B],Lattice[mod(j-1,B)+1,site+(j-1)÷B],Lattice[mod(k-1,B)+1,site+(k-1)÷B]} for (site,i,j,k,l) in Ind3)
    end

    return H
end;

function Uijkl(U::Dict{NTuple{4, Int64}, Float64},B,T,cdc)
    # input is dict with permutations i,j,k,l (in order Cdi Cdj Ck Cl, NOT Uijkl). Indices range over i,j,k,l = 1,...,r*B
    # At least one index in every tuple (i,j,k,l) has to be at site 0

    Ind = []
    for (i,j,k,l) in keys(U)
        if minimum((i,j,k,l)) > B
            error("At least one index in every tuple (i,j,k,l) has to be at site 0.")
        elseif length(unique((i, j, k, l))) != 4
            error("All indices must be different.")
        elseif !(U[(i,j,k,l)] ≈ U[(l,k,j,i)])
            @warn("U1111 is not Hermitian.")
        end
        for site in 1:T
            push!(Ind,(site,i,j,k,l))
        end
    end

    @tensor C[-1 -2 -3 -4; -5 -6 -7 -8] := cdc[-1 -2; -5 -6] * cdc[-3 -4; -7 -8]

    Lattice = InfiniteStrip(B,T*B)

    @tensor init_operator[-1; -2] := cdc[-1 2; 2 -2]
    H = @mpoham sum(0.0*init_operator{i} for i in vertices(Lattice))
    if !isempty(Ind)
        H += @mpoham sum(0.5*U[(i,j,k,l)]*C{Lattice[mod(i-1,B)+1,site+(i-1)÷B],Lattice[mod(l-1,B)+1,site+(l-1)÷B],Lattice[mod(j-1,B)+1,site+(j-1)÷B,],Lattice[mod(k-1,B)+1,site+(k-1)÷B]} for (site,i,j,k,l) in Ind)
    end

    return H
end

function hamiltonian(simul::Union{MB_Sim, MBC_Sim})
    t = simul.t
    u = simul.u
    J = simul.J
    U13_OS = simul.U13
    U112::Dict{NTuple{4, Int64}, Float64} = get(simul.kwargs, :U112, Dict{Tuple{Int, Int, Int, Int}, Float64}())
    U1111::Dict{NTuple{4, Int64}, Float64} = get(simul.kwargs, :U1111, Dict{Tuple{Int, Int, Int, Int}, Float64}())
    spin::Bool = get(simul.kwargs, :spin, false)

    Bands,width_t = size(t)
    Bands1,width_u = size(u)
    Bands2,width_J = size(J)
    Bands3,_ = size(U13_OS)
    U13_IS::Array{Float64, 3} = get(simul.kwargs, :U13_IS, zeros(Bands,Bands,4))
    if !(Bands == Bands1 == Bands2 == Bands3 == size(U13_IS)[1])
        return error("Number of bands is incosistent.")
    end

    if hasproperty(simul, :P)
        P = simul.P
        Q = simul.Q
        if iseven(P)
            T = Q
        else 
            T = 2*Q
        end
        cdc = Hopping(P,Q,spin)
        OSI = OSInteraction(P,Q,spin)
        n = Number(P,Q,spin)
    else
        T = 1
        cdc = Hopping()
        OSI = OSInteraction()
        n = Number()
    end

    Range_t = Int((width_t-Bands)/Bands)
    Range_u = Int((width_u-Bands)/Bands)
    Range_J = Int((width_J-Bands)/Bands)
    Range_U13 = Int((size(U13_IS)[2])/Bands)

    # Define matrices
    u_OB = zeros(Bands)
    for i in 1:Bands
        u_OB[i] = u[i,i]
    end
    if u_OB == zeros(Bands)
        @warn "No on-band interaction found. This may lead to too low contributions of other Hamiltonian terms."
    end
    t_OS = t[:,1:Bands]
    μ = zeros(Bands)
    for i in 1:Bands
        μ[i] = t_OS[i,i]
    end
    u_OS = u[:,1:Bands]
    for i in 1:Bands
        u_OS[i,i] = 0.0
    end
    J_OS = J[:,1:Bands]

    # Implement Hamiltonian OB
    H_total = OB_interaction(u_OB,T,OSI)

    if μ != zeros(Bands)
        H_total += Chem_pot(μ,T,n)
    end

    # Implement Hamiltonian OS
    for (m,o,f) in [(t_OS,cdc,OS_Hopping),(u_OS,n,Direct_OS),(J_OS,cdc,Exchange_OS),(U13_OS,cdc,Uijjj_OS)]
        if m != zeros(Bands,Bands)
            H_total += f(m,T,o)
        end
    end

    # Implement Hamiltonian IS
    for (m,range,o,f) in [(t,Range_t,cdc,IS_Hopping),(u,Range_u,n,Direct_IS),(J,Range_J,cdc,Exchange_IS),(U13_IS,Range_U13,cdc,Uijjj_IS)]
        for i in 1:range
            if m != U13_IS
                M = m[:,(Bands*i+1):(Bands*(i+1))]
                ZERO = zeros(Bands,Bands)
            else
                M = m[:,(Bands*(i-1)+1):(Bands*i),:]
                ZERO = zeros(Bands,Bands,4)
            end
            if M != ZERO
                H_total += f(M,i,T,o)
            end
        end
    end

    if !isempty(U112)
        H_total += Uijkk(U112::Dict{NTuple{4, Int64}, Float64},Bands,T,cdc)
    end

    if !isempty(U1111)
        H_total += Uijkl(U1111::Dict{NTuple{4, Int64}, Float64},Bands,T,cdc)
    end

    return H_total
end