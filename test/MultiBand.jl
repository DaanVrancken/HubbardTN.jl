println("
################
#  Multi-Band  #
################
")

tol = 1e-2


################
# Constructors #
################

@testset "ModelParams constructors" begin
    m1 = ModelParams([0.0 1.0; 1.0 0.0], Dict((1,1,1,1)=>5.0, (2,2,2,2)=>3.0))
    m2 = ModelParams([0.0 1.0; 1.0 0.0], [5.0 0.0; 0.0 3.0])
    m3 = ModelParams(2, Dict((1,2)=>1.0, (2,1)=>1.0, (1,1)=>0.0, (2,2)=>0.0), 
                     Dict((1,1,1,1)=>5.0, (2,2,2,2)=>3.0))
    @test m1.bands == m2.bands == m3.bands
    @test m1.t == m2.t == m3.t
    @test m1.U == m2.U == m3.U
end


#####################
# Hamiltonian terms #
#####################

symm = SymmetryConfig(U1Irrep, SU2Irrep, 2, (2,2))

t = [0.0 1.0 0.0 0.0; 1.0 0.0 0.5 0.0]
U =Dict((1,1,1,1)=>4.0,  # Direct on-site
        (2,2,2,2)=>4.0,
        (3,2,2,3)=>2.0,  # Direct inter-band
        (2,3,3,2)=>2.0,
        (1,2,1,2)=>1.0,  # Exchange
        (2,1,2,1)=>1.0,
        (1,1,2,2)=>1.0,  # Pair-exchange
        (2,2,1,1)=>1.0,
        (1,1,1,2)=>0.5,  # Bond-charge
        (1,1,2,1)=>0.5,
        (1,2,1,1)=>0.5,
        (2,1,1,1)=>0.5,
        (1,2,2,2)=>0.5,
        (2,1,2,2)=>0.5,
        (2,2,1,2)=>0.5,
        (2,2,2,1)=>0.5,
        (1,1,2,3)=>0.2,  # Three distinct sites
        (3,2,1,1)=>0.2,
        (1,2,3,4)=>0.1,  # Four distinct sites
        (4,3,2,1)=>0.1)

model = ModelParams(t , U)
calc = CalcConfig(symm, model)

E_norm = 0.38791

@testset "Hamiltonian terms" begin
    gs = compute_groundstate(calc; tol=tol/10)
    ψ₀ = gs["groundstate"]
    H = gs["ham"]
    E0 = expectation_value(ψ₀, H)/length(ψ₀)
    @test real(E0) ≈ E_norm atol=tol
end
