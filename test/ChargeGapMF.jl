println("
###########################
#  Charge Gap Mean Field  #
###########################
")

@testset "ChargeGapMF construction" begin
    # Set up model parameters
    ps = U1Irrep
    ss = Trivial
    cell_width = 2

    symm = SymmetryConfig(ps, ss, cell_width)

    t = Dict((1, 2) => 1.0, (2, 1) => 1.0, (2, 3) => 1.0, (3, 2) => 1.0, (3, 4) => 1.0, (4, 3) => 1.0)
    U = Dict((1, 1, 1, 1) => 4.0, (2, 2, 2, 2) => 4.0, (3, 3, 3, 3) => 4.0)

    model = HubbardParams(3, t, U)

    # Mean field term
    t_inter = Dict((1, 1) => 0.5, (3, 4) => 0.5, (1, 3) => 0.5)
    range   = 1
    beta    = ones(6, 6)

    term = ChargeGapMF(t_inter, range, beta, beta)

    # Hamiltonian
    calc = CalcConfig(symm, model, term)
    @test typeof(hamiltonian(calc)) <: InfiniteMPOHamiltonian
end

@testset "ChargeGapMF solve" begin
    # Construct calc
    symm  = SymmetryConfig(U1Irrep, Trivial, 2, 1//1)
    model = HubbardParams([0.0, 1.0], [4.0])

    t_inter = Dict((1, 1) => 1.0)
    range   = 1
    beta    = zeros(2,2)
    term    = ChargeGapMF(t_inter, range, beta, [1.0 0.0; 0.0 1.0])

    calc = CalcConfig(symm, model, term)

    # Solve
    gs = compute_groundstate(calc; svalue=3.0)
    sx = Sx(symm.particle_symmetry, symm.spin_symmetry; filling=symm.filling)
    @test real(expectation_value(gs["groundstate"], (1,) => sx)) ≈ -0.5 atol=1e-6
    @test real(expectation_value(gs["groundstate"], (2,) => sx)) ≈ -0.5 atol=1e-6
end