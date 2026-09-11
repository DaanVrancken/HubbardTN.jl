println("
########################
#  Chemical Potential  #
########################
")

tol = 1e-1

@testset "Find chemical potential" for spin_symmetry in [Trivial, U1Irrep, SU2Irrep]
    particle_symmetry = Trivial
    cell_width = 2

    symm = SymmetryConfig(particle_symmetry, spin_symmetry, cell_width)

    t = [0.0, 1.0] 
    U = [5.0]

    model = HubbardParams(t, U)
    calc = CalcConfig(symm, model)

    mu = find_chemical_potential(calc, 0.75; verbosity=0)
    @test mu ≈ 0.5364 atol=tol
end