println("
#################
#  Excitations  #
#################
")

tol = 1e-2

# Parameters
t = [0.0, 1.0];
u = [5.0];
filling = 2//3;


##########
# U1xSU2 #
##########

@testset "U1xSU2" begin
    symm = SymmetryConfig(U1Irrep, SU2Irrep, denominator(filling), filling)
    model = HubbardParams(t, u)
    calc = CalcConfig(symm, model)
    gs = compute_groundstate(calc; tol=tol/10)

    gap, kmin = compute_bandgap(gs, calc; resolution=5)
    @test gap >= 0.0
    @test 0.0 <= kmin <= π

    gap, kmin = compute_spingap(gs, calc; resolution=5)
    @test gap >= 0.0
    @test 0.0 <= kmin <= π

    gap, kmin = compute_pairing_energy(gs, calc; resolution=5)
    @test isreal(gap)
    @test 0.0 <= kmin <= π
end


###########
# U1xTriv #
###########

@testset "U1xTriv" begin
    symm = SymmetryConfig(U1Irrep, Trivial, denominator(filling), filling)
    model = HubbardParams(t, u)
    calc = CalcConfig(symm, model)
    gs = compute_groundstate(calc; tol=tol/10)

    gap, kmin = compute_bandgap(gs, calc; resolution=5)
    @test gap >= 0.0
    @test 0.0 <= kmin <= π

    gap, kmin = compute_pairing_energy(gs, calc; resolution=5)
    @test isreal(gap)
    @test 0.0 <= kmin <= π
end


###########
# TrivxU1 #
###########

t = [first(u)/2, 1.0];

@testset "TrivxU1" begin
    symm = SymmetryConfig(Trivial, U1Irrep, denominator(filling))
    model = HubbardParams(t, u)
    calc = CalcConfig(symm, model)
    gs = compute_groundstate(calc; tol=tol/10)

    gap, kmin = compute_spingap(gs, calc; resolution=5)
    @test gap >= 0.0
    @test 0.0 <= kmin <= π
end

#############
# FiniteMPS #
#############

# @testset "FiniteMPS" begin
#     cw = 6
#     symm = SymmetryConfig(U1Irrep, U1Irrep, cw, 1//1)
#     model = HubbardParams(t, u)
#     calc = CalcConfig(symm, model)
#     gs = compute_groundstate(calc; tol=tol/10, finite_mps=true)

#     gap, kmin = compute_spingap(gs, calc; resolution=5)
#     @test gap >= 0.0
#     @test 0.0 <= kmin <= π
# end