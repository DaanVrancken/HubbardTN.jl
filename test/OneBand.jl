println("
##############
#  One-Band  #
##############
")

tol = 1e-2


############################
# DEPENDENCE ON PARAMETERS #
############################

u_range = 0:1:4
E = zeros(length(u_range),1)

E_norm = [-1.26967, -1.03717, -0.84145, -0.68707, -0.57119]

@testset "Dependence on parameters" for u in u_range
    symm = SymmetryConfig(U1Irrep, SU2Irrep, 2, (1,1))
    model = ModelParams([0.0, 1.0], [Float64(u)])
    calc = CalcConfig(symm, model)
    gs = compute_groundstate(calc; tol=tol/10)
    ψ₀ = gs["groundstate"]
    H = gs["ham"]
    E0 = expectation_value(ψ₀, H)
    E[u+1] = sum(real(E0))/length(ψ₀)
    @test E[u+1] ≈ E_norm[u+1] atol=tol
end


#########################
# DEPENDENCE ON FILLING #
#########################

t = [0.0, 1.0];
u = [5.0];
P = [1, 1, 3];
Q = [2, 1, 2];
E = zeros(length(P),1);

E_norm = [-0.73923, -0.48432, 1.76080]

@testset "Dependence on filling" for i in eachindex(P)
    symm = SymmetryConfig(U1Irrep, U1Irrep, 2*Q[i], (P[i],Q[i]))
    model = ModelParams(t, u)
    calc = CalcConfig(symm, model)
    gs = compute_groundstate(calc; tol=tol/10)
    ψ₀ = gs["groundstate"]
    H = gs["ham"]
    E0 = expectation_value(ψ₀, H)
    E[i] = sum(real(E0))/length(ψ₀)
    @test E[i] ≈ E_norm[i] atol=tol
end


###############
# EXCITATIONS #
###############

symm = SymmetryConfig(U1Irrep, Trivial, 2, (1,1))
model = ModelParams([0.0, 1.0], [5.0])
calc = CalcConfig(symm, model)

gs = compute_groundstate(calc; tol=tol/10)
psi = gs["groundstate"]

Es_norm = [4.80946; 4.71406; 4.43375; 3.95171; 3.60006]

@testset "Excitations" begin 
    nums = 1;
    resolution = 5;
    momenta = collect(range(0, π, resolution))

    exc = compute_excitations(gs, momenta, [1,1])
    Es = exc["Es"]
    @test real(Es)≈Es_norm atol=tol
    @test imag(Es)≈zeros(size(Es)) atol=1e-8
end


#########
# Tools #
#########

@testset "Tools" begin
    D = dim_state(psi)
    @test typeof(D) == Vector{Int64}
    @test D > zeros(size(D))

    electron_number = sum(density_e(psi, symm))/length(psi)
    @test sum(electron_number)≈symm.filling[1]/symm.filling[2] atol=1e-8

    Nup, Ndown = density_spin(psi, symm)
    @test sum(Nup+Ndown)/length(psi)≈symm.filling[1]/symm.filling[2] atol=1e-8

    ms = calc_ms(psi, symm)
    @test typeof(ms) == Float64
end