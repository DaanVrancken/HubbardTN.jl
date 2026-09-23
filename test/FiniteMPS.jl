println("
################
#  Finite MPS  #
################
")

tol = 1e-2

s = 2.5
particle_symmetry = U1Irrep
spin_symmetry = U1Irrep

t = Dict((1,2)=>1.0, (2,1)=>1.0, (1,1)=>2.0)
U = Dict((1,1,1,1) => 4.0)
cell_width = 2
bands = 1

@testset "Compare infinite with finite MPS" begin
    symm = SymmetryConfig(particle_symmetry, spin_symmetry, cell_width)
    model = HubbardParams(bands, t, U)
    calc = CalcConfig(symm, model)
    gs = compute_groundstate(calc; svalue=s)
    ψ = gs["groundstate"]
    H = gs["ham"]

    E0 = expectation_value(ψ, H)
    Einf = sum(real(E0)) / length(H)

    symm = SymmetryConfig(particle_symmetry, spin_symmetry, cell_width*6)
    calc = CalcConfig(symm, model)

    gs = compute_groundstate(calc; svalue=s, finite_mps = true)
    ψ = gs["groundstate"]
    H = gs["ham"]

    E0 = expectation_value(ψ, H)
    E = sum(real(E0)) / length(H)
    @test Einf ≈ E atol=tol

    dim = dim_state(ψ)

    Ne = density_e(ψ, calc)
    @test sum(Ne) / length(Ne) ≈ 1.0 atol=tol

    ent = entanglement_spectrum(ψ, Int(bands*cell_width/2))
    println("Entanglement spectrum: \n")
    display(ent)
end

f = [1//2, 3//2];
E_norm = -1.768

@testset "Dependence on filling" for i in eachindex(f)
    symm = SymmetryConfig(U1Irrep, U1Irrep, cell_width*6, f[i])
    model = HubbardParams(bands, t, U)
    calc = CalcConfig(symm, model)
    gs = compute_groundstate(calc; svalue=s, finite_mps = true)
    ψ₀ = gs["groundstate"]
    H = gs["ham"]

    E0 = expectation_value(ψ₀, H)
    E = sum(real(E0))/length(ψ₀)
    @test E ≈ E_norm atol=tol
    
    Ne = density_e(ψ₀, calc)
    @test sum(Ne) / length(Ne) ≈ f[i] atol=tol
end