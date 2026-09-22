println("
################
#  Impurity MPS  #
################
")

tol = 1e-2

s = 2.5
filling = 1//1
cell_width = 12
bands = 1

particle_symmetry = U1Irrep
spin_symmetry = U1Irrep

t = Dict((1,2)=>1.0, (2,1)=>1.0, (1,1)=>2.0)
U = Dict((1,1,1,1) => 4.0)

# The indices below specify the site(s) where the impurity is located.
# t_imp and U_imp represent the difference between the impurity and bulk parameters.
t_imp = Dict{NTuple{2, Int64}, Float64}()
U_imp = Dict{NTuple{4, Int64}, Float64}((6,6,6,6) => 2.0)

symm = SymmetryConfig(particle_symmetry, spin_symmetry, cell_width, filling)
model = HubbardParams(bands, t, U)
calc = CalcConfig(symm, model, ImpurityTerm(t_imp, U_imp))
gs = compute_groundstate(calc; svalue=s, finite_mps=true)
ψ = gs["groundstate"]
H = gs["ham"]

E0 = expectation_value(ψ, H)
Ne = density_e(ψ, calc)
@test sum(Ne) / length(Ne) ≈ 1.0 atol=tol
@test maximum(Ne) - minimum(Ne) > tol