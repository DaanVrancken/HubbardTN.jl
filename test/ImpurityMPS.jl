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

t_tag = dict_tag(t)
U_tag = dict_tag(U)

t_imp = Dict((1,2)=>1.0, (2,1)=>1.0, (1,1)=>2.0)
U_imp = Dict((1,1,1,1) => 4.0, (1,2,2,1) => 2.0, (2,1,1,2) => 2.0)

t_tag_imp = dict_tag(t_imp)
U_tag_imp = dict_tag(U_imp)

symm = SymmetryConfig(particle_symmetry, spin_symmetry, cell_width, filling)
model = HubbardParams(bands, t, U, t_imp, U_imp)
calc = CalcConfig(symm, model)
gs = compute_groundstate(calc; svalue=s, finite_mps=true, imp_mps=true)
ψ = gs["groundstate"]
H = gs["ham"]

E0 = expectation_value(ψ, H)
Ne = density_e(ψ, calc)
@test sum(Ne) / length(Ne) ≈ 1.0 atol=tol
@test maximum(Ne) - minimum(Ne) > tol