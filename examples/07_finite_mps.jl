using HubbardTN
using TensorKit, MPSKit

s = 2.0
particle_symmetry = U1Irrep
spin_symmetry = U1Irrep
bands = 1
cell_width = 8

t = Dict((1,2)=>1.0, (2,1)=>1.0, (1,1)=>2.0)

U = Dict(
    (1,1,1,1) => 4.0)


model = HubbardParams(bands, t, U)
symm = SymmetryConfig(particle_symmetry, spin_symmetry, cell_width) #Length will be bands * cell_width
calc = CalcConfig(symm, model)

gs = compute_groundstate(calc; svalue=s, finite_mps = true)
ψ = gs["groundstate"]
H = gs["ham"]

E0 = expectation_value(ψ, H)
E = sum(real(E0)) / length(H)
println("Groundstate energy: ", E)
dim = dim_state(ψ)
println("Max bond dimension: ", maximum(dim))

Ne = density_e(ψ, calc)
println("Number of electrons per site: ", sum(Ne) / length(Ne))

ent = entanglement_spectrum(ψ, Int(bands*cell_width/2))
println("Entanglement spectrum: \n")
display(ent)