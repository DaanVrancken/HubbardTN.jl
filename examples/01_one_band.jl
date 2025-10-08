##################
# INITIALISATION #
##################

using HubbardTN
using MPSKit
using KrylovKit

# store calculations at
path = joinpath("data", "one_band")


#################
# DEFINE SYSTEM #
#################

s = 2.5             # Schmidt cut value, determines bond dimension.
P = 1;              # Filling of P/Q. P/Q = 1 is half-filling.
Q = 1;
bond_dim = 20;      # Initial bond dimension of the state. Impact on result is small as DMRG modifies it.

# Define hopping, direct interaction, and chemical potential.
t=[1.0, 0.1];
u=[8.0];
J=[0.0];
μ=0.0;

# Spin=false will use SU(2) spin symmetry, the exact spin configuration cannot be deduced.
Spin = false

model = OB_Sim(t, u, μ, J, P, Q, s, bond_dim; spin=Spin);


########################
# COMPUTE GROUNDSTATES #
########################

dictionary = produce_groundstate(model; force=false, path=path);
ψ₀ = dictionary["groundstate"];
H = dictionary["ham"];
E0 = expectation_value(ψ₀, H);
E = sum(real(E0))/length(H);

println("Groundstate energy: $E")
println("Bond dimension: $(dim_state(ψ₀))")


########################
# COMPUTE EXCITATIONS #
########################

resolution = 5;
momenta = range(0, π, resolution);
nums = 1;

exc = produce_excitations(model, momenta, nums; charges=[0,0.0,0], path=path);
Es = exc["Es"];
println("Excitation energies: ")
println(Es)

println("Exciton energy for s=$s: $(real(Es[1,1]))")

gap, k = produce_bandgap(model: path=path)

println("Band Gap for s=$s: $gap eV at momentum $k")