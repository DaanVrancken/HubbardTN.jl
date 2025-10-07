##################
# INITIALISATION #
##################

using HubbardTN
using MPSKit
using KrylovKit

# Extract name of the current file. Will be used as code name for the simulation.
name_jl = last(splitpath(Base.source_path()))
name = first(split(name_jl,"."))

# store calculations at
path = "sims"


#################
# DEFINE SYSTEM #
#################

s = 2.5             # Schmidt cut value, determines bond dimension.
P = 1;              # Filling of P/Q. P/Q = 1 is half-filling.
Q = 1;
bond_dim = 20;      # Initial bond dimension of the state. Impact on result is small as DMRG modifies it.

# Define hopping, direct and exchange interaction matrices.
t = [0.000 3.803 -0.548 0.000; 3.803 0.000 2.977 -0.501];
U = [10.317 6.264 0.000 0.000; 6.264 10.317 6.162 0.000];
J = [0.000 0.123 0.000 0.000; 0.123 0.000 0.113 0.000];

model = hf.MB_Sim(t, U, J, P, Q, s, bond_dim; code = name);


########################
# COMPUTE GROUNDSTATES #
########################

dictionary = hf.produce_groundstate(model; path=path);
ψ₀ = dictionary["groundstate"];
H = dictionary["ham"];
E0 = expectation_value(ψ₀, H);
E = sum(real(E0))./length(H);

println("Groundstate energy: $E")
println("Bond dimension: $(hf.dim_state(ψ₀))")


#######################
# COMPUTE EXCITATIONS #
#######################

resolution = 5;
momenta = range(0, π, resolution);
nums = 1;

exc = hf.produce_excitations(model, momenta, nums; charges=[0,0.0,0], path=path);
Es = exc["Es"];
println("Excitation energies: ")
println(Es)

println("Exciton energy for s=$s: $(real(Es[1,1]))")

gap, k = hf.produce_bandgap(model; path=path)

println("Band Gap for s=$s: $gap eV at momentum $k")
