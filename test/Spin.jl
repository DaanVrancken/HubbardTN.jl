println("
##########
#  Spin  #
##########
")

##################
# INITIALISATION #
##################

Force = true
tol = 1e-1

# store calculations at
path = joinpath("sims","Spin")

model1 = hf.OB_Sim([1.0],[8.0], 0.0,1,1,2.0;spin=true)

# Extract name of the current file. Will be used as code name for the simulation.
name_jl = last(splitpath(Base.source_path()))
name = first(split(name_jl,"."))

t_OS = [0.0 0.0; 0.0 0.0];
t_IS = [1.0 0.0; 0.0 1.0];
t = cat(t_OS,t_IS, dims=2)
U = [3.0 0.0; 0.0 3.0]
V = [0.0 0.0; 0.0 0.0]
u = cat(U,V,dims=2)
J = [0.0 0.0; 0.0 0.0]
P = 1;
Q = 1;
bond_dim = 20;
model2 = hf.MB_Sim(t, u, J, P, Q, 2.0, bond_dim; code = name, spin=true);


###############
# GROUNDSTATE #
###############

dictionary1 = hf.produce_groundstate(model1; force=Force, path=path);

dictionary2 = hf.produce_groundstate(model2; force=Force, path=path);

@testset "Groundstate" begin
    E_norm1 = -0.32637
    ψ1 = dictionary1["groundstate"];
    H1 = dictionary1["ham"];
    E01 = expectation_value(ψ1, H1);
    E1 = sum(real(E01))./length(H1);
    @test E1 ≈ E_norm1 atol=tol

    E_norm2 = -0.63093
    ψ2 = dictionary2["groundstate"];
    H2 = dictionary2["ham"];
    E02 = expectation_value(ψ2, H2);
    E2 = sum(real(E02))./length(H2);
    @test E2 ≈ E_norm2 atol=tol
end

###############
# EXCITATIONS #
###############

@testset "Excitations" begin 
    nums = 1;
    resolution = 5;
    momenta = range(0, π, resolution);

    exc = hf.produce_excitations(model1, momenta, nums; charges=[0,0.0,0], force=Force, path=path);
    Es = exc["Es"];
    @test imag(Es)≈zeros(size(Es)) atol=1e-8
end


#########
# Tools #
#########

@testset "Tools" begin
    N1 = hf.density_state(model1)
    Nup1, Ndown1 = hf.density_spin(model1)

    @test sum(N1)/2 ≈ sum(Nup1 + Ndown1)/length(dictionary1["ham"])

    N2 = hf.density_state(model2)
    Nup2, Ndown2 = hf.density_spin(model2)

    @test sum(N2)/2 ≈ sum(Nup2 + Ndown2)/length(dictionary1["ham"])
end