println("
###################
#  One-Band Chem  #
###################
")

Force = false
tol = 1e-1

# store calculations at
path = joinpath("sims","OBC")


###############
# GROUNDSTATE #
###############

@testset "Groundstate" begin
    P=1; 
    Q=1;

    model = OBC_Sim([1.0], [1.0], P/Q, 2.0; mu=false, verbosity_mu=1);

    E_norm = -1.03541433

    dictionary = produce_groundstate(model; force=Force, path=path);
    ψ₀ = dictionary["groundstate"];
    H = dictionary["ham"];

    Ne = density_state(ψ₀);
    E0 = sum(expectation_value(ψ₀, H)) + dictionary["μ"]*Ne;

    @test Ne≈P/Q atol=1e-6
    @test real(E0)≈E_norm atol=tol
end


###############
# EXCITATIONS #
###############

μ= 3.535;
u=[7.658];
t=[2.726];

model = OBC_Sim(t, u, μ, 2.0);

dictionary = produce_groundstate(model; force=Force, path=path);

@testset "Excitations" begin
    ψ₀ = dictionary["groundstate"];
    H = dictionary["ham"];
    Ne = density_state(ψ₀);
    E0 = sum(expectation_value(ψ₀, H)) + dictionary["μ"]*Ne

    Es_norm = [4.13541796; 2.8491043; -0.4113358; 2.89508166; 4.17185897]

    nums = 1;
    resolution = 5;
    momenta = range(0, π, resolution);

    exc = produce_excitations(model, momenta, nums; force=Force, charges=[1,0.5,1], path=path);
    Es = exc["Es"];
    @test real(Es)≈Es_norm atol=tol
    @test imag(Es)≈zeros(size(Es)) atol=1e-8
end


#########
# Tools #
#########

@testset "Tools" begin
    D = dim_state(dictionary["groundstate"])
    @test typeof(D) == Vector{Int64}
    @test D > zeros(size(D))
end
