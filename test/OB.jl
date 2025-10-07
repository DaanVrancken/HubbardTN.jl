println("
##############
#  One-Band  #
##############
")

Force = false
tol = 1e-2

# store calculations at
path = joinpath("sims","OB")


############################
# DEPENDENCE ON PARAMETERS #
############################

U_max = 2;

u_range = 0:1:U_max;
E = zeros(U_max+1,1);
P = 1;
Q = 1;

E_norm = [-1.2696767, -1.037173, -0.84163698]

@testset "Dependence on parameters" for u in u_range
    model = OB_Sim([1.0], [Float64(u)], 0.0, [0.0], P, Q, 2.0);
    dictionary = produce_groundstate(model; force=Force, path=path);
    ψ₀ = dictionary["groundstate"];
    H = dictionary["ham"];
    E0 = expectation_value(ψ₀, H);
    E[u+1] = sum(real(E0))/length(ψ₀);
    @test E[u+1] ≈ E_norm[u+1] atol=tol
end


#########################
# DEPENDENCE ON FILLING #
#########################

t = [1.0];
u = [5.0];
P = [1, 1, 3];
Q = [2, 1, 2];
E = zeros(length(P),1);

E_norm = [-0.73920032, -0.48460447, 1.76073968]

@testset "Dependence on filling" for i in eachindex(P)
    model = OB_Sim(t, u, 0.0, [0.0], P[i], Q[i], 2.0);
    dictionary = produce_groundstate(model; force=Force, path=path);
    ψ₀ = dictionary["groundstate"];
    H = dictionary["ham"];
    E0 = expectation_value(ψ₀, H);
    E[i] = sum(real(E0))/length(ψ₀);
    @test E[i] ≈ E_norm[i] atol=tol
end


###############
# EXCITATIONS #
###############

P=1;
Q=1;
u=[5.0];
t=[1.0];

model = OB_Sim(t, u, 0.0, [0.0], P, Q, 2.0);

dictionary = produce_groundstate(model; force=Force, path=path);

@testset "Excitations" begin 
    ψ₀ = dictionary["groundstate"];
    H = dictionary["ham"];
    E0 = expectation_value(ψ₀, H);
    E = sum(real(E0))/length(H)

    Es_norm = [-0.17257389; -0.2673373; -0.5489149; -1.04588404; -1.425526126]

    nums = 1;
    resolution = 5;
    momenta = range(0, π, resolution);

    exc = produce_excitations(model, momenta, nums; force=Force, charges=[1,0.5,-1], path=path);
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

    electron_number = sum(density_state(model; path=path))/2
    @test sum(electron_number)≈P/Q atol=1e-8
end