using HubbardTN

using Test, Printf
using MPSKit, KrylovKit, TensorKit

# check if user supplied args --group="..."
pat = r"(?:--group=)(\w+)"
arg_id = findfirst(contains(pat), ARGS)
const GROUP = if isnothing(arg_id)
    uppercase(get(ENV, "GROUP", "ALL"))
else
    uppercase(only(match(pat, ARGS[arg_id]).captures))
end

# Run test suite
println("Starting tests")
ti = time()

@time begin
    if GROUP == "ALL" || GROUP == "ONEBAND"
        @time include("OneBand.jl")
    end
    if GROUP == "ALL" || GROUP == "MULTIBAND"
        @time include("MultiBand.jl")
    end
    if GROUP == "ALL" || GROUP == "EXCITATIONS"
        @time include("Excitations.jl")
    end
    if GROUP == "ALL" || GROUP == "THREEBODY"
        @time include("ThreeBody.jl")
    end
    if GROUP == "ALL" || GROUP == "CHEMICALPOTENTIAL"
        @time include("ChemicalPotential.jl")
    end
    if GROUP == "ALL" || GROUP == "STAGGEREDFIELD"
        @time include("StaggeredField.jl")
    end
    if GROUP == "ALL" || GROUP == "SPINFIELD"
        @time include("SpinField.jl")
    end
    if GROUP == "ALL" || GROUP == "CHARGEGAPMF"
        @time include("ChargeGapMF.jl")
    end
    if GROUP == "ALL" || GROUP == "HOLSTEIN"
        @time include("Holstein.jl")
    end
    if GROUP == "ALL" || GROUP == "FINITEMPS"
        @time include("FiniteMPS.jl")
    end
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")

println("
Not included in tests:
- Several combinations of terms in Hamiltonian.
- PairGapMF
- compute_domainwall(), get_alpha(), get_beta(), density_correlations().
- Saving tools.
")