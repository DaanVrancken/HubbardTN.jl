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
    if GROUP == "ALL" || GROUP == "OneBand"
        @time include("OneBand.jl")
    end
    if GROUP == "ALL" || GROUP == "MultiBand"
        @time include("MultiBand.jl")
    end

end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")

println("
Not included in tests:
- Equivalence of different symmetries.
- find_chemical_potential().
- compute_domainwall().
- Saving tools.
")