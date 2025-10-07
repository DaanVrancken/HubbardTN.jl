using HubbardTN

using Test, Printf
using MPSKit, KrylovKit

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
    if GROUP == "ALL" || GROUP == "OB"
        @time include("OB.jl")
    end
    if GROUP == "ALL" || GROUP == "MB"
        @time include("MB.jl")
    end
    if GROUP == "ALL" || GROUP == "OBC"
        @time include("OBC.jl")
    end
    if GROUP == "ALL" || GROUP == "MBC"
        @time include("MBC.jl")
    end
    if GROUP == "ALL" || GROUP == "SPIN"
        @time include("Spin.jl")
    end
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")

println("
Not included in tests:
- Tools.
- Uiiij, Uijkk, Uijkl interaction.
")