using Documenter
using HubbardTN

# For local usage: dirty workaround for limited acces rights of Julia.
# Might give troubles when ran on other system -> remove following 2 lines.
# If necessary, delete old "build" file manaually
# build_path = joinpath(pwd(),"build")
# if isdir(build_path)
#     run(`powershell -Command "Remove-Item -Recurse -Force '$build_path'"`)
# end

makedocs(
    sitename = "HubbardTN.jl",
    pages = [
        "Home" => "index.md",
        "Library" => "Functions.md",
        "Examples" => "Examples.md",
    ],
    format = Documenter.HTML(inventory_version = "0.1.0"),
)

deploydocs(
    repo = "github.com/DaanVrancken/HubbardTN.jl.git",
    push_preview = true,
    devbranch = "master",
)