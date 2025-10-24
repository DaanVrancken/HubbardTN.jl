using Documenter
using HubbardTN

makedocs(
    sitename = "HubbardTN.jl",
    pages = [
        "Home" => "index.md",
        "Library" => "Functions.md",
        "Examples" => "Examples.md",
    ],
    format = Documenter.HTML(inventory_version = "0.2.0"),
)

deploydocs(
    repo = "github.com/DaanVrancken/HubbardTN.jl.git",
    push_preview = true,
    devbranch = "master",
)