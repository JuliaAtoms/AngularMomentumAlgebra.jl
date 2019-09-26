using Documenter
using AngularMomentumAlgebra

makedocs(
    modules = [AngularMomentumAlgebra],
    sitename = "AngularMomentumAlgebra",
    pages = [
        "Home" => "index.md",
        "Definitions" => "definitions.md",
        "Common routines" => "common.md",
        "Tensors" => "tensors.md",
        "Spherical Tensors" => "spherical_tensors.md",
        "Coulomb interaction" => "coulomb.md",
        "Multipole expansions" => "multipole_expansions.md",
        "Energy expressions" => "energy_expressions.md",
    ],
    assets = ["assets/latex.js"]
)

deploydocs(repo = "github.com/JuliaAtoms/AngularMomentumAlgebra.jl.git")
