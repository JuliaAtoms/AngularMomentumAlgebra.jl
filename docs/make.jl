using Documenter
using AngularMomentumAlgebra

makedocs(
    modules = [AngularMomentumAlgebra],
    sitename = "AngularMomentumAlgebra",
    pages = [
        "Home" => "index.md",
        "Common routines" => "common.md",
        "Clebsch–Gordan coefficients" => "clebsch_gordan.md",
        "Tensors" => "tensors.md",
        "Coulomb interaction" => "coulomb.md",
        "Multipole expansions" => "multipole_expansions.md",
        "Energy expressions" => "energy_expressions.md",
    ],
    assets = ["assets/latex.js"]
)

deploydocs(repo = "github.com/JuliaAtoms/AngularMomentumAlgebra.jl.git")
