using Documenter
using AngularMomentumAlgebra

makedocs(
    modules = [AngularMomentumAlgebra],
    format = Documenter.HTML(assets = ["assets/latex.js"],
                             mathengine = Documenter.MathJax()),
    sitename = "AngularMomentumAlgebra",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
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
    repo="https://github.com/JuliaAtoms/AngularMomentumalgebra.jl/blob/{commit}{path}#L{line}",
)

deploydocs(repo = "github.com/JuliaAtoms/AngularMomentumAlgebra.jl.git")
