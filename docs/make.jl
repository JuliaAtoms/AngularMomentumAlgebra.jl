using Documenter
using AngularMomentumAlgebra
using AtomicLevels
using HalfIntegers

DocMeta.setdocmeta!(AngularMomentumAlgebra, :DocTestSetup,
                    :(using AngularMomentumAlgebra, AtomicLevels, HalfIntegers, LinearAlgebra); recursive=true)

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
        "Orbitals" => "orbitals.md",
        "Tensors" => [
            "General tensors" => "tensors.md",
            "Tensor DSL" => "tensor_dsl.md",
            "Various tensors" => [
                "Angular momenta" => "angular_momenta.md",
                "Spherical tensors" => "spherical_tensors.md",
                "Gradients" => "gradients.md",
                "Multipole expansions" => "multipole_expansions.md",
                "Coulomb interaction" => "coulomb.md",
            ],
            "Tensor matrix elements" => "tensor_matrix_elements.md",
        ],
        "Energy expressions" => "energy_expressions.md",
    ],
    repo="https://github.com/JuliaAtoms/AngularMomentumalgebra.jl/blob/{commit}{path}#L{line}",
)

deploydocs(repo = "github.com/JuliaAtoms/AngularMomentumAlgebra.jl.git")
