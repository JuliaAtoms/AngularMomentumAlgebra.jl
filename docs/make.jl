using Documenter
using AngularMomentumAlgebra
using AtomicLevels
using HalfIntegers
using LinearAlgebra

DocMeta.setdocmeta!(AngularMomentumAlgebra, :DocTestSetup,
                    :(using AngularMomentumAlgebra, AtomicLevels, HalfIntegers, LinearAlgebra); recursive=true)

makedocs(
    modules = [AngularMomentumAlgebra],
    format = Documenter.HTML(
        mathengine = MathJax2(Dict(:TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(
                :defd => "≝",
                :abs => ["|#1|",1],
                :ket => ["|#1\\rangle",1],
                :bra => ["\\langle#1|",1],
                :braket => ["\\langle#1|#2\\rangle",2],
                :matrixel => ["\\langle#1|#2|#3\\rangle",3],
                :redmatrixel => ["\\langle#1||#2||#3\\rangle",3],
                :expect => ["\\langle#1\\rangle",1],
                :vec => ["\\mathbf{#1}",1],
                :mat => ["\\mathsf{#1}",1],
                :conj => ["#1^*",1],
                :im => "\\mathrm{i}",
                :tensor => ["\\hat{\\mathbf{#1}}",1],
                :operator => ["\\mathfrak{#1}",1],
                :Hamiltonian => "\\operator{H}",
                :hamiltonian => "\\operator{h}",
                :Lagrangian => "\\operator{L}",
                :fock => "\\operator{f}",
                :lagrange => ["\\epsilon_{#1}",1],
                :vary => ["\\delta_{#1}",1],
                :onebody => ["(#1|#2)",2],
                :twobody => ["[#1|#2]",2],
                :twobodydx => ["[#1||#2]",2],
                :direct => ["{\\operator{J}_{#1}}",1],
                :exchange => ["{\\operator{K}_{#1}}",1],
                :diff => ["\\mathrm{d}#1\\,",1],
                :angroot => ["{\\prod}_{#1}", 1],
                :wignerthreej => ["\\begin{pmatrix}#1\\end{pmatrix}", 1],
                :wignersixj => ["\\begin{Bmatrix}#1\\end{Bmatrix}", 1],
                :Heaviside => "\\Theta"
            ),
        ))),
    ),
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
    doctest = false,
)

deploydocs(repo = "github.com/JuliaAtoms/AngularMomentumAlgebra.jl.git",
           push_preview = true,)
