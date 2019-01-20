module AngularMomentumAlgebra

using AtomicLevels
using WignerSymbols
using Combinatorics
using LinearAlgebra
using Formatting
using UnicodeFun
using Parameters
using Printf

include("symbolics.jl")
include("j.jl")
include("common.jl")
include("clebsch_gordan.jl")
include("kronecker.jl")
include("tensors.jl")
include("slater_integrals.jl")
include("conjugate_orbitals.jl")
include("slater_determinants.jl")
include("one_two_body_operators.jl")
include("energy_expressions.jl")
include("couplings.jl")

# Forward declaration, to be implemented on case-by-case basis by
# other packages. Should return a
# Pair{symmetry,AbstractMatrix{EnergyExpression}} (for a single
# argument) or a vector of such (for a vector argument).
energy_expression(::Any; verbosity=0) =
    error("Not implemented")

export energy_expression

end # module
