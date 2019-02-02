module AngularMomentumAlgebra

using AtomicLevels
using WignerSymbols
using Combinatorics
using LinearAlgebra
using SparseArrays
using UnicodeFun
using Printf

include("conjugate_orbitals.jl")
include("one_body.jl")
include("repulsion_potentials.jl")
include("two_body.jl")
include("multipole_expansion.jl")
include("energy_expressions.jl")
include("invariant_sets.jl")

include("slater_determinants.jl")

# Forward declaration, to be implemented on case-by-case basis by
# other packages. Should return a
# Pair{symmetry,AbstractMatrix{EnergyExpression}} (for a single
# argument) or a vector of such (for a vector argument).
energy_expression(::Any; verbosity=0) =
    error("Not implemented")

export energy_expression

end # module
