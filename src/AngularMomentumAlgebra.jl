module AngularMomentumAlgebra

using AtomicLevels
using WignerSymbols
using Formatting
using UnicodeFun
using LinearAlgebra
using Parameters
using Printf

include("symbolics.jl")
include("j.jl")
include("common.jl")
include("clebsch_gordan.jl")
include("kronecker.jl")
include("tensors.jl")
include("slater_integrals.jl")
include("couplings.jl")

# Forward declaration, to be implemented on case-by-case basis by
# other packages. Should return a
# Pair{symmetry,AbstractMatrix{EnergyExpression}} (for a single
# argument) or a vector of such (for a vector argument).
energy_expression(::Any; verbosity=0) =
    error("Not implemented")

export energy_expression

end # module
