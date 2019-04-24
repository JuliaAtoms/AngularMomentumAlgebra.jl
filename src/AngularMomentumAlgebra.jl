module AngularMomentumAlgebra

using AtomicLevels
using EnergyExpressions
import EnergyExpressions: QuantumOperator, NBodyTermFactor, NBodyTerm,
    OrbitalMatrixElement, NBodyMatrixElement
using WignerSymbols
using UnicodeFun

include("common.jl")
include("clebsch_gordan.jl")
include("linear_combinations.jl")
include("tensors.jl")
include("coulomb.jl")
include("multipole_expansion.jl")
include("energy_expressions.jl")

export # Re-exports from EnergyExpressions.jl
    Conjugate, OrbitalOverlap, EnergyExpression,
    FieldFreeOneBodyHamiltonian, OneBodyHamiltonian, CoulombInteraction, CoulombPotential,
    # Own types
    CoulombInteractionMultipole, CoulombPotentialMultipole

end # module
