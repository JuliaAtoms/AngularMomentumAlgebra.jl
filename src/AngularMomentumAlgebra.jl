module AngularMomentumAlgebra

using AtomicLevels
using EnergyExpressions
import EnergyExpressions: QuantumOperator, NBodyTermFactor, NBodyTerm,
    OrbitalMatrixElement, NBodyMatrixElement
using LinearAlgebra
using SparseArrays
using WignerSymbols
using HalfIntegers
using UnicodeFun

include("common.jl")
include("clebsch_gordan.jl")
include("linear_combinations.jl")
include("orbitals.jl")
include("tensors.jl")
include("wigner_eckart.jl")
include("cartesian.jl")
include("angular_momenta.jl")
include("spherical_tensors.jl")
include("gradients.jl")
include("coulomb.jl")
include("multipole_expansion.jl")
include("energy_expressions.jl")

export # Re-exports from EnergyExpressions.jl
    Conjugate, OrbitalOverlap, EnergyExpression,
    FieldFreeOneBodyHamiltonian, OneBodyHamiltonian, CoulombInteraction, CoulombPotential,
    # Re-exports from LinearAlgebra
    ⋅, dot,
    # Own types
    CoulombInteractionMultipole, CoulombPotentialMultipole

end # module
