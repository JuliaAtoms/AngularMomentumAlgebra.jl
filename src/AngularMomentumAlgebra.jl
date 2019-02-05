module AngularMomentumAlgebra

using AtomicLevels
using EnergyExpressions
using WignerSymbols

include("multipole_expansion.jl")

# Forward declaration, to be implemented on case-by-case basis by
# other packages. Should return a
# Pair{symmetry,AbstractMatrix{EnergyExpression}} (for a single
# argument) or a vector of such (for a vector argument).
energy_expression(::Any; verbosity=0) =
    error("Not implemented")

export energy_expression,
    # Re-exports from EnergyExpressions.jl
    Conjugate,
    OneBodyHamiltonian, OneBodyIntegral,
    RepulsionPotential, DirectPotential, ExchangePotential, DirectExchangePotentials,
    GeneralRepulsionIntegral, DirectIntegral, ExchangeIntegral, TwoBodyIntegral,
    isdiagonal,
    OneBodyEnergyExpression, TwoBodyEnergyExpression,
    adjacent_in_coord!, diagonal_in_coord!,
    adjacent_in_coord, diagonal_in_coord,
    hamiltonian_matrix, one_body_hamiltonian_matrix, two_body_hamiltonian_matrix,
    coupled_states, invariant_sets


end # module
