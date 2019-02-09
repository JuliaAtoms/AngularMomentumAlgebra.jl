module AngularMomentumAlgebra

using AtomicLevels
using EnergyExpressions
import EnergyExpressions: QuantumOperator, NBodyTermFactor, NBodyTerm, OrbitalMatrixElement, NBodyMatrixElement
using WignerSymbols
using UnicodeFun

include("multipole_expansion.jl")

Base.Matrix(::QuantumOperator, csfs::Vector{<:CSF}, overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[]) =
    throw(ArgumentError("Derivation of energy expression for atomic CSFs not yet implemented"))

function Base.Matrix(op::QuantumOperator, spin_cfgs::VSC,
                     overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[]) where {VSC<:AbstractVector{<:Configuration{<:SpinOrbital}}}
    M = Matrix(op, SlaterDeterminant.(spin_cfgs), overlaps)
    m,n = size(M)
    Mm = zeros(eltype(M), m, n)
    for i = 1:m
        for j = 1:n
            Mm[i,j] = transform(multipole_expand, M[i,j])
        end
    end
    Mm
end

export # Re-exports from EnergyExpressions.jl
    Conjugate,
    OneBodyHamiltonian, CoulombInteraction

end # module
