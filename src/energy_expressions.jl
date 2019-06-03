Base.Matrix(::QuantumOperator, csfs::Vector{<:CSF},
            overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[]) =
    throw(ArgumentError("Derivation of energy expression for atomic CSFs not yet implemented"))

"""
    Matrix(op::QuantumOperator,
           spin_cfgs::Vector{<:Configuration{<:SpinOrbital}}[, overlaps])

Generate the energy-expression associated with the quantum operator
`op`, in the basis of the spin-configurations `spin_cfgs`, with an
optional set of orbital `overlaps`, specifying any desired
non-orthogonalities. The energy expression is generated in a
basis-agnostic way by EnergyExpressions.jl and then
multipole-expanded.
"""
function Base.Matrix(op::QuantumOperator, spin_cfgs::VSC,
                     overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[]) where {VSC<:AbstractVector{<:Configuration{<:SpinOrbital}}}
    M = Matrix(op, SlaterDeterminant.(spin_cfgs), overlaps)
    m,n = size(M)
    Mm = spzeros(eltype(M), m, n)
    for i = 1:m
        for j = 1:n
            me = transform(multipole_expand, M[i,j])
            iszero(me) && continue
            Mm[i,j] = me
        end
    end
    Mm
end
