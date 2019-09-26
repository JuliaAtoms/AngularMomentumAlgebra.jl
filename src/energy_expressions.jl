Base.Matrix(::QuantumOperator, csfs::Vector{<:CSF},
            overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[]) =
    throw(ArgumentError("Derivation of energy expression for atomic CSFs not yet implemented"))

"""
    integrate_spinor(me)

Perform the spin-angular integration of the matrix element `me`,
leaving only a radial integral multiplied by a spin-angular
coefficient. The spin-angular integral is dependent on the specific
combination of spin-orbitals and the operator (expressed as a tensor);
the default implementation is to leave `me` as-is, corresponding to a
spin-angular integral of unity.
"""
integrate_spinor(me) = me

"""
    integrate_spinor(me::OrbitalMatrixElement{2,<:Any,<:CoulombInteraction,<:Any})

Perform the spin-angular integration of the two-electron matrix
element `me`, by first multipole-expanding the Coulomb interaction and
then integrating all the resulting terms over the spin-angular
coordinates (see [`multipole_expand`](@ref)).
"""
integrate_spinor(me::OrbitalMatrixElement{2,<:Any,<:CoulombInteraction,<:Any}) =
    multipole_expand(me)

"""
    Matrix(op::QuantumOperator,
           spin_cfgs::Vector{<:Configuration{<:SpinOrbital}}[, overlaps])

Generate the energy-expression associated with the quantum operator
`op`, in the basis of the spin-configurations `spin_cfgs`, with an
optional set of orbital `overlaps`, specifying any desired
non-orthogonalities. The energy expression is generated in a
basis-agnostic way by EnergyExpressions.jl and each term is then
integrated over the spin-angular coordinates using
[`integrate_spinor`](@ref).
"""
function Base.Matrix(op::QuantumOperator, spin_cfgs::VSC,
                     overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[]) where {VSC<:AbstractVector{<:Configuration{<:SpinOrbital}}}
    M = Matrix(op, SlaterDeterminant.(spin_cfgs), overlaps)
    m,n = size(M)
    Mm = spzeros(eltype(M), m, n)
    for i = 1:m
        for j = 1:n
            # println(repeat("=", 100))
            # @show i,j
            me = transform(integrate_spinor, M[i,j])
            iszero(me) && continue
            Mm[i,j] = me
        end
    end
    Mm
end
