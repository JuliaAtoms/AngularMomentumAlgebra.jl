"""
    CoulombInteractionMultipole(k, g)

Represents the `k`th multipole of the multipole expansion of the
Coulomb interaction `g`.
"""
struct CoulombInteractionMultipole{G<:CoulombInteraction} <: TwoBodyOperator
    k::Int
    g::G
end

"""
    CoulombPotentialMultipole

Type alias for contraction of the
[`CoulombInteractionMultipole`](@ref) over one coordinate, thereby
forming a potential in the other coordinate.
"""
const CoulombPotentialMultipole{A,B} = ContractedOperator{1,2,1,A,<:CoulombInteractionMultipole,B}

Base.show(io::IO, ci::CoulombInteractionMultipole) = write(io, "ĝ", to_superscript(ci.k))

function Base.show(io::IO, me::OrbitalMatrixElement{2,A,<:CoulombInteractionMultipole,B}) where {A,B}
    if me.a == me.b # Direct interaction
        write(io, "F",to_superscript(me.o.k),"($(me.a[1]),$(me.a[2]))")
    elseif me.a[1] == me.b[2] && me.a[2] == me.b[1] # Exchange interaction
        write(io, "G",to_superscript(me.o.k),"($(me.a[1]),$(me.a[2]))")
    else # General case
        write(io, "R",to_superscript(me.o.k),"(", join(string.(me.a), ","), ";", join(string.(me.b), ","), ")")
    end
end

Base.show(io::IO, me::CoulombPotentialMultipole{A,B}) where {A,B}=
    write(io, "r⁻¹×Y",to_superscript(me.o.k),"($(me.a[1]),$(me.b[1]))")

radial_integral(a, (k,g)::Tuple{<:Integer,<:CoulombInteraction}, b) =
    OrbitalMatrixElement(a, CoulombInteractionMultipole(k,g), b)

"""
    CoulombTensor(k)

Construct a Coulomb interaction tensor of rank `k`.
"""
struct CoulombTensor{k} <: Tensor{k,'K'} end

"""
    system(::Type{CoulombTensor})

A Coulomb tensor only acts on the coordinates ``r``, ``\\theta`` and
``\\phi``.
"""
system(::Type{<:CoulombTensor}) = SpatialSubSystem()

@tensor(CoulombTensor{k} where k) do
    begin
        n′ ~ n # The Coulomb interaction couples orbitals of different
               # n, but there is no selection rule.
        ℓ′ ∈ abs(ℓ - k):2:(ℓ+k)
    end

    raw"""
    rme((n′,ℓ′), ::CoulombTensor{k}, (n,ℓ))

Computes the reduced matrix element of `𝐊`.
"""
    rme(ℓ′, SphericalTensor(k), ℓ)
end

"""
    ranks(a, ::Type{CoulombTensor}, b)

Return which tensor ranks for Coulomb tensors that fulfill the
triangle condition between spin-orbitals `a` and `b`.
"""
ranks(a::SpinOrbital, ::Type{CoulombTensor}, b::SpinOrbital) =
    triangle_range(a.orb.ℓ, b.orb.ℓ)

export CoulombInteractionMultipole, CoulombPotentialMultipole, CoulombTensor
