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

export CoulombInteractionMultipole, CoulombPotentialMultipole
