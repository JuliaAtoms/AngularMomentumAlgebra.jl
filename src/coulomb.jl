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
    system(::CoulombTensor)

A Coulomb tensor only acts on the coordinates ``r``, ``\\theta`` and
``\\phi``.
"""
system(::CoulombTensor) = SpatialSubSystem()

@doc raw"""
    RadialCoulombMatrixElement

This represents the matrix element of the radial component of the
Coulomb tensor operator:

```math
\tensor{K}^{(k)}(i) \defd
\left\{[1-\Heaviside(r_j-r_i)]r_i^k +
\frac{\Heaviside(r_j-r_i)}{r_i^{k+1}}\right\}
\tensor{C}^{(k)}(i),
\quad
i = 1,2,
\quad
j = 3-i,
```
"""
struct RadialCoulombMatrixElement{k} <: OneBodyOperator end

RadialCoulombMatrixElement(k) = RadialCoulombMatrixElement{k}()

function Base.:(*)(a::RadialCoulombMatrixElement{k₁}, b::RadialCoulombMatrixElement{k₂}) where {k₁,k₂}
    @assert k₁ == k₂
    CoulombInteractionMultipole(k₁)
end

@tensor(CoulombTensor{k} where k) do
    begin
        n′ ~ n # The Coulomb interaction couples orbitals of different
               # n, but there is no selection rule.
        ℓ′ ∈ abs(ℓ - k):2:(ℓ+k)
    end

    raw"""
    rme((n′,ℓ′), ::CoulombTensor{k}, (n,ℓ))

Computes the reduced matrix element of `𝐊` in terms of
[`RadialCoulombMatrixElement`](@ref).
"""
    rme(ℓ′, SphericalTensor(k), ℓ) # *RadialCoulombMatrixElement(k)
end

"""
    ranks(a, ::Type{CoulombTensor}, b)

Return which tensor ranks for Coulomb tensors that fulfill the
triangle condition between spin-orbitals `a` and `b`.
"""
ranks(a::SpinOrbital, ::Type{CoulombTensor}, b::SpinOrbital) =
    triangle_range(a.orb.ℓ, b.orb.ℓ)

export CoulombInteractionMultipole, CoulombPotentialMultipole, CoulombTensor
