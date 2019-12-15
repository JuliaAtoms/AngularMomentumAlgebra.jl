# * Spherical tensors

"""
    SphericalTensor(k)

Construct a spherical tensor of rank `k`.
"""
struct SphericalTensor{k} <: Tensor{k,'C'} end

"""
    rme(ℓ′,Cᵏ,ℓ)

Calculate the reduced matrix element `⟨ℓ′||C⁽ᵏ⁾||ℓ⟩` of the spherical
tensor of rank `k`. Condon–Shortley phase convention and using the
definition of Eq. (13.2.107) in Varshalovich (1988).
"""
rme(ℓ′::Real,Cᵏ::SphericalTensor,ℓ::Real) = ∏(ℓ)*clebschgordan(ℓ,0,rank(Cᵏ),0,ℓ′,0)
rme(a::Orbital,Cᵏ::SphericalTensor,b::Orbital) = rme(a.ℓ, Cᵏ, b.ℓ)

"""
    rme(o′, Cᵏ, o)

Calculate the reduced matrix element `⟨o′||C⁽ᵏ⁾||o⟩` of the spherical
tensor of rank `k`, between the relativistic spin-orbitals `o′` and
`o`. Since the spherical tensors only act on the quantum numbers `ℓ′`
and `ℓ`, the reduced matrix has to be evaluated via the uncoupling
formula given in Eqs. (13.1.40) & (13.2.5) in Varshalovich (1988),
together with reduced matrix element for the uncoupled angular momenta
`ℓ′` and `ℓ` given by Eq. (13.2.107) (ibid).
"""
function rme(o′::RelativisticOrbital,Cᵏ::SphericalTensor,o::RelativisticOrbital)
    ℓ′,ℓ = o′.ℓ,o.ℓ
    k = rank(Cᵏ)
    iseven(ℓ′+k+ℓ) || return zero(Float64)

    j′,j = o′.j,o.j
    s = half(1)

    powneg1(Int(j+ℓ′+s+k))*∏(j,j′)*wigner6j(ℓ,  k, ℓ′,
                                            j′, s, j)*rme(ℓ,Cᵏ,ℓ′)

    # # @warn "Rederive rme!"
    
    # powneg1(Int(j′+s))*∏(j,j′)*wigner3j(j′, k,  j,
    #                                     s,  0, -s)
end

couples(a::SpinOrbital{<:Orbital}, ::Type{SphericalTensor}, b::SpinOrbital{<:Orbital}) =
    a.m[2] == b.m[2]

jmⱼ(o′::SpinOrbital{<:Orbital}, Cᵏ::SphericalTensor, o::SpinOrbital{<:Orbital}) =
    o′.orb.ℓ,o′.m[1],o.orb.ℓ,o.m[1],o′.m[2]==o.m[2]

"""
    ranks(a, ::Type{SphericalTensor}, b)

Return which tensor ranks for spherical tensors that fulfill the
triangle condition between spin-orbitals `a` and `b`.
"""
ranks(a::SpinOrbital, ::Type{SphericalTensor}, b::SpinOrbital) =
    triangle_range(a.orb.ℓ, b.orb.ℓ)

module Dipoles
import ..cartesian_tensor_component, ..SphericalTensor

"""
    𝐫̂

The angular part of the dipole operator; the elements correspond to
`[x,y,z]`, i.e. the Cartesian tensor components. Can be entered as
`\\bfr\\hat`.

# Examples

```jldoctest
julia> using AngularMomentumAlgebra.Dipoles

julia> z = 𝐫̂[3]
𝐂⁽¹⁾₀

julia> wigner_eckart(0, 0, z, 1, 0)
0.5773502691896256
```

"""
const 𝐫̂ = [cartesian_tensor_component(SphericalTensor(1), c)
           for c in [:x, :y, :z]]
export 𝐫̂
end

export SphericalTensor, rme
