# * Spherical tensors

"""
    SphericalTensor(k)

Construct a spherical tensor of rank `k`.
"""
struct SphericalTensor{k} <: Tensor{k,'C'} end

"""
    system(::SphericalTensor)

A spherical tensor only acts on the coordinates ``\\theta`` and
``\\phi``.
"""
system(::SphericalTensor) = OrbitalAngularMomentumSubSystem()

@tensor(SphericalTensor{k} where k) do
    ℓ′ ∈ abs(ℓ - k):2:(ℓ+k)

    raw"""
    rme(ℓ′,𝐂̂ᵏ,ℓ)

Calculate the reduced matrix element of the spherical tensor of rank
`k`:

```math
\begin{aligned}
\redmatrixel{\ell'}{\tensor{C}^{(k)}}{\ell}
&=
\angroot{\ell}
C_{\ell 0;k,0}^{\ell'0} =
(-)^{\ell-k}
\angroot{\ell\ell'}
\wignerthreej{\ell&k&\ell'\\0&0&0}.
\end{aligned}
\tag{V13.2.107}
```
"""
    ∏(ℓ)*clebschgordan(ℓ,0,rank(k),0,ℓ′,0)
end

couples(a::SpinOrbital{<:Orbital}, ::Type{SphericalTensor}, b::SpinOrbital{<:Orbital}) =
    a.m[2] == b.m[2]

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
𝐂̂⁽¹⁾₀

julia> dot(SpinOrbital(o"1s", 0, half(1)), z, SpinOrbital(o"2p", 0, half(1)))
0.5773502691896256
```

"""
const 𝐫̂ = [cartesian_tensor_component(SphericalTensor(1), c)
           for c in [:x, :y, :z]]
export 𝐫̂
end

export SphericalTensor, rme
