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
    ∏(ℓ)*clebschgordan(ℓ,0,k,0,ℓ′,0)
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

# * Dipole tensors

"""
    Dipole()

Construct a dipole tensor
"""
struct Dipole <: Tensor{1,'D'} end

"""
    system(::Dipole)

A dipole tensor only acts on the coordinates ``r``, ``\\theta`` and
``\\phi``.
"""
system(::Dipole) = SpatialSubSystem()

@doc raw"""
    RadialMatrixElement()

This represents the matrix element of the radial component of the
dipole operator:

```math
\expect{r} =
\int_0^\infty\diff{r}r^2
\conj{\Psi}_{n'\ell'}(r)
r
\Psi_{n\ell}(r)
```
"""
struct RadialMatrixElement <: OneBodyOperator end

Base.show(io::IO, ::RadialMatrixElement) = write(io, "r")

@tensor(Dipole) do
    begin
        n′ ~ n # The dipole couples orbitals of different n, but
               # there is no selection rule.
        ℓ′ == ℓ ± 1
    end

    raw"""
    rme((n′,ℓ′), ::Dipole, (n,ℓ))

Computes the reduced matrix element of `𝐃` in terms of
[`RadialMatrixElement`](@ref).
"""
    rme(ℓ′, SphericalTensor(1), ℓ)*RadialMatrixElement()
end

module Dipoles
import ..cartesian_tensor_component, ..SphericalTensor, ..Dipole

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

julia> dot(SpinOrbital(o"2s", 0, half(1)), z, SpinOrbital(o"2p", 0, half(1)))
0.5773502691896256
```

"""
const 𝐫̂ = [cartesian_tensor_component(SphericalTensor(1), c)
           for c in [:x, :y, :z]]

"""
    𝐫

The dipole operator; the elements correspond to `[x,y,z]`, i.e. the
Cartesian tensor components. Can be entered as `\\bfr`.

# Examples

```jldoctest
julia> using AngularMomentumAlgebra.Dipoles

julia> z = 𝐫[3]
𝐃̂⁽¹⁾₀

julia> dot(SpinOrbital(o"1s", 0, half(1)), z, SpinOrbital(o"2p", 0, half(1)))
0.57735r

julia> dot(SpinOrbital(o"2s", 0, half(1)), z, SpinOrbital(o"2p", 0, half(1)))
0.57735r

julia> dot(SpinOrbital(o"2p", 0, half(1)), z, SpinOrbital(o"3d", 0, half(1)))
0.516398r
```
"""
const 𝐫 = [cartesian_tensor_component(Dipole(), c)
           for c in [:x, :y, :z]]
export 𝐫̂, 𝐫
end

export SphericalTensor, Dipole, rme
