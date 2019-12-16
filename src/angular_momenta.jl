abstract type AngularMomentum{label} <: Tensor{1,label} end

# * Orbital Angular Momentum

"""
    OrbitalAngularMomentum()

The angular momentum ``\\tensor{L}`` is associated with the spatial
coordinates ``\\theta`` and ``\\phi``.
"""
struct OrbitalAngularMomentum <: AngularMomentum{'L'} end
system(::OrbitalAngularMomentum) = OrbitalAngularMomentumSubSystem()

@doc raw"""
    rme(ℓ′, 𝐋̂, ℓ)

Calculate the reduced matrix element of the orbital angular momentum:

```math
\begin{equation}
\tag{V13.2.67}
\redmatrixel{\ell'}{\tensor{L}}{\ell} \defd
\delta_{\ell'\ell}\sqrt{\ell(\ell+1)(2\ell+1)}
\end{equation}
```

# Examples

```jldoctest
julia> rme(1, OrbitalAngularMomentum(), 1)
2.449489742783178

julia> rme(1, OrbitalAngularMomentum(), 2)
0
```
"""
function rme(ℓ′, ::OrbitalAngularMomentum, ℓ)
    @δ ℓ′,ℓ
    √(ℓ*(ℓ+1)*(2ℓ+1))
end

# * Spin Angular Momentum

"""
    SpinAngularMomentum()

The spin angular momentum ``\\tensor{S}`` is the intrinsic angular
momentum associated with the coordinate ``s``.
"""
struct SpinAngularMomentum <: AngularMomentum{'S'} end
system(::SpinAngularMomentum) = SpinSubSystem()

@doc raw"""
    rme(s′, ::SpinAngularMomentum, s)

Calculate the reduced matrix element of the spin angular momentum:

```math
\begin{equation}
\tag{V13.2.95}
\redmatrixel{s'}{\tensor{S}}{s} \defd
\delta_{ss'}
\sqrt{s(s+1)(2s+1)}
\end{equation}
```

# Examples

```jldoctest
julia> rme(half(1), SpinAngularMomentum(), half(1))
1.224744871391589

julia> rme(half(1), SpinAngularMomentum(), half(3))
0
```
"""
function rme(s′, ::SpinAngularMomentum, s)
    @δ s′,s
    √(s*(s+1)*(2s+1))
end

# * Total Angular Momentum

"""
    TotalAngularMomentum()

The total angular momentum ``\\tensor{J} = \\tensor{L} + \\tensor{S}``
results from the coupling of the orbital and spin angular momenta.
"""
struct TotalAngularMomentum <: AngularMomentum{'J'} end
system(::TotalAngularMomentum) = TotalAngularMomentumSubSystem()

@doc raw"""
    rme((ℓ′,s′,J′), ::TotalAngularMomentum, (ℓ,s,J))

Calculate the reduced matrix element of the total angular momentum:

```math
\begin{equation}
\tag{V13.2.38}
\redmatrixel{\ell's'J'}{\tensor{J}}{\ell s J} \defd
\delta_{\ell\ell'}\delta_{ss'}\delta_{JJ'}
\sqrt{J(J+1)(2J+1)}
\end{equation}
```

# Examples

```jldoctest
julia> rme((1,half(1),half(3)), TotalAngularMomentum(), (1,half(1),half(3)))
3.872983346207417

julia> rme((1,half(1),half(3)), TotalAngularMomentum(), (1,half(1),half(1)))
0
```

"""
function rme((ℓ′,s′,J′)::Tuple{<:Number, <:Number, <:Number},
             ::TotalAngularMomentum,
             (ℓ,s,J)::Tuple{<:Number, <:Number, <:Number})
    @δ ℓ′,ℓ s′,s J′,J
    √(J*(J+1)*(2J+1))
end

export OrbitalAngularMomentum, SpinAngularMomentum, TotalAngularMomentum
