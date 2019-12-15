"""
    Gradient()

Construct a gradient tensor.
"""
struct Gradient <: Tensor{1,'∇'} end

@doc raw"""
    RadialGradientMatrixElement(k)

This represents the matrix element of the radial component of the
gradient operator:

```math
\begin{equation}
\tag{V13.2.24}
\int_0^\infty\diff{r}r^2
\conj{\Psi}_{n'\ell'}(r)
\left(\partial_r+\frac{k}{r}\right)
\Psi_{n\ell}(r)
\end{equation}
```
"""
struct RadialGradientMatrixElement <: OneBodyOperator
    k::Int
end

function Base.show(io::IO, me::RadialGradientMatrixElement)
    if iszero(me.k)
        write(io, "∂ᵣ")
    else
        write(io, "(∂ᵣ ")
        write(io, me.k > 0 ? "+" : "-")
        write(io, " $(abs(me.k))/r)")
    end
end


function rme(ℓ′::Real,∇::Gradient,ℓ::Real)
    if ℓ′ == ℓ+1
        √(ℓ+1)*RadialGradientMatrixElement(-ℓ)
    elseif ℓ′==ℓ-1
        - √(ℓ)*RadialGradientMatrixElement(ℓ+1)
    else
        0
    end
end
rme(a::Orbital,∇::Gradient,b::Orbital) = rme(a.ℓ, ∇, b.ℓ)

jmⱼ(o′::SpinOrbital{<:Orbital}, ∇::Gradient, o::SpinOrbital{<:Orbital}) =
    o′.orb.ℓ,o′.m[1],o.orb.ℓ,o.m[1],o′.m[2]==o.m[2]

export Gradient
