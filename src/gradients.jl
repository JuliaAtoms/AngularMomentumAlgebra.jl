"""
    Gradient()

Construct a gradient tensor.
"""
struct Gradient <: Tensor{1,'∇'} end

"""
    system(::Type{Gradient})

The gradient only acts on the coordinates ``r``, ``\\theta``, and
``\\phi``.
"""
system(::Type{Gradient}) = SpatialSubSystem()

"""
    ReducedGradient()

Construct a gradient tensor acting on reduced wavefunctions.
"""
struct ReducedGradient <: Tensor{1,'∂'} end

"""
    system(::Type{ReducedGradient})

The reduced gradient only acts on the coordinates ``r``, ``\\theta``,
and ``\\phi``.
"""
system(::Type{ReducedGradient}) = SpatialSubSystem()

@doc raw"""
    RadialGradientMatrixElement(k)

This represents the matrix element of the radial component of the
gradient operator:

```math
\begin{equation}
\tag{V13.2.23}
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

@tensor(Gradient) do
    begin
        n′ ~ n # The gradient couples orbitals of different n, but
               # there is no selection rule.
        ℓ′ == ℓ ± 1
    end

    raw"""
    rme((n′,ℓ′), ::Gradient, (n,ℓ))

Computes the reduced matrix element of `∇` in terms of
[`RadialGradientMatrixElement`](@ref).
"""
    if ℓ′ == ℓ+1
        √(ℓ+1)*RadialGradientMatrixElement(-ℓ)
    elseif ℓ′==ℓ-1
        - √(ℓ)*RadialGradientMatrixElement(ℓ+1)
    end
end

@tensor(ReducedGradient) do
    begin
        n′ ~ n # The gradient couples orbitals of different n, but
               # there is no selection rule.
        ℓ′ == ℓ ± 1
    end

    raw"""
    rme((n′,ℓ′), ::ReducedGradient, (n,ℓ))

Computes the reduced matrix element of `∂` in terms of
[`RadialGradientMatrixElement`](@ref).
"""
    if ℓ′ == ℓ+1
        √(ℓ+1)*RadialGradientMatrixElement(-(ℓ+1))
    elseif ℓ′==ℓ-1
        - √(ℓ)*RadialGradientMatrixElement(ℓ)
    end
end

module LinearMomenta
import ..cartesian_tensor_component, ..Gradient, ..ReducedGradient

@doc raw"""
    𝐩

The linear momentum operator ``\vec{p}=-\im\nabla``; the elements
correspond to `[px,py,pz]`, i.e. the Cartesian tensor components. Can
be entered as `\bfp`.
"""
const 𝐩 = [-im*cartesian_tensor_component(Gradient(), c)
           for c in [:x, :y, :z]]

@doc raw"""
    𝐩̃

The linear momentum operator ``\vec{p}=-\im\nabla``, but _evaluated in
the basis of reduced wavefunctions_; the elements correspond to
`[px,py,pz]`, i.e. the Cartesian tensor components. Can be entered as
`\bfp\tilde`.
"""
const 𝐩̃ = [-im*cartesian_tensor_component(ReducedGradient(), c)
           for c in [:x, :y, :z]]
export 𝐩, 𝐩̃
end

export Gradient, ReducedGradient
