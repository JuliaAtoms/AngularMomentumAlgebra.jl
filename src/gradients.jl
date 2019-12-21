"""
    Gradient()

Construct a gradient tensor.
"""
struct Gradient <: Tensor{1,'∇'} end

"""
    system(::Gradient)

The gradient only acts on the coordinates ``r``, ``\\theta``, and
``\\phi``.
"""
system(::Gradient) = SpatialSubSystem()

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

export Gradient
