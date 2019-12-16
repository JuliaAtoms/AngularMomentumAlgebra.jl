# [Gradients](@id gradients)

The gradient is a rank-1 tensor, the reduced matrix element of which
is given by

```math
\begin{equation}
\tag{V13.2.21}
\redmatrixel{n'\ell'}{\tensor{\nabla}^{(1)}}{n\ell} =
\sqrt{\ell+1}
A_{n'\ell'n\ell}\delta_{\ell',\ell+1} +
\sqrt{\ell}
B_{n'\ell'n\ell}\delta_{\ell',\ell-1},
\end{equation}
```

where

```math
\begin{equation}
\tag{V13.2.22}
\begin{aligned}
A_{n'\ell'n\ell} &\defd
\int_0^\infty\diff{r}r^2
\conj{\Psi}_{n'\ell'}(r)
\left(\partial_r-\frac{\ell}{r}\right)
\Psi_{n\ell}(r),\\
B_{n'\ell'n\ell} &\defd
\int_0^\infty\diff{r}r^2
\conj{\Psi}_{n'\ell'}(r)
\left(\partial_r+\frac{\ell+1}{r}\right)
\Psi_{n\ell}(r).
\end{aligned}
\end{equation}
```

## Example

```jldoctest
julia> using AngularMomentumAlgebra, AtomicLevels

julia> orbitals = sos"k[s-d]"
18-element Array{SpinOrbital{Orbital{Symbol},Tuple{Int64,HalfIntegers.Half{Int64}}},1}:
 ks₀α
 ks₀β
 kp₋₁α
 kp₋₁β
 kp₀α
 kp₀β
 kp₁α
 kp₁β
 kd₋₂α
 kd₋₂β
 kd₋₁α
 kd₋₁β
 kd₀α
 kd₀β
 kd₁α
 kd₁β
 kd₂α
 kd₂β

julia> a,b,c = orbitals[[3,9,13]]
3-element Array{SpinOrbital{Orbital{Symbol},Tuple{Int64,HalfIntegers.Half{Int64}}},1}:
 kp₋₁α
 kd₋₂α
 kd₀α

julia> ∂x = cartesian_tensor_component(Gradient(), :x)
- 0.707107 𝛁̂⁽¹⁾₁ + 0.707107 𝛁̂⁽¹⁾₋₁

julia> dot(a, ∂x, b)
0.447214(∂ᵣ + 3/r)

julia> dot(b, ∂x, a)
0.447214(∂ᵣ - 1/r)

julia> dot(a, ∂x, c)
- 0.182574(∂ᵣ + 3/r)

julia> dot(c, ∂x, a)
- 0.182574(∂ᵣ - 1/r)
```

## Reference

```@docs
Gradient
system(::Gradient)
AngularMomentumAlgebra.RadialGradientMatrixElement
```
