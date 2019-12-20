# [Spherical tensors](@id tensors_spherical_tensors)

```@meta
DocTestSetup = quote
    using AngularMomentumAlgebra
    using AtomicLevels
    using LinearAlgebra
end
```

The spherical tensors are related to the spherical harmonics as

```math
\tensor{C}^{(k)}_q \defd
\sqrt{\frac{4\pi}{2k+1}}
Y^k_q.
\tag{V5.1.7}
```

```@docs
SphericalTensor
system(::SphericalTensor)
rme(ℓ′::Real,𝐂̂ᵏ::SphericalTensor,ℓ::Real)
AngularMomentumAlgebra.couplings(tensor::SphericalTensor,ℓ)
AngularMomentumAlgebra.ranks
```

### Dipole operator

The dipole operator is a rank-1 Cartesian tensor that may be expressed
using the rank-1 spherical tensor:

```math
\hat{\vec{r}} \equiv
\begin{bmatrix}\hat{x}\\\hat{y}\\\hat{z}\end{bmatrix}
\equiv
\begin{bmatrix}
\frac{1}{\sqrt{2}}[-\tensor{C}^{(1)}_1 + \tensor{C}^{(1)}_{-1}]\\
\frac{\im}{\sqrt{2}}[\tensor{C}^{(1)}_1 + \tensor{C}^{(1)}_{-1}]\\
\tensor{C}^{(1)}_0
\end{bmatrix}
```

This submodule exists only as a shortcut to [Cartesian tensor
components](@ref).

```@meta
CurrentModule = AngularMomentumAlgebra
```

```@docs
Dipoles.𝐫̂
```

```@meta
CurrentModule = nothing
DocTestSetup = nothing
```
