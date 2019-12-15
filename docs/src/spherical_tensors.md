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
```

The reduced matrix element of the spherical tensor is given by

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
For the reduced matrix element of the same spherical tensor in the
``\ell s j`` basis, we must use the uncoupling formula together with
the above formula in the uncoupled basis:
```math
\begin{aligned}
\redmatrixel{\ell' s' j'}{\tensor{C}^{(k)}}{\ell s j}
&=
\delta_{ss'}(-)^{j+\ell'+s+k}
\angroot{jj'}
\wignersixj{\ell&s&j\\j'&k&\ell'}
\redmatrixel{\ell'}{\tensor{C}^{(k)}}{\ell}\\
&=
(-)^{\ell'+\ell+s+j}
\angroot{\ell\ell'jj'}
\wignersixj{\ell&s&j\\j'&k&\ell'}
\wignerthreej{\ell&k&\ell'\\0&0&0}.
\end{aligned}
\tag{V13.\{1.40,2.5\}}
```

```@docs
rme
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
