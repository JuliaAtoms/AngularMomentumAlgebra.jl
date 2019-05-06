# Tensors

```@meta
DocTestSetup = quote
    using AngularMomentumAlgebra
    using AtomicLevels
    using LinearAlgebra
end
```

```@docs
Tensor
TensorComponent
AngularMomentumAlgebra.LinearCombinationTensor
```

## The Wigner–Eckart theorem

The Wigner–Eckart theorem states that the matrix element of a tensor
component $\tensor{T}^{(k)}_q$ can be evaluated as

$$\begin{equation}
\matrixel{n'j'm'}{\tensor{T}^{(k)}_q}{njm}=
(-)^{j'-m'}
\begin{pmatrix}
j'&k&j\\
-m'&q&m
\end{pmatrix}
\redmatrixel{n'j'}{\tensor{T}^{(k)}}{nj},
\tag{V13.1.2}
\end{equation}$$

where the _reduced matrix element_
$\redmatrixel{n'j'}{\tensor{T}^{(k)}}{nj}$ does not depend on
$m,m'$.

```@docs
wigner_eckart
```

## [Spherical tensors](@id tensors_spherical_tensors)

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
\begin{pmatrix}
\ell&k&\ell'\\0&0&0
\end{pmatrix}.
\end{aligned}
\tag{V13.2.107}
```

```@docs
rme
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
