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
-m'&\kappa&m
\end{pmatrix}
\redmatrixel{n'j'}{\tensor{T}^{(k)}_q}{nj},
\tag{V13.1.1}
\end{equation}$$

where the _reduced matrix element_
$\redmatrixel{n'j'}{\tensor{T}^{(k)}_q}{nj}$ does not depend on
$m,m'$.

```@docs
wigner_eckart
```

## [Spherical tensors](@id tensors_spherical_tensors)

The spherical tensors are related to the spherical harmonics as

$$\begin{equation}
\tensor{C}^{(k)}_q \defd
\sqrt{\frac{4\pi}{2k+1}}
Y^k_q.
\tag{V5.1.7}
\end{equation}$$

```@docs
SphericalTensor
```

The reduced matrix element of the spherical tensor is given by

$$\begin{equation}
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
\end{equation}$$

```@docs
rme
```

```@meta
 DocTestSetup = nothing
```
