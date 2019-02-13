# Tensors

```@meta
DocTestSetup = quote
    using AngularMomentumAlgebra
    using AtomicLevels
end
```

```@docs
Tensor
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
