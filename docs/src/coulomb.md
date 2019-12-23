# Coulomb interaction

```@meta
DocTestSetup = quote
    using AngularMomentumAlgebra
    using AtomicLevels
end
```

The Coulomb interaction between two electrons at coordinates $1$ and
$2$, respectively, can in spherical coordinates be multipole-expanded
as:

```math
\begin{equation}
\frac{1}{r_{12}} =
\sum_{k=0}^\infty
\frac{r_<^k}{r_>^{k+1}}
P_k(\cos\theta)
\end{equation}
```

which by the _addition theorem_ for the spherical harmonics can be
further expanded as

```math
\begin{equation}
\frac{1}{r_{12}} =
\sum_{kq}
\frac{4\pi}{2k+1}
\frac{r_<^k}{r_>^{k+1}}
\conj{{Y^k_q}}(1)Y^k_q(2)\equiv
\sum_k
\frac{r_<^k}{r_>^{k+1}}
\tensor{C}^k(1)\cdot\tensor{C}^k(2),
\tag{V5.17.9}
\end{equation}
```

where we in the last step have used the definition of the [Spherical
tensors](@ref tensors_spherical_tensors). By introducing the
help-tensor

```math
\begin{equation}
\tensor{K}^{(k)}(i) \defd
\left\{[1-\Heaviside(r_j-r_i)]r_i^k +
\frac{\Heaviside(r_j-r_i)}{r_i^{k+1}}\right\}
\tensor{C}^{(k)}(i),
\quad
i = 1,2,
\quad
j = 3-i,
\end{equation}
```

where

```math
\begin{equation}
\Heaviside(x) = \begin{cases}
0, & x < 0\\
1, & x > 0,
\end{cases}
\end{equation}
```

is the Heaviside function, we can rewrite the Coulomb interaction as

```math
\begin{equation}
\frac{1}{r_{12}}=
\sum_k
\tensor{K}^{(k)}(1)\cdot\tensor{K}^{(k)}(2).
\tag{V5.17.9*}
\end{equation}
```

```@docs
CoulombInteractionMultipole
CoulombPotentialMultipole
CoulombTensor
system(::CoulombTensor)
AngularMomentumAlgebra.RadialCoulombMatrixElement
rme((n′,ℓ′), ::CoulombTensor, (n,ℓ))
AngularMomentumAlgebra.couplings(tensor::CoulombTensor, (n,ℓ)::Tuple{<:Number, <:Number})
```

```@meta
 DocTestSetup = nothing
```
