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

$$\begin{equation}
\frac{1}{r_{12}} =
\sum_{k=0}^\infty
\frac{r_<^k}{r_>^{k+1}}
P_k(\cos\theta)
\end{equation}$$

which by the _addition theorem_ for the spherical harmonics can be
further expanded as

$$\begin{equation}
\frac{1}{r_{12}} =
\sum_{kq}
\frac{4\pi}{2k+1}
\frac{r_<^k}{r_>^{k+1}}
\conj{{Y^k_q}}(1)Y^k_q(2)\equiv
\sum_{kq}
\frac{r_<^k}{r_>^{k+1}}
\tensor{C}^k_q(1)\cdot\tensor{C}^k_q(2),
\tag{V5.17.9}
\end{equation}$$

where we in the last step have used the definition of the [Spherical
tensors](@ref tensors_spherical_tensors).

```@docs
CoulombInteractionMultipole
CoulombPotentialMultipole
```

```@meta
 DocTestSetup = nothing
```
