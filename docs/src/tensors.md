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

## Product tensors

```@docs
TensorProduct
```

### Scalar product

The matrix element of scalar product of two tensors acting on
different coordinates is given by

$$\begin{equation}
\begin{aligned}
&\matrixel{n_aj_am_a;n_bj_bm_b}{[\tensor{P}^{(k)}(1)\cdot\tensor{Q}^{(k)}(2)]}{n_cj_cm_c;n_dj_dm_d}\\
=&
\frac{1}{\angroot{j_aj_b}}
\sum_\alpha(-)^{-\alpha}
C_{j_cm_c;k,\alpha}^{j_am_a}
C_{j_dm_d;k,-\alpha}^{j_bm_b}\\
&\times
\redmatrixel{n_aj_a}{\tensor{P}^{(k)}(1)}{n_cj_c}
\redmatrixel{n_bj_b}{\tensor{Q}^{(k)}(2)}{n_dj_d}
\end{aligned}
\tag{V13.1.26}
\end{equation}$$

Since the [Clebsch–Gordan coefficients](@ref) can be rewritten using 3j
symbols and the 3j symbols vanish unless $m_c + \alpha - m_3 = m_d -
\alpha - m_b = 0$, we have

$$\begin{equation}
\alpha = m_a - m_c = m_d-m_b
\implies
-\alpha + m_a + m_b = m_b + m_c.
\end{equation}$$

```@docs
TensorScalarProduct
dot(::Tensor, ::Tensor)
AngularMomentumAlgebra.integrate_spinors((a,b), X::TensorScalarProduct, (c,d))
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

```@meta
CurrentModule = nothing
DocTestSetup = nothing
```
