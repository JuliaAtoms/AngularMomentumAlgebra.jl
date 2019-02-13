# Multipole expansions

```@meta
CurrentModule = AngularMomentumAlgebra
DocTestSetup = quote
    using AngularMomentumAlgebra
    using AtomicLevels
end
```

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
multipole_expand_scalar_product
multipole_expand
```

## Spherical tensors
For spherical tensors, we can insert the expression for the reduced
matrix elements (V13.2.107) into (V13.1.26) above, to get

$$\begin{equation}
\begin{aligned}
\implies
&\matrixel{\ell_am_a;\ell_bm_b}{[\tensor{C}^{(k)}(1)\cdot\tensor{C}^{(k)}(2)]}{\ell_cm_c;\ell_dm_d}\\
=&
\frac{1}{\angroot{\ell_a\ell_b}}
\sum_\alpha(-)^{-\alpha}
C_{\ell_cm_c;k,\alpha}^{\ell_am_a}
C_{\ell_dm_d;k,-\alpha}^{\ell_bm_b}
\redmatrixel{\ell_a}{\tensor{C}^{(k)}(1)}{\ell_c}
\redmatrixel{\ell_b}{\tensor{C}^{(k)}(2)}{\ell_d}\\
=&
\frac{\angroot{\ell_c\ell_d}}{\angroot{\ell_a\ell_b}}
C_{\ell_c0;k0}^{\ell_a0}
C_{\ell_d0;k0}^{\ell_b0}
\sum_\alpha(-)^{-\alpha}
C_{\ell_cm_c;k,\alpha}^{\ell_am_a}
C_{\ell_dm_d;k,-\alpha}^{\ell_bm_b}\\
=&\begin{aligned}[t]
\sum_\alpha&
(-)^{-\alpha + m_a + m_b}
\angroot{\ell_a\ell_b\ell_c\ell_d}\\
&\times
\begin{pmatrix}
\ell_c&k&\ell_a\\
m_c&\alpha&-m_a
\end{pmatrix}
\begin{pmatrix}
\ell_d&k&\ell_b\\
m_d&-\alpha&-m_b
\end{pmatrix}\\
&\times
\begin{pmatrix}
\ell_c&k&\ell_a\\0&0&0
\end{pmatrix}
\begin{pmatrix}
\ell_d&k&\ell_b\\0&0&0
\end{pmatrix}.
\end{aligned}
\end{aligned}
\end{equation}$$

```@meta
 DocTestSetup = nothing
```
