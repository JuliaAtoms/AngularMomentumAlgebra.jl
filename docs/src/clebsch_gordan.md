# Clebsch–Gordan coefficients

The Clebsch–Gordan coefficients are related to the 3j symbols as

$$\begin{equation}
C_{j_1m_1j_2m_2}^{j_3m_3} \equiv
\braket{j_1m_1j_2m_2}{j_3m_3} =
(-)^{j_1-j_2+m_3}\angroot{j_3}
\begin{pmatrix}
j_1&j_2&j_3\\
m_1&m_2&-m_3
\end{pmatrix}.
\tag{V8.1.12}
\end{equation}$$

```@meta
DocTestSetup = quote
    using AngularMomentumAlgebra
    using AtomicLevels
end
```

```@docs
clebsch_gordan_condon_shortley
```

```@meta
 DocTestSetup = nothing
```
