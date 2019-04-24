# Common routines

```@meta
CurrentModule = AngularMomentumAlgebra
DocTestSetup = quote
    using AngularMomentumAlgebra
    using AtomicLevels
end
```

A commonly occurring factor in angular momentum algebra is

$$\begin{equation}
\angroot{j_1j_2...j_n}
\defd[(2j_1+1)(2j_2+1)...(2j_n+1)]^{1/2}.
\tag{V13.1.3½}
\end{equation}$$

```@docs
∏
triangle_range
powneg1
jmⱼ
spin
LinearCombination
@linearly_combinable
```

```@meta
 DocTestSetup = nothing
```
