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
system(::Type{<:Tensor})
AngularMomentumAlgebra.OrbitalRadialOverlap
TensorOperator
many_electron_scalar_product
```

## Cartesian tensor components

The transformation of the tensor components of a rank-1 tensor from
the Cartesian basis to the "natural" basis is given by:

```math
\begin{equation}
M(+1,0,-1 \leftarrow x, y, z)=
\begin{bmatrix}
-\frac{1}{\sqrt{2}}&\frac{\im}{\sqrt{2}}&\cdot\\
\cdot&\cdot&1\\
\frac{1}{\sqrt{2}}&\frac{\im}{\sqrt{2}}&\cdot
\end{bmatrix}
\end{equation}
```

```@docs
cartesian_tensor_component
```

## Product tensors

```@docs
TensorProduct
```

### Scalar product

```@docs
TensorScalarProduct
dot(::Tensor, ::Tensor)
```

```@meta
CurrentModule = nothing
DocTestSetup = nothing
```
