# DSL for specifying new tensors

When specifying a new tensor type, a few functions need to be
provided:

- `AngularMomentumAlgebra.system(::TensorType)` which returns the
  [`AngularMomentumAlgebra.System`](@ref) the tensor acts on,
- `Base.iszero` which indicates whether a reduced matrix element of
  the tensor vanishes without actually computing it,
- `couplings` which generates all sets of quantum numbers the tensor
  couples to from a given set of quantum numbers,
- `rme` which computes the actual reduced matrix element.

```@docs
AngularMomentumAlgebra.@tensor
AngularMomentumAlgebra.strip′
AngularMomentumAlgebra.identify_quantum_numbers
AngularMomentumAlgebra.recursepm
AngularMomentumAlgebra.parse_selection_rule
AngularMomentumAlgebra.generate_iszero
AngularMomentumAlgebra.generate_rme
AngularMomentumAlgebra.generate_couplings
```
