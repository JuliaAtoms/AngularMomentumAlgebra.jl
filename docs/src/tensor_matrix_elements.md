# Tensor matrix elements

The matrix element of any irreducible tensor operator of rank ``k`` is
defined as

```math
\begin{equation}
\tag{V13.1.1}
\matrixel{n'j'm'}{\tensor{T}^{(k)}_q}{njm} \defd
\int\diff{\tau}
\conj{\Psi}_{n'j'm'}
\tensor{T}^{(k)}_q
\Psi_{njm}.
\end{equation}
```

For simple quantum systems, it is most easily evaluated using [the
Wigner–Eckart theorem](@ref wigner_eckart); however, in most cases,
the quantum system and/or the tensor is composed of multiple parts,
which complicates the situation. The table below shows the most common
cases:

| Basis \ Tensor | Acts on entire system                                | Acts on subsystems                     |
|----------------|------------------------------------------------------|----------------------------------------|
| **Uncoupled**  | Transform to coupled basis: ``\eqref{eqn:coupling}`` | [Direct evaluation](@ref)              |
| **Coupled**    | Wigner–Eckart: ``\eqref{eqn:wigner-eckart}``         | Uncoupling: ``\eqref{eqn:uncoupling}`` |

## Interface

There are two interfaces provided for computation of matrix element
of tensor operators:

- A low-level interface `matrix_element((γj′, m′), 𝐓ᵏq, (γj′, m′))`
  (and friends), where `γj′`, `m′`, `γj′`, and `m′` are quantum
  numbers.
- A high-level interface `dot(a, 𝐓ᵏq, b)`, where `a` and `b` are
  `SpinOrbital`s from
  [AtomicLevels.jl](https://github.com/JuliaAtoms/AtomicLevels.jl.git). The
  high-level interface dispatches as appropriate to the low-level
  interface, depending on whether `a` and `b` are expressed in the
  coupled basis or not, and which part of the quantum system `𝐓ᵏq`
  acts on.

### High-level interface

There are two main functions in the high-level interface:
- `dot(a, 𝐓ᵏq, b)` for one-body interactions
- `dot((a,b), 𝐓ᵏq, (c,d))` for two-body interactions (e.g. [Coulomb
  interaction](@ref)).

```@docs
LinearAlgebra.dot(a::SpinOrbital, 𝐓ᵏq::TensorComponent, b::SpinOrbital)
```

#### Intermediate-level interface

`dot` dispatches to this level, passing the [`system`](@ref) of the
tensor operator considered as the first argument. At this level, the
spin-orbitals are translated into quantum numbers, employed by the
low-level interface.

##### Coupled orbitals

In the case of coupled orbitals, `RelativisticOrbital`s in the
nomenclature of AtomicLevels.jl, if the operator acts on the entire
system (or at least the total angular momentum), the Wigner–Eckart
theorem ``\eqref{eqn:wigner-eckart}`` can be applied. If however, the
operators acts on a subsystem, the uncoupling formula
``\eqref{eqn:uncoupling}`` has to be employed.

```@docs
matrix_element(::Union{FullSystem,TotalAngularMomentumSubSystem}, a::SpinOrbital{<:RelativisticOrbital}, 𝐓ᵏq::TensorComponent, b::SpinOrbital{<:RelativisticOrbital})
matrix_element(system, a::SpinOrbital{<:RelativisticOrbital}, 𝐓ᵏq::TensorComponent, b::SpinOrbital{<:RelativisticOrbital})
matrix_element(::Tuple{S,S}, a::SpinOrbital{<:RelativisticOrbital}, X::TensorScalarProduct, b::SpinOrbital{<:RelativisticOrbital}) where {S<:Union{FullSystem,TotalAngularMomentumSubSystem}}
matrix_element(systems::Tuple{S,S}, a::SpinOrbital{<:RelativisticOrbital}, X::TensorScalarProduct, b::SpinOrbital{<:RelativisticOrbital}) where {S<:SubSystem}
matrix_element(systems::Tuple{<:SubSystem,<:SubSystem}, a::SpinOrbital{<:RelativisticOrbital}, X::TensorScalarProduct, b::SpinOrbital{<:RelativisticOrbital})
```

- [x] ``\matrixel{a}{\tensor{L}^2}{b}``, ``\matrixel{a}{\tensor{S}^2}{b}``
  Scalar product tensors acting on one coordinate only
- [x] ``\matrixel{a}{\tensor{L}\cdot\tensor{S}}{b}`` Scalar product tensor acting on two
  different coordinates of the same orbital
- [ ] ``\matrixel{ab}{\tensor{C}^{(k)}\cdot\tensor{C}^{(k)}}{cd}`` Scalar
  product tensor acting on (subsystems of) two different orbitals

##### Uncoupled orbitals

```@docs
matrix_element(system::Union{FullSystem,TotalAngularMomentumSubSystem}, a::SpinOrbital{<:Orbital}, 𝐓ᵏq::TensorComponent, b::SpinOrbital{<:Orbital})
matrix_element(system, a::SpinOrbital{<:Orbital}, 𝐓ᵏq::TensorComponent, b::SpinOrbital{<:Orbital})
matrix_element(systems::Tuple{S,S}, a::SpinOrbital{<:Orbital}, X::TensorScalarProduct, b::SpinOrbital{<:Orbital}) where {S<:Union{FullSystem,TotalAngularMomentumSubSystem}}
matrix_element(systems::Tuple{S,S}, a::SpinOrbital{<:Orbital}, X::TensorScalarProduct, b::SpinOrbital{<:Orbital}) where {S<:SubSystem}
matrix_element(systems::Tuple{<:SubSystem,<:SubSystem}, a::SpinOrbital{<:Orbital}, X::TensorScalarProduct, b::SpinOrbital{<:Orbital})
```

## Tensor acts on entire system

### [The Wigner–Eckart theorem](@id wigner_eckart)

```@docs
matrix_element((γj′, m′), Tᵏq::TensorComponent, (γj, m))
```

### Product tensors

```@docs
matrix_element((γj′, m′), X::TensorScalarProduct, (γj, m))
```

### Uncoupled basis functions

For a tensor operator that depends on all coordinates, its matrix
element in the uncoupled basis are computed via a basis transform to
the coupled basis:

```@docs
matrix_element((γj₁′, m₁′), (γj₂′, m₂′), 𝐓ᵏq, (γj₁, m₁), (γj₂, m₂))
```

## Tensor acts on subsystems

### Uncoupling

When the tensor operator is reducible and only acts on one part of the
quantum system, in the coupled basis we employ the following
uncoupling formula:

```@docs
matrix_element_via_uncoupling
```

### Direct evaluation

### Product tensors

```@docs
matrix_element2((γj₁′, m₁′), (γj₂′, m₂′), X::TensorScalarProduct, (γj₁, m₁), (γj₂, m₂))
matrix_element2((γj₁′, γj₂′, j′, m′), X::TensorScalarProduct, (γj₁, γj₂, j, m))
```

## Old stuff

```@docs
AngularMomentumAlgebra.complementary_space_factor
AngularMomentumAlgebra.rme_j′j
```

The matrix element of a scalar product of two tensors depends on
whether the scalar product is an irreducible tensor or not.


## Reference

```@docs
dot(a, X::TensorScalarProduct, b)
dot((a,b)::Tuple, X::TensorScalarProduct, (c,d)::Tuple)
```
