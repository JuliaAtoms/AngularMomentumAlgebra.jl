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

#### Coupled orbitals

In the case of coupled orbitals, `RelativisticOrbital`s in the
nomenclature of AtomicLevels.jl, if the operator acts on the entire
system (or at least the total angular momentum), the Wigner–Eckart
theorem ``\eqref{eqn:wigner-eckart}`` can be applied. If however, the
operators acts on a subsystem, the uncoupling formula
``\eqref{eqn:uncoupling}`` has to be employed.

```@docs
LinearAlgebra.dot(a::SpinOrbital{<:RelativisticOrbital}, 𝐓ᵏq::TensorComponent, b::SpinOrbital{<:RelativisticOrbital})
matrix_element(::Union{FullSystem,TotalAngularMomentumSubSystem}, a::SpinOrbital{<:RelativisticOrbital}, 𝐓ᵏq::TensorComponent, b::SpinOrbital{<:RelativisticOrbital})
matrix_element(system, a::SpinOrbital{<:RelativisticOrbital}, 𝐓ᵏq::TensorComponent, b::SpinOrbital{<:RelativisticOrbital})
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

```math
\begin{equation}
\tag{V13.1.39}
\matrixel{n_1'j_1'm_1';n_2'j_2'm_2'}{\tensor{T}^{(k)}_q(1)}{n_1j_1m_1;n_2j_2m_2}
=
\delta_{n_2'n_2}\delta_{j_2'j_2}\delta_{m_2'm_2}
\matrixel{n_1'j_1'm_1'}{\tensor{T}^{(k)}_q(1)}{n_1j_1m_1}
\end{equation}
```

#### Product tensors

##### Same coordinate

The matrix element of a scalar product of two tensors acting on
the same coordinate (of a subsystem) is given by

```math
\begin{equation}
\begin{aligned}
&\matrixel{n_1'j_1'n_2'j_2'j'm'}{[\tensor{P}^{(k)}(1)\cdot\tensor{Q}^{(k)}(1)]}{n_1j_1n_2j_2jm}\\
=&\delta_{n_2'n_2}\delta_{j_2'j_2}\delta_{j_1'j_1}\delta_{j'j}\delta_{m'm}
\frac{1}{\angroot{j_1}^2}
(-)^{-j_1}\\
&\times
\sum_{JN}(-)^J
\redmatrixel{n_1'j_1}{\tensor{P}^{(k)}(1)}{NJ}
\redmatrixel{NJ}{\tensor{Q}^{(k)}(1)}{n_1j_1}.
\end{aligned}
\tag{V13.1.43}
\end{equation}
```

Apart from the additional factor ``\delta_{n_2'n_2}\delta_{j_2'j_2}``,
this expression is equivalent to
``\eqref{eqn:scalar-product-tensor-matrix-element}``.

##### Different coordinates

```@docs
matrix_element((γj₁′, m₁′), (γj₂′, m₂′), X::TensorScalarProduct, (γj₁, m₁), (γj₂, m₂))
```

In the coupled basis, the equivalent formula is

```math
\begin{equation}
\begin{aligned}
&\matrixel{n_1'j_1'n_2'j_2'j'm'}{[\tensor{P}^{(k)}(1)\cdot\tensor{Q}^{(k)}(2)]}{n_1j_1n_2j_2jm}\\
=& \delta_{j'j}\delta_{m'm}
(-)^{j+j_1+j_2'}
\wignersixj{j_1'&j_1&k\\j_2&j_2'&j}\\
&\times
\redmatrixel{n_1'j_1'}{\tensor{P}^{(k)}(1)}{n_1j_1}
\redmatrixel{n_2'j_2'}{\tensor{Q}^{(k)}(2)}{n_2j_2}.
\end{aligned}
\tag{V13.1.29}
\end{equation}
```

### Old stuff

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
