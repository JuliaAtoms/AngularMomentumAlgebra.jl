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

##### Different coordinates

The matrix element of a scalar product of two tensors acting on
different coordinates is given by (in the uncoupled basis)

```math
\begin{equation}
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
\end{equation}
```

Since the [Clebsch–Gordan coefficients](@ref) can be rewritten using 3j
symbols and the 3j symbols vanish unless $m_c + \alpha - m_3 = m_d -
\alpha - m_b = 0$, we have

```math
\alpha = m_a - m_c = m_d-m_b
\implies
-\alpha + m_a + m_b = m_b + m_c.
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

This case occurs in two-body interactions, such as the [Coulomb
interaction](@ref), where ``a,b`` and ``c,d`` are pairs of orbitals
and the scalar product tensor is a multipole expansion in [Spherical
tensors](@ref tensors_spherical_tensors), but also in the case of the
operator ``\tensor{L}\cdot\tensor{S}`` and coordinates ``1`` and ``2``
correspond to orbital and spin angular momenta, respectively.

We can verify this using the classical result known from spin–orbit
splitting:

```math
\begin{aligned}
J^2 &= (\tensor{L}+\tensor{S})^2 = L^2 + 2\tensor{L}\cdot\tensor{S} + S^2\\
\implies
\expect{\tensor{L}\cdot\tensor{S}} &=
\frac{1}{2}(\expect{J^2} - \expect{L^2} - \expect{S^2}) =
\frac{1}{2}[J(J+1) - L(L+1) - S(S+1)]
\end{aligned}
```

In the uncoupled basis, ``J`` is not a good quantum number (it is not
a constant of motion), except for _pure states_, i.e. those with
maximal ``\abs{m_\ell + m_s}``:

```jldoctest
julia> using AngularMomentumAlgebra, AtomicLevels, HalfIntegers

julia> a = SpinOrbital(o"2p", 1, half(1))
2p₁α

julia> 𝐋 = OrbitalAngularMomentum()
𝐋̂⁽¹⁾

julia> 𝐒 = SpinAngularMomentum()
𝐒̂⁽¹⁾

julia> X = 𝐋⋅𝐒
(𝐋̂⁽¹⁾⋅𝐒̂⁽¹⁾)

julia> J = half(3) # Set manually since in uncoupled basis
3/2

julia> L = a.orb.ℓ
1

julia> S = half(1)
1/2

julia> dot(a, X, a), 1/2*(J*(J+1)-L*(L+1)-S*(S+1))
(0.5, 0.5)
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
