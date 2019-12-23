# * Tensor acts on entire system
# ** Wigner–Eckart

@doc raw"""
    matrix_element((γj′, m′), Tᵏq::TensorComponent, (γj, m))

Calculate the matrix element `⟨γ′j′m′|Tᵏq|γjm⟩` via Wigner–Eckart's
theorem:


```math
\begin{equation}
\begin{aligned}
\matrixel{\gamma'j'm'}{\tensor{T}^{(k)}_q}{\gamma jm}
&\defd
(-)^{2k} \frac{1}{\angroot{j'}}
C_{jm;kq}^{j'm'}
\redmatrixel{\gamma' j'}{\tensor{T}^{(k)}}{\gamma j} \\
&=
(-)^{j'-m'}
\begin{pmatrix}
j'&k&j\\
-m'&q&m
\end{pmatrix}
\redmatrixel{\gamma'j'}{\tensor{T}^{(k)}}{\gamma j},
\end{aligned}
\label{eqn:wigner-eckart}
\tag{V13.1.2}
\end{equation}
```

where the _reduced matrix element_
$\redmatrixel{n'j'}{\tensor{T}^{(k)}}{nj}$ does not depend on
$m,m'$. `j′` and `j` are the total angular momenta with `m′` and `m`
being their respective projections. `γ′` and `γ` are all other quantum
numbers needed to fully specify the quantum system; their presence
depend on the quantum system.

# Examples

```jldoctest
julia> matrix_element((2, 1), TensorComponent(OrbitalAngularMomentum(), 1), (2, 0))
-1.7320508075688774

julia> matrix_element((0, 0), TensorComponent(SphericalTensor(1), 0), (1, 0))
0.5773502691896256

julia> matrix_element(((1,half(1),half(1)), -half(1)),
                     TensorComponent(TotalAngularMomentum(), -1),
                     ((1,half(1),half(1)), half(1)))
0.7071067811865475
```
"""
function matrix_element((γj′, m′), Tᵏq::TensorComponent, (γj, m))
    Tᵏ = parent(Tᵏq)
    r = rme(γj′, Tᵏ, γj)
    iszero(r) && return 0
    j′,j = last(γj′), last(γj)
    c = powneg1(Int(j′-m′))*wigner3j(j′, rank(Tᵏ), j,
                                     -m′, component(Tᵏq), m)
    iszero(c) && return 0
    c*r
end

# ** Product tensors

@doc raw"""
    matrix_element((γj′, m′), X::TensorScalarProduct, (γj, m))

Calculate the matrix element of a scalar product tensor according to:

```math
\begin{equation}
\begin{aligned}
\matrixel{n'j'm'}{[\tensor{P}^{(k)}\cdot\tensor{Q}^{(k)}]}{njm}
=&
\delta_{jj'}\delta_{mm'}
\frac{1}{\angroot{j}^2}\\
&\times\sum_{n_1j_1}
(-)^{-j+j_1}
\redmatrixel{n'j}{\tensor{P}^{(k)}}{n_1j_1}
\redmatrixel{n_1j_1}{\tensor{Q}^{(k)}}{nj}
\end{aligned}
\label{eqn:scalar-product-tensor-matrix-element}
\tag{V13.1.11}
\end{equation}
```

The permissible values of ``n_1j_1`` in the summation are found using
[`AngularMomentumAlgebra.couplings`](@ref); it is assumed that the
summation only consists of a finite amount of terms and that

```math
\redmatrixel{n'j}{\tensor{P}^{(k)}}{n_1j_1}\neq0
\iff
\redmatrixel{n_1j_1}{\tensor{P}^{(k)}}{n'j}\neq0,
```

i.e. that ``\tensor{P}^{(k)}`` is (skew)symmetric.

# Examples

```jldoctest
julia> 𝐒 = SpinAngularMomentum()
𝐒̂⁽¹⁾

julia> 𝐒² = 𝐒⋅𝐒
(𝐒̂⁽¹⁾⋅𝐒̂⁽¹⁾)

julia> matrix_element((half(1), half(1)),
                      𝐒², (half(1), half(1)))
0.7499999999999998

julia> half(1)*(half(1)+1) # S(S+1)
0.75

julia> 𝐉 = TotalAngularMomentum()
𝐉̂⁽¹⁾

julia> 𝐉² = 𝐉⋅𝐉
(𝐉̂⁽¹⁾⋅𝐉̂⁽¹⁾)

julia> matrix_element(((1, half(1), half(3)), half(3)),
                      𝐉², ((1, half(1), half(3)), half(3)))
3.7500000000000004

julia> half(3)*(half(3)+1) # J(J+1)
3.75
```
"""
function matrix_element((γj′, m′), X::TensorScalarProduct, (γj, m))
    j′, j = last(γj′), last(γj)
    @δ j′,j m′,m

    T,U = X.T,X.U
    c = 0.0

    # We assume that iszero(⟨n′j′||𝐓̂⁽ᵏ⁾||n₁j₁⟩) ⇔
    # iszero(⟨n₁j₁||𝐓̂⁽ᵏ⁾||n′j′⟩).
    Tγj₁ = couplings(T, γj)
    Uγj₁ = couplings(U, γj)
    γj₁s = map(((Tγj₁, Uγj₁),) -> ∩(Tγj₁, Uγj₁), zip(Tγj₁, Uγj₁))
    for γj₁ ∈ Iterators.product(γj₁s...)
        length(γj₁) == 1 && (γj₁ = first(γj₁))
        j₁ = last(γj₁)
        c += powneg1(Int(-j+j₁))*rme(γj′,T,γj₁)*rme(γj₁,U,γj)
    end

    c/∏(j)^2
end

# ** Uncoupled basis functions

@doc raw"""
    matrix_element((γj₁′, m₁′), (γj₂′, m₂′), 𝐓ᵏq, (γj₁, m₁), (γj₂, m₂))

Compute the matrix element of the irreducible tensor `𝐓ᵏq` acting on
coordinates `1` and `2`, by first coupling `γj₁′m₁′γj₂′m₂′` and
`γj₁m₁γj₂m₂` to all permissible `j′m′` and `jm`, respectively,
according to

```math
\begin{equation}
\begin{aligned}
&\matrixel{γ_1'j_1'm_1';γ_2'j_2'm_2'}{\tensor{P}^{(k)}_q(1,2)}{γ_1j_1m_1;γ_2j_2m_2} \\
=& (-)^{2k}
\frac{1}{\angroot{j'}}
\sum_{jmj'm'}
C_{j_1m_1;j_2m_2}^{jm}
C_{j_1'm_1';j_2'm_2'}^{j'm'}
C_{jm;kq}^{j'm'}\\
&\times
\redmatrixel{γ_1'j_1'γ_2'j_2'j'}{\tensor{P}^{(k)}(1,2)}{γ_1j_1γ_2j_2j} \\
\equiv&
\sum_{jmj'm'}
C_{j_1m_1;j_2m_2}^{jm}
C_{j_1'm_1';j_2'm_2'}^{j'm'}
\matrixel{γ_1'j_1'γ_2'j_2'j'm'}{\tensor{P}^{(k)}(1,2)}{γ_1j_1γ_2j_2jm}
\end{aligned}
\tag{V13.1.23}
\label{eqn:coupling}
\end{equation}
```

The non-vanishing terms of the sum are found using
[`AngularMomentumAlgebra.couplings`](@ref).

# Examples

```jldoctest
julia> 𝐉 = TotalAngularMomentum()
𝐉̂⁽¹⁾

julia> 𝐉₀ = TensorComponent(𝐉, 0)
𝐉̂⁽¹⁾₀

julia> matrix_element((1,1), (half(1),half(1)),
                      𝐉₀, (1,1), (half(1), half(1)))
1.5

julia> matrix_element((1,-1), (half(1),half(1)),
                      𝐉₀, (1,-1), (half(1), half(1)))
-0.4999999999999999

julia> 𝐉₁ = TensorComponent(𝐉, 1)
𝐉̂⁽¹⁾₁

julia> matrix_element((1,1), (half(1),half(1)),
                      𝐉₁, (1,0), (half(1), half(1)))
-1.0

julia> 𝐉² = 𝐉⋅𝐉
(𝐉̂⁽¹⁾⋅𝐉̂⁽¹⁾)

julia> matrix_element((1,1), (half(1),half(1)),
                      𝐉², (1,1), (half(1), half(1)))
3.7500000000000004

julia> half(3)*(half(3)+1) # J(J+1)
3.75
```
"""
function matrix_element((γj₁′, m₁′), (γj₂′, m₂′), 𝐓ᵏq, (γj₁, m₁), (γj₂, m₂))
    j₁′,j₂′,j₁,j₂ = last(γj₁′), last(γj₂′), last(γj₁), last(γj₂)
    j′s,m′ = couplings(j₁′, m₁′, j₂′, m₂′)
    js,m = couplings(j₁, m₁, j₂, m₂)

    γj′ = (γj₁′..., γj₂′...)
    γj = (γj₁..., γj₂...)

    v = 0.0
    for j′ ∈ j′s
        c′ = clebschgordan(j₁′, m₁′, j₂′, m₂′, j′, m′)
        for j ∈ js
            me = matrix_element(((γj′...,j′), m′), 𝐓ᵏq, ((γj...,j), m))
            iszero(me) && continue
            c = c′*clebschgordan(j₁, m₁, j₂, m₂, j, m)
            v += c*me
        end
    end
    v
end

# * Tensor acts on subsystems

# ** Uncoupling of coupled basis states

@doc raw"""
    matrix_element_via_uncoupling((γj₁′, γj₂′, j′, m′), 𝐓ᵏq, (γj₁, γj₂, j, m))

Compute the matrix element of the tensor `𝐓ᵏq` which acts on
coordinate `1` only in the coupled basis, by employing the uncoupling
formula

```math
\begin{equation}
\begin{aligned}
&\matrixel{γ_1'j_1'γ_2'j_2'j'm'}{\tensor{T}^{(k)}_q(1)}{γ_1j_1γ_2j_2jm}\\
=& \delta_{j_2'j_2}\delta_{γ_2'γ_2}
(-)^{j+j_1'+j_2-k}
\angroot{j}
C_{jm;kq}^{j'm'}\\
&\times
\wignersixj{j_1&j_2&j\\j'&k&j_1'}
\redmatrixel{γ_1'j_1'}{\tensor{T}^{(k)}(1)}{γ_1j_1}
\end{aligned}
\tag{V13.1.40}
\label{eqn:uncoupling}
\end{equation}
```

# Examples

```jldoctest
julia> 𝐋₀ = TensorComponent(OrbitalAngularMomentum(), 0)
𝐋̂⁽¹⁾₀

julia> matrix_element_via_uncoupling((1, half(1), half(3), half(3)),
                                     𝐋₀, (1, half(1), half(3), half(3)))
0.9999999999999999

julia> matrix_element((1, 1), 𝐋₀, (1, 1)) # For comparison
0.9999999999999999

julia> 𝐒₀ = TensorComponent(SpinAngularMomentum(), 0)
𝐒̂⁽¹⁾₀

julia> matrix_element_via_uncoupling((half(1), 1, half(3), half(3)),
                                     𝐒₀, (half(1), 1, half(3), half(3)))
0.49999999999999994

julia> matrix_element((half(1),half(1)), 𝐒₀, (half(1),half(1)))
0.49999999999999994
```

"""
function matrix_element_via_uncoupling((γj₁′, γj₂′, j′, m′), 𝐓ᵏq, (γj₁, γj₂, j, m))
    @δ γj₂′,γj₂

    𝐓ᵏ = parent(𝐓ᵏq)
    r = rme(γj₁′, 𝐓ᵏ, γj₁)
    iszero(r) && return 0

    j₁′,j₂′,j₁,j₂ = last(γj₁′), last(γj₂′), last(γj₁), last(γj₂)
    k = rank(𝐓ᵏ)
    q = component(𝐓ᵏq)
    c = clebschgordan(j,m,k,q,j′,m′)
    iszero(c) && return 0
    w6j = wigner6j(j₁, j₂, j,
                   j′, k,  j₁′)
    iszero(w6j) && return 0

    powneg1(Int(j+j₁′+j₂-k))*∏(j)*c*w6j*r
end


# ** Uncoupled basis states

@doc raw"""
    matrix_element2((γj₁′, m₁′), (γj₂′, m₂′), X::TensorScalarProduct, (γj₁, m₁), (γj₂, m₂))

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
\redmatrixel{n_bj_b}{\tensor{Q}^{(k)}(2)}{n_dj_d} \\
\equiv&
\sum_\alpha
(-)^{-\alpha}
\matrixel{n_aj_am_a}{\tensor{P}^{(k)}_{\alpha}(1)}{n_cj_cm_c}
\matrixel{n_bj_bm_b}{\tensor{Q}^{(k)}_{-\alpha}(2)}{n_dj_dm_d}
\end{aligned}
\tag{V13.1.26}
\label{eqn:scalar-product-tensor-matrix-element-diff-coords-uncoupled}
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

This case occurs in two-body interactions, such as the [Coulomb
interaction](@ref), where ``a,b`` and ``c,d`` are pairs of orbitals
and the scalar product tensor is a term in the multipole expansion in
terms of [Spherical tensors](@ref tensors_spherical_tensors):

```jldoctest
julia> 𝐂⁰ = SphericalTensor(0)
𝐂̂⁽⁰⁾

julia> matrix_element2((0, 0), (0, 0), 𝐂⁰⋅𝐂⁰, (0,0), (0, 0)) # ⟨1s₀,1s₀|𝐂⁰⋅𝐂⁰|1s₀,1s₀⟩
1.0

julia> 𝐂¹ = SphericalTensor(1)
𝐂̂⁽¹⁾

julia> matrix_element2((0, 0), (1, 0), 𝐂¹⋅𝐂¹, (1,0), (2, 0)) # ⟨1s₀,2p₀|𝐂¹⋅𝐂¹|2p₀,3d₀⟩
0.29814239699997186

julia> matrix_element2((0, 0), (1, 1), 𝐂¹⋅𝐂¹, (1,0), (2, 1)) # ⟨1s₀,2p₁|𝐂¹⋅𝐂¹|2p₀,3d₁⟩
0.25819888974716104
```

but also in the case of the operator ``\tensor{L}\cdot\tensor{S}`` and
coordinates ``1`` and ``2`` correspond to orbital and spin angular
momenta, respectively. We can verify this using the classical result
known from spin–orbit splitting:

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
julia> 𝐋 = OrbitalAngularMomentum()
𝐋̂⁽¹⁾

julia> 𝐒 = SpinAngularMomentum()
𝐒̂⁽¹⁾

julia> X = 𝐋⋅𝐒
(𝐋̂⁽¹⁾⋅𝐒̂⁽¹⁾)

julia> matrix_element2((1, 1), (half(1), half(1)),
                      X, (1,1), (half(1), half(1)))
0.4999999999999999

julia> 1/2*(half(3)*(half(3)+1)-1*(1+1)-half(1)*(half(1)+1)) # 1/2(J(J+1)-L(L+1)-S(S+1))
0.5
```
"""
function matrix_element2((γj₁′, m₁′), (γj₂′, m₂′), X::TensorScalarProduct, (γj₁, m₁), (γj₂, m₂))
    T,U = X.T,X.U
    k = rank(T)

    α = Int(m₁′-m₁)
    (α != m₂-m₂′ || abs(α) > k) && return 0

    powneg1(-α)*
    matrix_element((γj₁′, m₁′), TensorComponent(T,α), (γj₁, m₁))*
    matrix_element((γj₂′, m₂′), TensorComponent(U,-α), (γj₂, m₂))
end

@doc raw"""
    matrix_element2((γj₁′, γj₂′, j′, m′), X::TensorScalarProduct, (γj₁, γj₂, j, m))

The matrix element of a scalar product of two tensors acting on
different coordinates is given by (in the coupled basis)

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
\label{eqn:scalar-product-tensor-matrix-element-diff-coords-coupled}
\end{equation}
```
"""
function matrix_element2((γj₁′, γj₂′, j′, m′), X::TensorScalarProduct, (γj₁, γj₂, j, m))
    @δ j′,j m′,m
    j₁′,j₂′,j₁,j₂ = last(γj₁′), last(γj₂′), last(γj₁), last(γj₂)

    T,U = X.T,X.U
    Tr = rme(γj₁′, T, γj₁)
    iszero(Tr) && return 0
    Ur = rme(γj₂′, U, γj₂)
    iszero(Ur) && return 0

    powneg1(Int(j+j₁+j₂′))*Tr*Ur*wigner6j(j₁′, j₁, rank(T),
                                          j₂, j₂′, j)
end

# * Tensor matrix elements in orbital basis

"""
    dot(a::SpinOrbital,
        𝐓ᵏq::Union{TensorComponent,TensorScalarProduct},
        b::SpinOrbital)

Compute the matrix element `⟨a|𝐓ᵏq|b⟩` in the basis of spin-orbitals,
dispatching to the correct low-level function `matrix_element`,
depending on the value of `system(𝐓ᵏq)`.

# Examples with coupled orbitals

```jldoctest
julia> a,b,c,d = (SpinOrbital(ro"2p", half(3)),
                  SpinOrbital(ro"2p", half(1)),
                  SpinOrbital(ro"2s", half(1)),
                  SpinOrbital(ro"3d", half(3)))
(2p(3/2), 2p(1/2), 2s(1/2), 3d(3/2))

julia> 𝐋, 𝐒, 𝐉 = OrbitalAngularMomentum(), SpinAngularMomentum(), TotalAngularMomentum()
(𝐋̂⁽¹⁾, 𝐒̂⁽¹⁾, 𝐉̂⁽¹⁾)

julia> 𝐋², 𝐒², 𝐉² = 𝐋⋅𝐋, 𝐒⋅𝐒, 𝐉⋅𝐉
((𝐋̂⁽¹⁾⋅𝐋̂⁽¹⁾), (𝐒̂⁽¹⁾⋅𝐒̂⁽¹⁾), (𝐉̂⁽¹⁾⋅𝐉̂⁽¹⁾))

julia> dot(a, cartesian_tensor_component(𝐉, :x), b),
           1/2*√((3/2+1/2+1)*(3/2-1/2)) # 1/2√((J+M+1)*(J-M))
(0.8660254037844386, 0.8660254037844386)

julia> dot(a, cartesian_tensor_component(𝐉, :z), a), a.m[1]
(1.5, 3/2)

julia> dot(a, TensorComponent(𝐋, 0), a)
0.9999999999999999

julia> dot(c, cartesian_tensor_component(Gradient(), :x), a)
- 0.408248(∂ᵣ + 2/r)

julia> dot(c, cartesian_tensor_component(SphericalTensor(1), :x), a)
-0.40824829046386296

julia> orbitals = rsos"2[p]"
6-element Array{SpinOrbital{RelativisticOrbital{Int64},Tuple{Half{Int64}}},1}:
 2p-(-1/2)
 2p-(1/2)
 2p(-3/2)
 2p(-1/2)
 2p(1/2)
 2p(3/2)

julia> map(o -> dot(o, 𝐉², o), orbitals)
6-element Array{Float64,1}:
 0.7499999999999998
 0.7499999999999998
 3.7500000000000004
 3.7500000000000004
 3.7500000000000004
 3.7500000000000004

julia> 1/2*(1/2+1),3/2*(3/2+1) # J(J+1)
(0.75, 3.75)

julia> dot(a, 𝐋², a), 1*(1+1)
(2.0, 2)

julia> dot(d, 𝐋², d), 2*(2+1)
(5.999999999999999, 6)

julia> dot(a, 𝐒², a), 1/2*(1/2+1)
(0.7499999999999998, 0.75)

julia> dot(a, 𝐋⋅𝐒, a)
0.4999999999999999
```

# Examples with uncoupled orbitals

```jldoctest
julia> a,b,c,d = (SpinOrbital(o"2p", 1, half(1)),
                  SpinOrbital(o"2p", -1, half(1)),
                  SpinOrbital(o"2s", 0, half(1)),
                  SpinOrbital(o"3d", 2, -half(1)))
(2p₁α, 2p₋₁α, 2s₀α, 3d₂β)

julia> 𝐋,𝐒,𝐉 = OrbitalAngularMomentum(),SpinAngularMomentum(),TotalAngularMomentum()
(𝐋̂⁽¹⁾, 𝐒̂⁽¹⁾, 𝐉̂⁽¹⁾)

julia> 𝐋²,𝐒²,𝐉² = 𝐋⋅𝐋,𝐒⋅𝐒,𝐉⋅𝐉
((𝐋̂⁽¹⁾⋅𝐋̂⁽¹⁾), (𝐒̂⁽¹⁾⋅𝐒̂⁽¹⁾), (𝐉̂⁽¹⁾⋅𝐉̂⁽¹⁾))

julia> dot(a, cartesian_tensor_component(𝐉, :z), a), sum(a.m)
(1.5, 3/2)

julia> dot(a, TensorComponent(𝐋, 0), a)
0.9999999999999999

julia> dot(d, TensorComponent(𝐋, 0), d)
1.9999999999999998

julia> dot(d, TensorComponent(𝐒, 0), d)
-0.49999999999999994

julia> dot(c, cartesian_tensor_component(Gradient(), :x), a)
- 0.408248(∂ᵣ + 2/r)

julia> dot(c, cartesian_tensor_component(SphericalTensor(1), :x), a)
-0.40824829046386285

julia> orbitals = sos"2[p]"
6-element Array{SpinOrbital{Orbital{Int64},Tuple{Int64,Half{Int64}}},1}:
 2p₋₁α
 2p₋₁β
 2p₀α
 2p₀β
 2p₁α
 2p₁β

julia> # Only 2p₋₁β and 2p₁α are pure states, with J = 3/2 => J(J + 1) = 3.75
       map(o -> dot(o, 𝐉², o), orbitals)
6-element Array{Float64,1}:
 1.7499999999999998
 3.7500000000000004
 2.75
 2.75
 3.7500000000000004
 1.7499999999999998

julia> dot(a, 𝐋², a), 1*(1+1)
(2.0, 2)

julia> dot(d, 𝐋², d), 2*(2+1)
(5.999999999999999, 6)

julia> dot(a, 𝐒², a), half(1)*(half(1)+1)
(0.7499999999999998, 0.75)

julia> dot(a, 𝐋⋅𝐒, a)
0.4999999999999999
```
"""
LinearAlgebra.dot(a::SpinOrbital,
                  𝐓ᵏq::Union{TensorComponent,TensorScalarProduct},
                  b::SpinOrbital) =
                      matrix_element(system(𝐓ᵏq), a, 𝐓ᵏq, b)

LinearAlgebra.dot(o′, lct::LinearCombinationTensor, o) =
    sum(filter!(!iszero, [c*dot(o′, Tᵏq, o) for (Tᵏq,c) in lct]))

# ** Coupled orbitals

"""
    matrix_element(::Union{FullSystem,TotalAngularMomentumSubSystem},
                   a::SpinOrbital{<:RelativisticOrbital},
                   𝐓ᵏq::TensorComponent,
                   b::SpinOrbital{<:RelativisticOrbital})

The matrix element of a tensor acting on the full system or the total
angular momentum, evaluated in the basis of coupled spin-orbitals, is
simply computed using the Wigner–Eckart theorem
``\\eqref{eqn:wigner-eckart}``.
"""
matrix_element(system::Union{FullSystem,TotalAngularMomentumSubSystem},
               a::SpinOrbital{<:RelativisticOrbital},
               𝐓ᵏq::TensorComponent,
               b::SpinOrbital{<:RelativisticOrbital}) =
                   matrix_element(quantum_numbers(system, a),
                                  𝐓ᵏq, quantum_numbers(system, b))

"""
    matrix_element(system,
                   a::SpinOrbital{<:RelativisticOrbital},
                   𝐓ᵏq::TensorComponent,
                   b::SpinOrbital{<:RelativisticOrbital})

The matrix element of a tensor acting on `system`, which is a
subsystem, evaluated in the basis coupled spin-orbitals, needs to be
computed via the uncoupling formula ``\\eqref{eqn:uncoupling}``.
"""
function matrix_element(system,
                        a::SpinOrbital{<:RelativisticOrbital},
                        𝐓ᵏq::TensorComponent,
                        b::SpinOrbital{<:RelativisticOrbital})
    j′, m′ = a.orb.j, a.m[1]
    j, m = b.orb.j, b.m[1]
    γj₁′, γj₁ = first.(quantum_numbers(system, a, b))
    γj₂′, γj₂ = first.(other_quantum_numbers(system, a, b))
    matrix_element_via_uncoupling((γj₁′, γj₂′, j′, m′), 𝐓ᵏq, (γj₁, γj₂, j, m))
end

"""
    matrix_element(::Tuple{S,S},
                   a::SpinOrbital{<:RelativisticOrbital},
                   X::TensorScalarProduct,
                   b::SpinOrbital{<:RelativisticOrbital}) where {S<:Union{FullSystem,TotalAngularMomentumSubSystem}}

The matrix element of a tensor scalar product, where both factors act
on the full system or the total angular momentum, evaluated in the
basis of coupled spin-orbitals, is computed using
``\\eqref{eqn:scalar-product-tensor-matrix-element}``.
"""
matrix_element(systems::Tuple{S,S},
               a::SpinOrbital{<:RelativisticOrbital},
               X::TensorScalarProduct,
               b::SpinOrbital{<:RelativisticOrbital}) where {S<:Union{FullSystem,TotalAngularMomentumSubSystem}} =
                   matrix_element(quantum_numbers(systems[1], a),
                                  X, quantum_numbers(systems[1], b))

@doc raw"""
    matrix_element(systems::Tuple{S,S},
                   a::SpinOrbital{<:RelativisticOrbital},
                   X::TensorScalarProduct,
                   b::SpinOrbital{<:RelativisticOrbital}) where {S<:SubSystem}

The matrix element of a tensor scalar product, where both factors act
on the same subsystem, evaluated in the
basis of coupled spin-orbitals, is computed using

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
"""
function matrix_element(systems::Tuple{S,S},
                        a::SpinOrbital{<:RelativisticOrbital},
                        X::TensorScalarProduct,
                        b::SpinOrbital{<:RelativisticOrbital}) where {S<:SubSystem}
    γj₁′, γj₁ = first.(quantum_numbers(systems[1], a, b))
    γj₂′, γj₂ = first.(other_quantum_numbers(systems[1], a, b))
    @δ γj₂′,γj₂ a.orb.j,b.orb.j
    matrix_element((γj₁′, a.m[1]), X, (γj₁, b.m[1]))
end

"""
    matrix_element(systems::Tuple{<:SubSystem,<:SubSystem},
                   a::SpinOrbital{<:RelativisticOrbital},
                   X::TensorScalarProduct,
                   b::SpinOrbital{<:RelativisticOrbital})

The matrix element of a tensor scalar product, where the factors act
on different subsystems, evaluated in the basis of coupled
spin-orbitals, is computed using
``\\eqref{eqn:scalar-product-tensor-matrix-element-diff-coords-coupled}``.
"""
function matrix_element(systems::Tuple{<:SubSystem,<:SubSystem},
                        a::SpinOrbital{<:RelativisticOrbital},
                        X::TensorScalarProduct,
                        b::SpinOrbital{<:RelativisticOrbital})
    γj₁′, γj₁ = first.(quantum_numbers(systems[1], a, b))
    γj₂′, γj₂ = first.(quantum_numbers(systems[2], a, b))
    matrix_element2((γj₁′, γj₂′, a.orb.j, a.m[1]), X, (γj₁, γj₂, b.orb.j, b.m[1]))
end

# ** Uncoupled orbitals

"""
    matrix_element(::Union{FullSystem,TotalAngularMomentumSubSystem},
                   a::SpinOrbital{<:Orbital},
                   𝐓ᵏq::TensorComponent,
                   b::SpinOrbital{<:Orbital})

The matrix element of a tensor acting on the full system or the total
angular momentum, evaluated in the basis of uncoupled spin-orbitals,
is computed by transforming to the coupled system via
``\\eqref{eqn:coupling}``.
"""
matrix_element(system::Union{FullSystem,TotalAngularMomentumSubSystem},
               a::SpinOrbital{<:Orbital},
               𝐓ᵏq::TensorComponent,
               b::SpinOrbital{<:Orbital}) =
    matrix_element(quantum_numbers(system, a)..., 𝐓ᵏq, quantum_numbers(system, b)...)

@doc raw"""
    matrix_element(system,
                   a::SpinOrbital{<:Orbital},
                   𝐓ᵏq::TensorComponent,
                   b::SpinOrbital{<:Orbital})

The matrix element of a tensor acting on `system`, which is a
subsystem, evaluated in the basis uncoupled spin-orbitals, is given by

```math
\begin{equation}
\begin{aligned}
&\matrixel{n_1'j_1'm_1';n_2'j_2'm_2'}{\tensor{T}^{(k)}_q(1)}{n_1j_1m_1;n_2j_2m_2} \\
=&
\delta_{n_2'n_2}\delta_{j_2'j_2}\delta_{m_2'm_2} \\
&\matrixel{n_1'j_1'm_1'}{\tensor{T}^{(k)}_q(1)}{n_1j_1m_1}.
\end{aligned}
\label{eqn:tensor-matrix-element-subsystem-uncoupled}
\tag{V13.1.39}
\end{equation}
```
"""
function matrix_element(system,
                        a::SpinOrbital{<:Orbital},
                        𝐓ᵏq::TensorComponent,
                        b::SpinOrbital{<:Orbital})
    γjm₂′,γjm₂ = other_quantum_numbers(system, a, b)
    @δ γjm₂′,γjm₂
    matrix_element(quantum_numbers(system, a), 𝐓ᵏq, quantum_numbers(system, b))
end

"""
    matrix_element(::Tuple{S,S},
                   a::SpinOrbital{<:Orbital},
                   X::TensorScalarProduct,
                   b::SpinOrbital{<:Orbital}) where {S<:Union{FullSystem,TotalAngularMomentumSubSystem}}

The matrix element of a tensor scalar product, where both factors act
on the full system or the total angular momentum, evaluated in the
basis of uncoupled spin-orbitals, is computed by transforming to the
coupled system via ``\\eqref{eqn:coupling}``, which then dispatches to
``\\eqref{eqn:scalar-product-tensor-matrix-element}``.
"""
matrix_element(systems::Tuple{S,S},
               a::SpinOrbital{<:Orbital},
               X::TensorScalarProduct,
               b::SpinOrbital{<:Orbital}) where {S<:Union{FullSystem,TotalAngularMomentumSubSystem}} =
                   matrix_element(quantum_numbers(systems[1], a)...,
                                  X, quantum_numbers(systems[1], b)...)

"""
    matrix_element(systems::Tuple{S,S},
                   a::SpinOrbital{<:Orbital},
                   X::TensorScalarProduct,
                   b::SpinOrbital{<:Orbital}) where {S<:SubSystem}

The matrix element of a tensor scalar product, where both factors act
on the same subsystem, evaluated in the basis of uncoupled
spin-orbitals, is just a special case of
``\\eqref{eqn:tensor-matrix-element-subsystem-uncoupled}`` combined
with ``\\eqref{eqn:scalar-product-tensor-matrix-element}``.
"""
function matrix_element(systems::Tuple{S,S},
                        a::SpinOrbital{<:Orbital},
                        X::TensorScalarProduct,
                        b::SpinOrbital{<:Orbital}) where {S<:SubSystem}
    γjm₂′,γjm₂ = other_quantum_numbers(systems[1], a, b)
    @δ γjm₂′,γjm₂
    matrix_element(quantum_numbers(systems[1], a),
                   X, quantum_numbers(systems[1], b))
end

"""
    matrix_element(systems::Tuple{<:SubSystem,<:SubSystem},
                   a::SpinOrbital{<:Orbital},
                   X::TensorScalarProduct,
                   b::SpinOrbital{<:Orbital})

The matrix element of a tensor scalar product, where the factors act
on different subsystems, evaluated in the basis of uncoupled
spin-orbitals, is computed using
``\\eqref{eqn:scalar-product-tensor-matrix-element-diff-coords-uncoupled}``.
"""
function matrix_element(systems::Tuple{<:SubSystem,<:SubSystem},
                        a::SpinOrbital{<:Orbital},
                        X::TensorScalarProduct,
                        b::SpinOrbital{<:Orbital})
    γjm₁′, γjm₁ = quantum_numbers(systems[1], a, b)
    γjm₂′, γjm₂ = quantum_numbers(systems[2], a, b)
    matrix_element2(γjm₁′, γjm₂′, X, γjm₁, γjm₂)
end

# * Old stuff

complementary_space_factor(::Union{FullSystem,TotalAngularMomentumSubSystem}, _, _, _) = 1

@doc raw"""
    complementary_space_factor(::SpatialSubSystems, k,
                               o′::SpinOrbital{<:Orbital},
                               o::SpinOrbital{<:Orbital})

Tensors acting on the spatial coordinates only, are diagonal in spin
space:

```math
\begin{equation}
\tag{V13.2.3}
\matrixel{n'\ell'm_\ell';s'm_s'}{\tensor{M}^{(k)}_q(r,\theta,\phi)}{n\ell m_\ell;sm_s} =
\delta_{ss'}\delta_{m_sm_s'}
\matrixel{n'\ell'm_\ell'}{\tensor{M}^{(k)}_q(r,\theta,\phi)}{n\ell m_\ell}
\end{equation}
```
"""
function complementary_space_factor(::SpatialSubSystems, k,
                                    o′::SpinOrbital{<:Orbital},
                                    o::SpinOrbital{<:Orbital})
    (s′,mₛ′),(s,mₛ) = quantum_numbers(SpinSubSystem(), o′, o)
    @δ s,s′ mₛ,mₛ′
end

@doc raw"""
    complementary_space_factor(::SpinSubSystem, k,
                               o′::SpinOrbital{<:Orbital},
                               o::SpinOrbital{<:Orbital})

Tensors acting on the spin coordinates only, are diagonal in the spatial coordinates:

```math
\begin{equation}
\tag{V13.2.4}
\matrixel{n'\ell'm_\ell';s'm_s'}{\tensor{N}^{(k)}_q(\xi)}{n\ell m_\ell;sm_s} =
\delta_{\ell\ell'}\delta_{m_\ell m_\ell'}
\matrixel{s'm_s'}{\tensor{N}^{(k)}_q(\xi)}{sm_s}
\end{equation}
```
"""
function complementary_space_factor(::SpinSubSystem, k,
                                    o′::SpinOrbital{<:Orbital},
                                    o::SpinOrbital{<:Orbital})
    (ℓ′,mℓ′),(ℓ,mℓ) = quantum_numbers(OrbitalAngularMomentumSubSystem(), o′, o)
    @δ ℓ,ℓ′ mℓ,mℓ′
end

@doc raw"""
    complementary_space_factor(::SpatialSubSystems, k,
                               o′::SpinOrbital{<:RelativisticOrbital},
                               o::SpinOrbital{<:RelativisticOrbital})

Tensors acting on the spatial coordinates only, are diagonal in spin
space; in the coupled basis, the following uncoupling formula must be
employed:

```math
\begin{equation}
\tag{V13.2.5}
\redmatrixel{n'\ell's'J'}{\tensor{M}^{(k)}}{n\ell s J} =
\delta_{ss'}
(-)^{J+\ell'+s'+k}
\angroot{JJ'}
\wignersixj{\ell&s&J\\J'&k&\ell'}
\redmatrixel{n'\ell'}{\tensor{M}^{(k)}}{n\ell}
\end{equation}
```
"""
function complementary_space_factor(::SpatialSubSystems, k,
                                    o′::SpinOrbital{<:RelativisticOrbital},
                                    o::SpinOrbital{<:RelativisticOrbital})
    ((ℓ′,s′,J′),mₛ′),((ℓ,s,J),mₛ) = quantum_numbers(TotalAngularMomentumSubSystem(), o′, o)
    @δ s,s′
    powneg1(Int(J+ℓ′+s′+k))*∏(J,J′)*wigner6j(ℓ,  s, J,
                                             J′, k, ℓ′)
end

@doc raw"""
    complementary_space_factor(::SpinSubSystem, k,
                               o′::SpinOrbital{<:RelativisticOrbital},
                               o::SpinOrbital{<:RelativisticOrbital})

Tensors acting on the spin coordinates only, are diagonal in the
spatial coordinates; in the coupled basis, the following uncoupling
formula must be employed:

```math
\begin{equation}
\tag{V13.2.6}
\redmatrixel{n'\ell's'J'}{\tensor{N}^{(k)}}{n\ell s J} =
\delta_{\ell\ell'}
(-)^{J'+\ell+s+k}
\angroot{JJ'}
\wignersixj{s&\ell&J\\J'&k&s'}
\redmatrixel{n's'}{\tensor{N}^{(k)}}{ns}
\end{equation}
```
"""
function complementary_space_factor(::SpinSubSystem, k,
                                    o′::SpinOrbital{<:RelativisticOrbital},
                                    o::SpinOrbital{<:RelativisticOrbital})
    ((ℓ′,s′,J′),mₛ′),((ℓ,s,J),mₛ) = quantum_numbers(TotalAngularMomentumSubSystem(), o′, o)
    @δ ℓ,ℓ′
    powneg1(Int(J′+ℓ+s+k))*∏(J,J′)*wigner6j(s,  ℓ, J,
                                            J′, k, s′)
end

total_system(s, _) = s
total_system(_, ::SpinOrbital{<:RelativisticOrbital}) = TotalAngularMomentumSubSystem()

"""
    rme_j′j(o′::SpinOrbital, T::Tensor, o::SpinOrbital)

Return the reduced matrix element `⟨o′||T||o⟩` and the two angular
momenta pertaining to the tensor `T`, along with their projections.
"""
function rme_j′j(o′::SpinOrbital, T::Tensor, o::SpinOrbital)
    o′, T, o
    s = system(T)
    f = complementary_space_factor(s, rank(T), o′, o)
    (γ′,_),(γ,_) = quantum_numbers(s, o′, o)
    # The j′,m′,j,m appearing in the prefactor of the Wigner–Eckart
    # theorem are the total angular momenta and their projections for
    # coupled spin-orbitals and the angular momenta and their
    # projections of the subsystem acted upon by the tensor T for
    # uncoupled spin-orbitals.
    (γ̃′,m′),(γ̃,m) = quantum_numbers(total_system(s, o′), o′, o)
    j′,j = last(γ̃′),last(γ̃)
    (iszero(f) ? 0 : f*rme(γ′, T, γ)),(j′,m′),(j,m)
end

couples(o′::SpinOrbital, T::Tensor, o::SpinOrbital) =
    !iszero(first(rme_j′j(o′, T, o)))

"""
    dot((a,b), X, (c,d))

Perform the spin-angular integration of the scalar-product tensor
`X≡(T⁽ᵏ⁾⋅U⁽ᵏ⁾)`, where `T` acts on the coordinate of orbitals `a` &
`c` and similarly, `U` acts on the coordinate of orbitals `b` &
`d`, according to Eq. (13.1.26) of Varshalovich (1988).

# Examples

```jldoctest
julia> a = SpinOrbital(ro"1s", half(1))
1s(1/2)

julia> b = SpinOrbital(ro"3d", half(1))
3d(1/2)

julia> 𝐂 = SphericalTensor(2)
𝐂̂⁽²⁾

julia> X = dot(𝐂,𝐂)
(𝐂̂⁽²⁾⋅𝐂̂⁽²⁾)

julia> dot((a,b), X, (b,a))
0.12000000000000001
```
"""
function LinearAlgebra.dot((a,b)::Tuple, X::TensorScalarProduct, (c,d)::Tuple)
    T,U = X.T,X.U
    k = rank(T)

    Tr,(ja,ma),(jc,mc) = rme_j′j(a,T,c)
    iszero(Tr) && return 0
    Ur,(jb,mb),(jd,md) = rme_j′j(b,U,d)
    iszero(Ur) && return 0

    α = Int(ma-mc)
    (α != md-mb || abs(α) > k) && return 0

    inv(∏(ja,jb)) * powneg1(-α) *
        clebschgordan(jc, mc, k, α, ja, ma) *
        clebschgordan(jd, md, k, -α, jb, mb) *
        Tr*Ur
end

export matrix_element, matrix_element2, matrix_element_via_uncoupling
