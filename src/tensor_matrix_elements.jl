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
function matrix_element((γj′, m′)::Tuple{<:Any, <:Number},
                        Tᵏq::TensorComponent,
                        (γj, m)::Tuple{<:Any, <:Number})
    Tᵏ = parent(Tᵏq)
    r = rme(γj′, Tᵏ, γj)
    iszero(r) && return 0
    j′,j = last(γj′), last(γj)
    c = powneg1(Int(j′-m′))*wigner3j(j′, rank(Tᵏ), j,
                                     -m′, component(Tᵏq), m)
    iszero(c) && return 0
    c*r
end

matrix_element(γj′m′::Tuple{<:Any, <:Number},
               lct::LinearCombinationTensor,
               γjm::Tuple{<:Any, <:Number}) =
    sum(filter!(!iszero, [c*matrix_element(γj′m′, Tᵏq, γjm) for (Tᵏq,c) in lct]))

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
    matrix_element((γj₁′, γj₂′, j′, m′), 𝐓ᵏq, (γj₁, γj₂, j, m))

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

julia> matrix_element((1, half(1), half(3), half(3)),
                      𝐋₀, (1, half(1), half(3), half(3)))
0.9999999999999999

julia> matrix_element((1, 1), 𝐋₀, (1, 1)) # For comparison
0.9999999999999999

julia> 𝐒₀ = TensorComponent(SpinAngularMomentum(), 0)
𝐒̂⁽¹⁾₀

julia> matrix_element((half(1), 1, half(3), half(3)),
                      𝐒₀, (half(1), 1, half(3), half(3)))
0.49999999999999994

julia> matrix_element((half(1),half(1)), 𝐒₀, (half(1),half(1)))
0.49999999999999994
```

"""
function matrix_element((γj₁′, γj₂′, j′, m′)::Tuple{<:Any, <:Any, <:Number, <:Number},
                        𝐓ᵏq::TensorComponent,
                        (γj₁, γj₂, j, m)::Tuple{<:Any, <:Any, <:Number, <:Number})
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

function matrix_element((γj₁′, γj₂′, j′, m′)::Tuple{<:Any, Tuple{}, <:Number, <:Number},
                        𝐓ᵏq::TensorComponent,
                        (γj₁, γj₂, j, m)::Tuple{<:Any, Tuple{}, <:Number, <:Number})
    last(γj₁′) == j′ && last(γj₁) == j ||
        throw(ArgumentError("Angular momentum mismatch"))
    matrix_element((γj₁′, m′), 𝐓ᵏq, (γj₁, m))
end


# ** Uncoupled basis states

@doc raw"""
    matrix_element2(γjm₁′, γjm₂′, X::TensorScalarProduct, γjm₁, γjm₂)

The matrix element of a scalar product of two tensors acting on
different coordinates is given by (in the uncoupled basis)

```math
\begin{equation}
\begin{aligned}
&\matrixel{\gamma_1'j_1'm_1';\gamma_2'j_2'm_2'}{[\tensor{P}^{(k)}(1)\cdot\tensor{Q}^{(k)}(2)]}{\gamma_1j_1m_1;\gamma_2j_2m_2}\\
=&
\frac{1}{\angroot{j_1'j_2'}}
\sum_\alpha(-)^{-\alpha}
C_{j_1m_1;k,\alpha}^{j_1'm_1'}
C_{j_2m_2;k,-\alpha}^{j_2'm_2'}\\
&\times
\redmatrixel{\gamma_1'j_1'}{\tensor{P}^{(k)}(1)}{\gamma_1j_1}
\redmatrixel{\gamma_2'j_2'}{\tensor{Q}^{(k)}(2)}{\gamma_2j_2} \\
\equiv&
\sum_\alpha
(-)^{-\alpha}
\matrixel{\gamma_1'j_1'm_1'}{\tensor{P}^{(k)}_{\alpha}(1)}{\gamma_1j_1m_1}
\matrixel{\gamma_2'j_2'm_2'}{\tensor{Q}^{(k)}_{-\alpha}(2)}{\gamma_2j_2m_2}
\end{aligned}
\tag{V13.1.26}
\label{eqn:scalar-product-tensor-matrix-element-diff-coords-uncoupled}
\end{equation}
```

Since the [Clebsch–Gordan coefficients](@ref) can be rewritten using 3j
symbols and the 3j symbols vanish unless $m_1 + \alpha - m_1' = m_2 -
\alpha - m_2' = 0$, we have

```math
\alpha = m_1' - m_1 = m_2-m_2'
```

This case occurs in two-body interactions, such as the [Coulomb
interaction](@ref), where ``1',1`` and ``2',2`` are pairs of orbitals
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
julia> X = OrbitalAngularMomentum()⋅SpinAngularMomentum()
(𝐋̂⁽¹⁾⋅𝐒̂⁽¹⁾)

julia> matrix_element2((1, 1), (half(1), half(1)),
                      X, (1,1), (half(1), half(1)))
0.4999999999999999

julia> 1/2*(half(3)*(half(3)+1)-1*(1+1)-half(1)*(half(1)+1)) # 1/2(J(J+1)-L(L+1)-S(S+1))
0.5
```
"""
function matrix_element2(γjm₁′, γjm₂′, X::TensorScalarProduct, γjm₁, γjm₂)
    T,U = X.T,X.U
    k = rank(T)

    α = Int(last(γjm₁′)-last(γjm₁))
    (α != last(γjm₂)-last(γjm₂′) || abs(α) > k) && return 0

    powneg1(-α)*
    matrix_element(γjm₁′, TensorComponent(T,α), γjm₁)*
    matrix_element(γjm₂′, TensorComponent(U,-α), γjm₂)
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

```jldoctest
julia> X = OrbitalAngularMomentum()⋅SpinAngularMomentum()
(𝐋̂⁽¹⁾⋅𝐒̂⁽¹⁾)

julia> matrix_element2((1, half(1), half(3), half(3)),
                       X, (1, half(1), half(3), half(3)))
0.4999999999999999
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

julia> # Same as previous, but with spin down
       dot(a, TensorComponent(𝐋, 0), SpinOrbital(o"2p", 1, -half(1)))
0

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

"""
    dot((a,b)::Tuple, X::TensorScalarProduct, (c,d)::Tuple)

Compute the matrix element `⟨a(1)b(2)|𝐓(1)⋅𝐔(2)|c(1)d(2)⟩` in the
basis of spin-orbitals, dispatching the correct low-level function
`matrix_element` depending on the value of `system(X)`, where `X` is
the scalar product tensor.

# Examples

```jldoctest
julia> 𝐊⁰,𝐊² = CoulombTensor(0),CoulombTensor(2)
(𝐊̂⁽⁰⁾, 𝐊̂⁽²⁾)

julia> a,b = SpinOrbital(o"1s", 0, half(1)),SpinOrbital(o"3d", 0, half(1))
(1s₀α, 3d₀α)

julia> dot((a,b), 𝐊⁰⋅𝐊⁰, (a,b))
1.0

julia> dot((a,b), 𝐊²⋅𝐊², (b,a))
0.19999999999999998

julia> a,b = SpinOrbital(ro"1s", half(1)),SpinOrbital(ro"3d", half(1))
(1s(1/2), 3d(1/2))

julia> dot((a,b), 𝐊⁰⋅𝐊⁰, (a,b))
1.0

julia> dot((a,b), 𝐊²⋅𝐊², (b,a))
0.12
```
"""
LinearAlgebra.dot((a,b)::Tuple, X::TensorScalarProduct, (c,d)::Tuple) =
    matrix_element(system(X), (a,b), X, (c,d))

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
    matrix_element((γj₁′, γj₂′, j′, m′), 𝐓ᵏq, (γj₁, γj₂, j, m))
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
function matrix_element((s,_)::Tuple{S,S},
                        a::SpinOrbital{<:RelativisticOrbital},
                        X::TensorScalarProduct,
                        b::SpinOrbital{<:RelativisticOrbital}) where {S<:SubSystem}
    γj₁′, γj₁ = first.(quantum_numbers(s, a, b))
    γj₂′, γj₂ = first.(other_quantum_numbers(s, a, b))
    @δ γj₂′,γj₂ a.orb.j,b.orb.j
    matrix_element((γj₁′, a.m[1]), X, (γj₁, b.m[1]))
end

"""
    matrix_element((s₁,s₂)::Tuple{<:SubSystem,<:SubSystem},
                   a::SpinOrbital{<:RelativisticOrbital},
                   X::TensorScalarProduct,
                   b::SpinOrbital{<:RelativisticOrbital})

The matrix element of a tensor scalar product, where the factors act
on different subsystems, evaluated in the basis of coupled
spin-orbitals, is computed using
``\\eqref{eqn:scalar-product-tensor-matrix-element-diff-coords-coupled}``.
"""
function matrix_element((s₁,s₂)::Tuple{<:SubSystem,<:SubSystem},
                        a::SpinOrbital{<:RelativisticOrbital},
                        X::TensorScalarProduct,
                        b::SpinOrbital{<:RelativisticOrbital})
    γj₁′, γj₁ = first.(quantum_numbers(s₁, a, b))
    γj₂′, γj₂ = first.(quantum_numbers(s₂, a, b))
    matrix_element2((γj₁′, γj₂′, a.orb.j, a.m[1]), X, (γj₁, γj₂, b.orb.j, b.m[1]))
end

"""
    matrix_element((s₁,s₂)::Tuple{<:System,<:System},
                   (a,b)::Tuple{<:SpinOrbital{<:RelativisticOrbital},
                                <:SpinOrbital{<:RelativisticOrbital}},
                   X::TensorScalarProduct,
                   (c,d)::Tuple{<:SpinOrbital{<:RelativisticOrbital},
                                <:SpinOrbital{<:RelativisticOrbital}})

The matrix element of a tensor scalar product, where the factors act
on the orbital pairs `a`,`c` and `b`,`d`, respectively, evaluated in
the basis of coupled spin-orbitals, is computed using
``\\eqref{eqn:scalar-product-tensor-matrix-element-diff-coords-uncoupled}``
together with ``\\eqref{eqn:uncoupling}`` applied to each matrix
element between single orbitals; the individual spin-orbitals couple
``\\ell`` and ``s`` to form a total ``j``, but they are not further
coupled to e.g. a term.
"""
function matrix_element((s₁,s₂)::Tuple{<:System,<:System},
                        (a,b)::Tuple{<:SpinOrbital{<:RelativisticOrbital},
                                     <:SpinOrbital{<:RelativisticOrbital}},
                        X::TensorScalarProduct,
                        (c,d)::Tuple{<:SpinOrbital{<:RelativisticOrbital},
                                     <:SpinOrbital{<:RelativisticOrbital}})
    γjm₁′ = (first(quantum_numbers(s₁, a)), first(other_quantum_numbers(s₁, a)),
             a.orb.j, a.m[1])
    γjm₂′ = (first(quantum_numbers(s₂, b)), first(other_quantum_numbers(s₂, b)),
             b.orb.j, b.m[1])
    γjm₁ = (first(quantum_numbers(s₁, c)), first(other_quantum_numbers(s₁, c)),
            c.orb.j, c.m[1])
    γjm₂ = (first(quantum_numbers(s₂, d)), first(other_quantum_numbers(s₂, d)),
            d.orb.j, d.m[1])

    matrix_element2(γjm₁′, γjm₂′, X, γjm₁, γjm₂)
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
function matrix_element((s,_)::Tuple{S,S},
                        a::SpinOrbital{<:Orbital},
                        X::TensorScalarProduct,
                        b::SpinOrbital{<:Orbital}) where {S<:SubSystem}
    @δ other_quantum_numbers(s, a, b)
    matrix_element(quantum_numbers(s, a),
                   X, quantum_numbers(s, b))
end

"""
    matrix_element((s₁,s₂)::Tuple{<:SubSystem,<:SubSystem},
                   a::SpinOrbital{<:Orbital},
                   X::TensorScalarProduct,
                   b::SpinOrbital{<:Orbital})

The matrix element of a tensor scalar product, where the factors act
on different subsystems, evaluated in the basis of uncoupled
spin-orbitals, is computed using
``\\eqref{eqn:scalar-product-tensor-matrix-element-diff-coords-uncoupled}``.
"""
function matrix_element((s₁,s₂)::Tuple{<:SubSystem,<:SubSystem},
                        a::SpinOrbital{<:Orbital},
                        X::TensorScalarProduct,
                        b::SpinOrbital{<:Orbital})
    γjm₁′, γjm₁ = quantum_numbers(s₁, a, b)
    γjm₂′, γjm₂ = quantum_numbers(s₂, a, b)
    matrix_element2(γjm₁′, γjm₂′, X, γjm₁, γjm₂)
end

"""
    matrix_element((s₁,s₂)::Tuple{<:System,<:System},
                   (a,b)::Tuple{<:SpinOrbital{<:Orbital},
                                <:SpinOrbital{<:Orbital}},
                   X::TensorScalarProduct,
                   (c,d)::Tuple{<:SpinOrbital{<:Orbital},
                                <:SpinOrbital{<:Orbital}})

The matrix element of a tensor scalar product, where the factors act
on the orbital pairs `a`,`c` and `b`,`d`, respectively, evaluated in
the basis of uncoupled spin-orbitals, is computed using
``\\eqref{eqn:scalar-product-tensor-matrix-element-diff-coords-uncoupled}``.
"""
function matrix_element((s₁,s₂)::Tuple{<:System,<:System},
                        (a,b)::Tuple{<:SpinOrbital{<:Orbital},
                                     <:SpinOrbital{<:Orbital}},
                        X::TensorScalarProduct,
                        (c,d)::Tuple{<:SpinOrbital{<:Orbital},
                                     <:SpinOrbital{<:Orbital}})
    @δ other_quantum_numbers(s₁, a, c) other_quantum_numbers(s₂, b, d)
    γjm₁′, γjm₁ = quantum_numbers(s₁, a, c)
    γjm₂′, γjm₂ = quantum_numbers(s₂, b, d)
    matrix_element2(γjm₁′, γjm₂′, X, γjm₁, γjm₂)
end

export matrix_element, matrix_element2
