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
\tag{V13.1.11}
\end{equation}
```

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
    # 𝐓ᵏ = parent(𝐓ᵏq)
    # k = rank(𝐓ᵏ)
    # q = component(𝐓ᵏq)
    for j′ ∈ j′s
        # c′ = clebschgordan(j₁′, m₁′, j₂′, m₂′, j′, m′)/∏(j′)
        c′ = clebschgordan(j₁′, m₁′, j₂′, m₂′, j′, m′)
        for j ∈ js
            # c = c′*clebschgordan(j₁, m₁, j₂, m₂, j, m)*clebschgordan(j, m, k, q, j′, m′)
            # r = rme((γj′...,j′), 𝐓ᵏ, (γj...,j))
            # iszero(r) && continue
            # v += c*r
            me = matrix_element(((γj′...,j′), m′), 𝐓ᵏq, ((γj...,j), m))
            iszero(me) && continue
            c = c′*clebschgordan(j₁, m₁, j₂, m₂, j, m)
            v += c*me
        end
    end

    # powneg1(2rank(𝐓ᵏ))*v
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

# """
#     dot(o′, T, o)

# Calculates the matrix element `⟨o′|T|o⟩` using [`matrix_element`](@ref).
# """
# LinearAlgebra.dot(o′::SpinOrbital, T::TensorComponent, o::SpinOrbital) =
#     matrix_element(o′, T, o)

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

# function LinearAlgebra.dot(a::SpinOrbital{<:Orbital}, X::TensorScalarProduct, b::SpinOrbital{<:Orbital})
#     T,U = X.T,X.U
#     # @show a, T, U, b
#     sT = system(T)
#     sU = system(U)
#     # @show sT, sU
#     sT != sU && return dot((a,a), X, (b,b))
#     0
# end

function LinearAlgebra.dot(a::SpinOrbital, T::Union{TensorComponent,TensorScalarProduct}, b::SpinOrbital)
    @show a, T, b
    0
end

LinearAlgebra.dot(o′, lct::LinearCombinationTensor, o) =
    sum(filter!(!iszero, [c*dot(o′, Tᵏq, o) for (Tᵏq,c) in lct]))

export matrix_element, matrix_element_via_uncoupling
