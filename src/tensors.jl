"""
    Tensor{k,label}

Abstract base for any tensor of rank `k`.
"""
abstract type Tensor{k,label} end

(::Type{T})(k::Int) where {T<:Tensor} = T{k}()

LinearAlgebra.rank(::Tensor{k}) where k = k

couples(a, ::Type{<:Tensor}, b) = true

"""
    system(::Tensor)

A general tensor acts on the full system, i.e. all coordinates.
"""
system(::Tensor) = FullSystem()

# This is only true for SO(3)
components(::Tensor{k}) where k = -k:k

function Base.show(io::IO, ::Tensor{k,label}) where {k,label}
    write(io,to_boldface(label)*"̂")
    write(io,"⁽",to_superscript(k),"⁾")
end

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
\delta_{ss'}\delta{m_sm_s'}
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
\delta_{\ell\ell'}\delta{m_\ell m_\ell'}
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
    # coupled spin-orbtials and the angular momenta and their
    # projections of the subsystem acted upon by the tensor T for
    # uncoupled spin-orbitals.
    (γ̃′,m′),(γ̃,m) = quantum_numbers(total_system(s, o′), o′, o)
    j′,j = last(γ̃′),last(γ̃)
    (iszero(f) ? 0 : f*rme(γ′, T, γ)),(j′,m′),(j,m)
end

couples(o′::SpinOrbital, T::Tensor, o::SpinOrbital) =
    !iszero(first(rme_j′j(o′, T, o)))

# * Tensor components

"""
    TensorComponent(tensor, q)

Represents the `q`th component of a `tensor`; `abs(q) ≤ rank(tensor)`.
"""
struct TensorComponent{T<:Tensor}
    tensor::T
    q::Int
    function TensorComponent(tensor::T, q::Int) where {T<:Tensor}
        k = rank(tensor)
        q ∈ components(tensor) ||
            throw(ArgumentError("Tensor component $(q) not possible for tensor $(tensor) of rank $k"))
        new{T}(tensor, q)
    end
end

function Base.show(io::IO, Tq::TensorComponent)
    show(io, Tq.tensor)
    write(io, to_subscript(Tq.q))
end

Base.parent(Tᵏq::TensorComponent) = Tᵏq.tensor
component(Tᵏq::TensorComponent) = Tᵏq.q

# * Tensor products

@doc raw"""
    TensorProduct{K}(T, U)

A tensor of rank `K` formed from the product of the tensors `T` and
`U`, according to

```math
\begin{equation}
\tag{V3.1.20}
\tensor{X}^{(K)}_Q \equiv
\{\tensor{T}^{(k_1)}\tensor{U}^{(k_2)}\}^{(K)}_Q \defd
\tensor{T}^{(k_1)}_{q_1}
\tensor{U}^{(k_2)}_{q_2}
C_{k_1q_1k_2q_2}^{KQ}
\end{equation}
```
"""
struct TensorProduct{K,A<:Tensor,B<:Tensor} <: Tensor{K,'X'}
    T::A
    U::B
end
TensorProduct(K::Int, T::A, U::B) where {A<:Tensor, B<:Tensor} =
    TensorProduct{K,A,B}(T, U)

function Base.show(io::IO, X::TensorProduct{K}) where K
    write(io,"{")
    show(io, X.T)
    show(io, X.U)
    write(io,"}⁽",to_superscript(K),"⁾")
end

@doc raw"""
    TensorScalarProduct(T, U)

A tensor of rank 0 formed from the product of the tensors `T` and
`U` (which have to have the same rank), according to

```math
\begin{equation}
\tag{V3.1.30,35}
(\tensor{T}^{(k)} \cdot
\tensor{U}^{(k)}) \defd
(-)^k\angroot{k}
\{\tensor{T}^{(k)}
\tensor{U}^{(k)}\}^{(0)}_0 \equiv
(-)^q
\tensor{T}^{(k)}_{q}
\tensor{U}^{(k)}_{-q}
\end{equation}
```
"""
struct TensorScalarProduct{A<:Tensor,B<:Tensor} <: Tensor{0,'X'}
    T::A
    U::B
    function TensorScalarProduct(T::A, U::B) where {A<:Tensor,B<:Tensor}
        rank(T) == rank(U) ||
            throw(ArgumentError("Cannot form scalar product of tensors of differing ranks"))
        new{A,B}(T,U)
    end
end

function Base.show(io::IO, X::TensorScalarProduct)
    write(io,"(")
    show(io, X.T)
    write(io,"⋅")
    show(io, X.U)
    write(io,")")
end

"""
    dot(T::Tensor, U::Tensor)

Form the scalar product of the two tensors `T` and `U`, which need to
have the same rank.

# Examples

```jldoctest
julia> SphericalTensor(4)⋅SphericalTensor(4)
(𝐂̂⁽⁴⁾⋅𝐂̂⁽⁴⁾)
```
"""
LinearAlgebra.dot(T::Tensor, U::Tensor) =
    TensorScalarProduct(T, U)

"""
    integrate_spinors((a,b), X, (c,d))

Perform the spin-angular integration of the scalar-product tensor
`X≡(T⁽ᵏ⁾⋅U⁽ᵏ⁾)`, where `T` acts on the coordinate of orbitals `a` &
`c` and similarly, `U` acts on the coordinate of orbitals `b` &
`d`, according to Eq. (13.1.26) of Varshalovich (1988).
"""
function integrate_spinors((a,b), X::TensorScalarProduct, (c,d))
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

# * Linear combination of tensors

"""
    LinearCombinationTensor

Represents a linear combination of tensor components.
"""
const LinearCombinationTensor{T<:Tensor,N<:Number} = LinearCombination{<:TensorComponent{<:T},N}

@linearly_combinable TensorComponent

export Tensor, TensorComponent,
    TensorProduct, TensorScalarProduct,
    system
