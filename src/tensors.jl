"""
    Tensor{k,label}

Abstract base for any tensor of rank `k`.
"""
abstract type Tensor{k,label} end

(::Type{T})(k::Int) where {T<:Tensor} = T{k}()

LinearAlgebra.rank(::Tensor{k}) where k = k
couples(a, ::Type{<:Tensor}, b) = true

# This is only true for SO(3)
components(::Tensor{k}) where k = -k:k

function Base.show(io::IO, ::Tensor{k,label}) where {k,label}
    write(io,to_boldface(label))
    write(io,"⁽",to_superscript(k),"⁾")
end

jmⱼ(o′::SpinOrbital{<:RelativisticOrbital}, Tᵏ::Tensor, o::SpinOrbital{<:RelativisticOrbital}) =
    o′.orb.j,o′.m[1],o.orb.j,o.m[1],true

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

jmⱼ(o′::SpinOrbital, Tᵏq::TensorComponent, o::SpinOrbital) =
    jmⱼ(o′, Tᵏq.tensor, o)

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
(𝐂⁽⁴⁾⋅𝐂⁽⁴⁾)
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

    ja,ma,jc,mc,Tdiag = jmⱼ(a,T,c)
    jb,mb,jd,md,Udiag = jmⱼ(b,U,d)

    Tdiag && Udiag || return zero(Float64)

    α = Int(ma-mc)
    (α != md-mb || abs(α) > k) && return zero(Float64)

    inv(∏(ja,jb)) * powneg1(-α) *
        clebschgordan(jc, mc, k, α, ja, ma) *
        clebschgordan(jd, md, k, -α, jb, mb) *
        rme(a.orb, X.T, c.orb) *
        rme(b.orb, X.U, d.orb)
end

# * Linear combination of tensors

"""
    LinearCombinationTensor

Represents a linear combination of tensor components.
"""
const LinearCombinationTensor{T<:Tensor,N<:Number} = LinearCombination{<:TensorComponent{<:T},N}

@linearly_combinable TensorComponent

# * Wigner--Eckart

"""
    wigner_eckart(j′, m′, T⁽ᵏ⁾q, j, m)

Computes the (spin-angular part of the) matrix element
`⟨n′j′m′|Tᵏq|njm⟩`, where `T⁽ᵏ⁾q` is the `q`th component of a tensor
of rank `k`, using the definition of Eq. (13.1.2) in Varshalovich (1988).
"""
function wigner_eckart(j′, m′, Tᵏq::TensorComponent, j, m)
    Tᵏ,q = Tᵏq.tensor, Tᵏq.q
    k = rank(Tᵏ)
    powneg1(Int(j′-m′))*wigner3j(j′, k, j,
                                 -m′, q, m)*rme(j′, Tᵏ, j)
end

wigner_eckart(j′, m′, lct::LinearCombinationTensor, j, m) =
    sum(c*wigner_eckart(j′, m′, Tᵏq, j, m) for (Tᵏq,c) in lct)

"""
    wigner_eckart(o′, T⁽ᵏ⁾q, o)

Computes the (spin-angular part of the) matrix element `⟨o′|Tᵏq|o⟩`,
where `T⁽ᵏ⁾q` is the `q`th component of a tensor of rank `k`, using
the definition of Eq. (13.1.2) in Varshalovich (1988).
"""
function wigner_eckart(o′::SpinOrbital, Tᵏq::TensorComponent, o::SpinOrbital)
    Tᵏ = Tᵏq.tensor
    r = rme(o′.orb, Tᵏ, o.orb)
    iszero(r) && return r
    k,q = rank(Tᵏ), Tᵏq.q
    j′,m′,j,m,isdiagonal = jmⱼ(o′,Tᵏq,o)
    isdiagonal || return zero(r)
    powneg1(Int(j′-m′))*wigner3j(j′, k, j,
                                 -m′, q, m)*r
end

wigner_eckart(o′, lct::LinearCombinationTensor, o) =
    sum(c*wigner_eckart(o′, Tᵏq, o) for (Tᵏq,c) in lct)

"""
    dot(o′, T, o)

Calculates the matrix element `⟨o′|T|o⟩` using [`wigner_eckart`](@ref).
"""
LinearAlgebra.dot(o′::SpinOrbital, T::Union{TensorComponent,LinearCombinationTensor}, o::SpinOrbital) =
    wigner_eckart(o′, T, o)

export Tensor, TensorComponent,
    TensorProduct, TensorScalarProduct,
    wigner_eckart
