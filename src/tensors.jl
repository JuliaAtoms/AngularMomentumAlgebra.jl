"""
    Tensor

Abstract base for any tensor.
"""
abstract type Tensor end

couples_spin(a, ::Type{Tensor}, b) = true

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
        abs(q) ≤ k ||
            throw(ArgumentError("Tensor component $(q) not possible for tensor $(tensor) of rank $k"))
        new{T}(tensor, q)
    end
end

# Type{T}(k::Int, q::Int) where {T<:Tensor} =
#     TensorComponent(T(k), q)

function Base.show(io::IO, Tq::TensorComponent)
    show(io, Tq.tensor)
    write(io, to_subscript(Tq.q))
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
of rank `k`, using the definition of Eq. (13.1.1) in Varshalovich (1988).
"""
function wigner_eckart(j′, m′, Tkq::TensorComponent, j, m)
    Tk,q = Tkq.tensor, Tkq.q
    k = rank(Tk)
    powneg1(j′-m′)*wigner3j(j′, k, j,
                            -m′, q, m)*rme(j′, Tk, j)
end

wigner_eckart(j′, m′, lct::LinearCombinationTensor, j, m) =
    sum(c*wigner_eckart(j′, m′, Tkq, j, m) for (Tkq,c) in lct)

# * Spherical tensors

"""
    SphericalTensor(k)

Construct a spherical tensor of rank `k`.
"""
struct SphericalTensor <: Tensor
    k::Int
end

function Base.show(io::IO, C::SphericalTensor)
    write(io,to_boldface("C"))
    write(io,"⁽",to_superscript(C.k),"⁾")
end

"""
    rank(C::SphericalTensor)

Returns the rank of the tensor `C`.
"""
LinearAlgebra.rank(C::SphericalTensor) = C.k

"""
    rme(ℓ′,Cᵏ,ℓ)

Calculate the reduced matrix element `⟨ℓ′||C⁽ᵏ⁾||ℓ⟩` of the spherical
tensor of rank `k`. Condon–Shortley phase convention and using the
definition of Eq. (13.2.107) in Varshalovich (1988).
"""
rme(ℓ′,Cᵏ::SphericalTensor,ℓ) = powneg1(ℓ-Cᵏ.k)*∏(ℓ,ℓ′)*wigner3j(ℓ,Cᵏ.k,ℓ′,
                                                                 0,0,0)

couples_spin(a, ::Type{SphericalTensor}, b) = a == b

export Tensor, TensorComponent, SphericalTensor, rme, wigner_eckart
