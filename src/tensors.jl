"""
    Tensor

Abstract base for any tensor.
"""
abstract type Tensor end

couples_spin(a, ::Type{Tensor}, b) = true

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
    rme(ℓ′,Cᵏ,ℓ)

Calculate the reduced matrix element `⟨ℓ′||C⁽ᵏ⁾||ℓ⟩` of the spherical
tensor of rank `k`. Condon–Shortley phase convention and using the
definition of Eq. (13.2.107) in Varshalovich (1988).
"""
rme(ℓ′,Cᵏ::SphericalTensor,ℓ) = powneg1(ℓ-Cᵏ.k)*∏(ℓ,ℓ′)*wigner3j(ℓ,Cᵏ.k,ℓ′,
                                                                 0,0,0)

couples_spin(a, ::Type{SphericalTensor}, b) = a == b

export Tensor, SphericalTensor, rme
