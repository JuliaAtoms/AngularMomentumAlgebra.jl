wigner_eckart(::FullSystem, ::SpinOrbital{<:Orbital}, ::TensorComponent, ::SpinOrbital{<:Orbital}) =
    throw(ArgumentError("Wigner–Eckart not applicable to the full system for uncoupled orbitals"))

"""
    wigner_eckart(o′::SpinOrbital, Tᵏq::TensorComponent, o::SpinOrbital)

Calculate the matrix element `⟨o′|Tᵏq|o⟩` via Wigner–Eckart's theorem,
automatically applying the appropriate uncoupling formulas (V13.2.5,6)
of Varshalovich (1988) if `o′` and `o` are coupled spin-orbitals and
`Tᵏ` acts on one part of the system only.
"""
function wigner_eckart(o′::SpinOrbital, Tᵏq::TensorComponent, o::SpinOrbital)
    Tᵏ = parent(Tᵏq)
    r,(j′,m′),(j,m) = rme_j′j(o′, Tᵏ, o)
    iszero(r) && return 0
    c = powneg1(Int(j′-m′))*wigner3j(j′, rank(Tᵏ), j,
                                     -m′, component(Tᵏq), m)
    iszero(c) && return 0
    c*r
end

wigner_eckart(o′, lct::LinearCombinationTensor, o) =
    sum(filter!(!iszero, [c*wigner_eckart(o′, Tᵏq, o) for (Tᵏq,c) in lct]))

"""
    dot(o′, T, o)

Calculates the matrix element `⟨o′|T|o⟩` using [`wigner_eckart`](@ref).
"""
LinearAlgebra.dot(o′::SpinOrbital, T::Union{TensorComponent,LinearCombinationTensor}, o::SpinOrbital) =
    wigner_eckart(o′, T, o)

export wigner_eckart
