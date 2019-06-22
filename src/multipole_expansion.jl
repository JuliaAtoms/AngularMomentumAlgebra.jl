# * Multipole expansion

"""
    multipole_expand_scalar_product(a, b, P, Q, c, d)

Multipole expand the matrix element `⟨ab|P⋅Q|cd⟩`, where the tensor
`P` acts on orbitals `a` & `c`, and the tensor `Q` acts on orbitals
`b` & `d`. The definition is taken from Eq. (13.1.26) of Varshalovich
(1988).
"""
function multipole_expand_scalar_product(a, b, ::Type{P}, ::Type{Q}, c, d,
                                         f::Int=1) where {P<:Tensor,Q<:Tensor}
    multipole_terms = Pair{Int,Float64}[]

    couples_spin(spin(a), P, spin(c)) && couples_spin(spin(b), Q, spin(d)) ||
        return multipole_terms

    ja,ma = jmⱼ(a)
    jb,mb = jmⱼ(b)
    jc,mc = jmⱼ(c)
    jd,md = jmⱼ(d)

    ks = triangle_range(ja,jc) ∩ triangle_range(jb,jd)

    for k in ks
        α = ma-mc
        (α != md-mb || abs(α) > k) && continue

        v = inv(∏(ja,jb)) * powneg1(-α) *
            clebsch_gordan_condon_shortley(jc, mc, k, α, ja, ma) *
            clebsch_gordan_condon_shortley(jd, md, k, -α, jb, mb) *
            rme(ja, P(k), jc) *
            rme(jb, Q(k), jd)

        push!(multipole_terms, k => f*v)
    end

    multipole_terms
end

# function multipole_expand_coulomb(a::SpinOrbital, b::SpinOrbital, c::SpinOrbital, d::SpinOrbital,
#                                   f::Int=1)
#     multipole_terms = Pair{Int,Float64}[]
#     a.spin == c.spin && b.spin == d.spin || return multipole_terms

#     ℓa,ℓb,ℓc,ℓd = a.orb.ℓ,b.orb.ℓ,c.orb.ℓ,d.orb.ℓ
#     ma,mb,mc,md = a.mℓ,b.mℓ,c.mℓ,d.mℓ

#     ks = triangle_range(ℓa,ℓc) ∩ triangle_range(ℓb,ℓd)

#     for k in ks
#         # v = 0
#         # # println("k = $k")
#         # for α = -k:k
#         #     v += ((isodd(-α + ma + mb) ? -1 : 1) *
#         #           wigner3j(ℓc, k, ℓa,
#         #                    mc, α, -ma) *
#         #           wigner3j(ℓd, k, ℓb,
#         #                    md, -α, -mb))
#         #     # v += (isodd(α) ? -1 : 1) *
#         #     #     clebschgordan(ℓc,mc, k,α, ℓa,ma) *
#         #     #     clebschgordan(ℓd,md, k,-α, ℓb,mb)
#         #     # println((α, v))
#         # end
#         # # v *= rme_spherical(ℓa, k, ℓc)*rme_spherical(ℓb, k, ℓd)
#         # v *= (∏(ℓa,ℓb,ℓc,ℓd)*
#         #       wigner3j(ℓc, k, ℓa,
#         #                0,  0, 0)*
#         #       wigner3j(ℓd, k, ℓb,
#         #                0,  0, 0))
#         # iszero(v) || push!(multipole_terms, k => f*v)

#         α = ma-mc
#         (α != md-mb || abs(α) > k) && continue

#         push!(multipole_terms, k => f*((-1)^(mb+mc)*∏(ℓa,ℓb,ℓc,ℓd)*
#                                       wigner3j(ℓc, k, ℓa,
#                                                mc, α, -ma)*
#                                       wigner3j(ℓd, k, ℓb,
#                                                md, -α, -mb)*
#                                       wigner3j(ℓc, k, ℓa,
#                                                0,  0, 0)*
#                                       wigner3j(ℓd, k, ℓb,
#                                                0,  0, 0)))
#     end

#     multipole_terms
# end

"""
    multipole_expand(integral::OrbitalMatrixElement{2,A,<:CoulombInteraction,B})

Multipole-expand the two-body integral resulting from the Coulomb
repulsion between two electrons.
"""
function multipole_expand(integral::OrbitalMatrixElement{2,A,<:CoulombInteraction,B}) where {A,B}
    (a,b),g,(c,d) = integral.a,integral.o,integral.b

    terms = NBodyTerm[]

    # for (k,coeff) in multipole_expand_coulomb(a,b,c,d)
    for (k,coeff) in multipole_expand_scalar_product(a,b, SphericalTensor, SphericalTensor, c,d)
        push!(terms, NBodyTerm([OrbitalMatrixElement(A[a,b],
                                                     CoulombInteractionMultipole(k,g),
                                                     B[c,d])], coeff))
    end

    NBodyMatrixElement(terms)
end

"""
    multipole_expand(integral::NBodyTermFactor)

Dummy method that returns `integral` unchanged, used for all
`NBodyTermFactor`s that are _not_ to be multipole-expanded.
"""
multipole_expand(integral::NBodyTermFactor) = NBodyMatrixElement([integral])

export multipole_expand_scalar_product, multipole_expand
