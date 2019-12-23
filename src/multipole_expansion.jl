# * Multipole expansion

"""
    multipole_expand_scalar_product(a, b, P, Q, c, d)

Multipole-expand the matrix element `⟨ab|P⋅Q|cd⟩`, where the tensor
`P` acts on orbitals `a` & `c`, and the tensor `Q` acts on orbitals
`b` & `d`. The definition is taken from Eq. (13.1.26) of Varshalovich
(1988).
"""
function multipole_expand_scalar_product(a, b, P, Q, c, d, f::Int=1)
    multipole_terms = Pair{Int,Float64}[]

    couples(a, P, c) && couples(b, Q, d) ||
        return multipole_terms

    ks = ranks(a,P,c) ∩ ranks(b,Q,d)

    for k in ks
        X = P(k)⋅Q(k)
        v = dot((a,b), X, (c,d))
        iszero(v) && continue

        push!(multipole_terms, k => f*v)
    end

    multipole_terms
end

"""
    multipole_expand(integral::OrbitalMatrixElement{2,A,<:CoulombInteraction,B})

Multipole-expand the two-body integral resulting from the Coulomb
repulsion between two electrons.
"""
function multipole_expand(integral::OrbitalMatrixElement{2,<:Any,<:CoulombInteraction,<:Any})
    (a,b),g,(c,d) = integral.a,integral.o,integral.b

    terms = NBodyTerm[]

    for (k,coeff) in multipole_expand_scalar_product(a, b, CoulombTensor, CoulombTensor, c, d)
        push!(terms, NBodyTerm([radial_integral(integral.a, (k,g), integral.b)], coeff))
    end

    NBodyMatrixElement(terms)
end

# """
#     multipole_expand(integral::NBodyTermFactor)

# Dummy method that returns `integral` unchanged, used for all
# `NBodyTermFactor`s that are _not_ to be multipole-expanded.
# """
# multipole_expand(integral::NBodyTermFactor) = NBodyMatrixElement([integral])

export multipole_expand_scalar_product, multipole_expand
