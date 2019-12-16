"""
    triangle_range(a,b)

Find all (even) `k` such that `|a-b| ≤ k ≤ a + b`. This is useful when
expanding matrix elements of tensors between angular momenta `a` and
`b` in multipoles `k`; `triangle_range` can then be used to decided
which multipole terms are required.
"""
function triangle_range(a,b)
    kmin = abs(a-b)
    if !iseven(kmin + a + b)
        kmin += 1
    end
    kmin:2:(a+b)
end

"""
    ∏(ℓs...)

Calculates `√((2ℓ₁+1)(2ℓ₂+1)...(2ℓₙ+1))`, which is a common factor in
angular momentum algebra.
"""
function ∏(ℓs...)
    v = 1
    for ℓ in ℓs
        v *= 2ℓ + 1
    end
    √v
end

"""
    powneg1(k) = (-)ᵏ

Calculates powers of negative unity for integer `k`.
"""
powneg1(k::Integer) = isodd(k) ? -1 : 1

"""
    jmⱼ(o::SpinOrbital)

Return the angular momentum and its projection on the z axis of the
spin-orbital `o`.
"""
jmⱼ(o::SpinOrbital{<:Orbital}) = o.orb.ℓ, o.m[1]
jmⱼ(o::SpinOrbital{<:RelativisticOrbital}) = o.orb.j, o.m[1]

"""
    spin(o::SpinOrbital)

Return the spin of the spin-orbital `o`.
"""
spin(o::SpinOrbital{<:Orbital}) = o.m[2]

"""
    @δ (a,b)[, (c,d) ...]

Kronecker ``\\delta_{ab}\\delta_{cd}...`` that tests each pair of
values for equality and quick-returns `0` at the first inequality.
"""
macro δ(vars...)
    code = map(vars) do v
        v.head == :tuple && length(v.args) == 2 ||
            throw(ArgumentError("Expected tuple of 2 variables to compare"))
        a,b = v.args
        :($(esc(a)) == $(esc(b)) || return 0)
    end
    Expr(:block, code..., :(1))
end
