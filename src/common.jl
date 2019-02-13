"""
    triangle_range(a,b)

Find all (even) `k` such that `|a-b| ≤ k ≤ a + b`. This is useful when
expanding matrix elements of tensors between angular momenta `a` and
`b` in multipoles `k`; `triangl_range` can then be used to decided
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
jmⱼ(o::SpinOrbital) = o.orb.ℓ, o.mℓ

"""
    spin(o::SpinOrbital)

Return the spin of the spin-orbital `o`.
"""
spin(o::SpinOrbital) = o.spin
