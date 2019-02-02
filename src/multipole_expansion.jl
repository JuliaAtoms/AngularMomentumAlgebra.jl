# * Multipole expansion

#=

Equations (3.42–44) of

- Lindgren, I. (1986). Atomic many-body theory. Berlin New York:
  Springer-Verlag.

\[\langle ab|r_{12}^{-1}|cd\rangle =
\delta(m^a_s,m^c_s)\delta(m^b_s,m^d_s)
\sum_k X^k(ab,cd)\]
where
\[X^k(ab,cd)=X(k,\ell_a\ell_b\ell_c\ell_d)R^k(ab,cd),\]
\[X(k,\ell_a\ell_b\ell_c\ell_d)=
(-)^k
\langle\ell_a||\mathbf{C}^{(k)}||\ell_c\rangle
\langle\ell_b||\mathbf{C}^{(k)}||\ell_d\rangle,\]
and
\[R^k(ab,cd)=
\iint\mathrm{d}{r_1}\mathrm{d}r_2
P_a^*(r_1)P_b^*(r_2)
\frac{r^k_<}{r^{k+1}_>}
P_c(r_1)P_d(r_2).\]

The two-body integral type is defined as
\[[ab||cd] \equiv
\langle ab|r_{12}^{-1}|cd\rangle
-\langle ab|r_{12}^{-1}|db\rangle.\]

The reduced matrix element is given by [Eq. (2.127), ibid.]:
\[\langle\ell||\mathbf{C}^{(k)}||\ell'\rangle \equiv
(-)^\ell\sqrt{(2\ell+1)(2\ell'+1)}
\begin{pmatrix}\ell&k&\ell'\\0&0&0\end{pmatrix}.\]
=#

function triangle_range(a,b)
    kmin = abs(a-b)
    if !iseven(kmin + a + b)
        kmin += 1
    end
    kmin:2:(a+b)
end
rme(ℓ,k,ℓ′) = (-1)^ℓ*√((2ℓ+1)*(2ℓ′+1))*wigner3j(ℓ,k,ℓ′,0,0,0)

function multipole_expand(a::SpinOrbital, b::SpinOrbital, c::SpinOrbital, d::SpinOrbital,
                          f::Int=1)
    multipole_terms = Pair{Int,Float64}[]
    a.spin == c.spin && b.spin == d.spin || return multipole_terms

    ℓa,ℓb,ℓc,ℓd = a.orb.ℓ,b.orb.ℓ,c.orb.ℓ,d.orb.ℓ
    ks = triangle_range(ℓa,ℓc) ∩ triangle_range(ℓb,ℓd)
    for k in ks
        push!(multipole_terms, k => f*rme(ℓa,k,ℓc)*rme(ℓb,k,ℓd))
    end

    multipole_terms
end

function multipole_expand(integral::TwoBodyIntegral)
    a,b,c,d = integral.a,integral.b,integral.c,integral.d

    direct_terms = multipole_expand(a,b,c,d)
    exchange_terms = multipole_expand(a,b,d,c,-1)

    direct_terms,exchange_terms
end

export multipole_expand
