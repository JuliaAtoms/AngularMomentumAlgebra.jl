using Symbolics
import AngularMomentumAlgebra: @new_number

# * OneBodyHamiltonian operator

struct OneBodyHamiltonian{O} <: Symbolic
    orb::O
end
Base.:(==)(a::OneBodyHamiltonian, b::OneBodyHamiltonian) = a.orb == b.orb
function Base.show(io::IO, h::OneBodyHamiltonian)
    h.orb isa Bra && show(io, h.orb)
    write(io, "ĥ")
    h.orb isa Bra || show(io, h.orb)
end

Base.iszero(h::OneBodyHamiltonian) = h.orb == 0

# * One-body integral

struct OneBodyIntegral{A,B} <: Symbolic
    a::A
    b::B
end
Base.:(==)(a::OneBodyIntegral, b::OneBodyIntegral) = a.a == b.a && a.b == b.b
isdiagonal(I::OneBodyIntegral) = I.a == I.b

function Base.show(io::IO, I::OneBodyIntegral)
    write(io, "I(")
    show(io, I.a)
    if !isdiagonal(I)
        write(io, ", ")
        show(io, I.b)
    end
    write(io, ")")
end

Base.diff(I::OneBodyIntegral{A,B}, orb::O) where {A,B,O,II} =
    OneBodyHamiltonian(I.b == orb ? Bra(I.a) : 0)

Base.diff(I::OneBodyIntegral{A,B}, orb::O, occ::II) where {A,B,O,II} =
    diff(I, orb)/occ

Base.diff(I::OneBodyIntegral{A,B}, corb::Conjugate{O}) where {A,B,O,II} =
    OneBodyHamiltonian(I.a == corb.orbital ? Ket(I.b) : 0)

Base.diff(I::OneBodyIntegral{A,B}, corb::Conjugate{O}, occ::II) where {A,B,O,II} =
    diff(I, corb)/occ

# * Repulsion potentials

struct RepulsionPotential{kind,A,B} <: Symbolic
    a::A
    b::B
end
RepulsionPotential{kind}(a::A, b::B) where {kind,A,B} =
    RepulsionPotential{kind,A,B}(a,b)

Base.:(==)(a::RepulsionPotential{kind}, b::RepulsionPotential{kind′}) where {kind,kind′} =
    kind == kind′ && a.a == b.a && a.b == b.b ||
    a.a == a.b == b.a == b.b

const DirectPotential{A,B} = RepulsionPotential{:direct,A,B}
const ExchangePotential{A,B} = RepulsionPotential{:exchange,A,B}

function Base.show(io::IO, J::DirectPotential)
    write(io, "Ĵ{")
    show(io, J.a)
    write(io, ";")
    show(io, J.b)
    write(io, "}")
end

function Base.show(io::IO, K::ExchangePotential)
    write(io, "K̂{")
    show(io, K.a)
    write(io, ";")
    show(io, K.b)
    write(io, "}")
end

struct DirectExchangePotentials{A,B,O}
    a::A
    b::B
    o::O # The orbital o feels the potentials formed by a & b
end
Base.:(==)(a::DirectExchangePotentials,b::DirectExchangePotentials) =
    a.a == b.a && a.b == b.b && a.o == b.o

Base.zero(::Type{DirectExchangePotentials}) =
    DirectExchangePotentials(0,0,0)

Base.iszero(JK::DirectExchangePotentials) =
    JK.a == JK.b == JK.o == 0

isdiagonal(dxp::DirectExchangePotentials) =
    dxp.a == dxp.b

function Base.show(io::IO, JK::DirectExchangePotentials)
    write(io, "[Ĵ{")
    show(io, JK.a)
    write(io, ";")
    show(io, JK.b)
    write(io, "}")
    write(io, " - ")
    write(io, "K̂{")
    show(io, JK.a)
    write(io, ";")
    show(io, JK.b)
    write(io, "}]")
    show(io, JK.o)
end

# * Repulsion integrals

abstract type AbstractRepulsionIntegral{A,B,C,D} <: Symbolic end

struct GeneralRepulsionIntegral{A,B,C,D} <: AbstractRepulsionIntegral{A,B,C,D}
    a::A
    b::B
    c::C
    d::D
end
Base.:(==)(a::GeneralRepulsionIntegral, b::GeneralRepulsionIntegral) =
    a.a == b.a && a.b == b.b && a.c == b.c && a.d == b.d

Base.diff(::GeneralRepulsionIntegral, orb, occ) =
    error("Not implemented")

struct DirectIntegral{A,B} <: AbstractRepulsionIntegral{A,B,A,B}
    a::A
    b::B
end
Base.:(==)(a::DirectIntegral, b::DirectIntegral) =
    a.a == b.a && a.b == b.b

isdiagonal(F::DirectIntegral) = F.a == F.b

function Base.diff(F::DirectIntegral{A,B}, orb::O, occ::I=1) where {A,B,O,I}
    orb ∉ [F.a, F.b] && return 0
    other = F.a == orb ? F.b : F.a
    inv(occ)*DirectPotential(other,other)*Bra(orb)
end

function Base.diff(F::DirectIntegral{A,B}, corb::Conjugate{O}, occ::I=1) where {A,B,O,I}
    corb.orbital ∉ [F.a, F.b] && return 0
    other = F.a == corb.orbital ? F.b : F.a
    inv(occ)*DirectPotential(other,other)*Ket(corb.orbital)
end

struct ExchangeIntegral{A,B} <: AbstractRepulsionIntegral{A,B,B,A}
    a::A
    b::B
end
Base.:(==)(a::ExchangeIntegral, b::ExchangeIntegral) =
    a.a == b.a && a.b == b.b

function Base.diff(G::ExchangeIntegral{A,B}, orb::O, occ::I=1) where {A,B,O,I}
    orb ∉ [G.a, G.b] && return 0
    other = G.a == orb ? G.b : G.a
    inv(occ)*ExchangePotential(other,other)*Bra(orb)
end

function Base.diff(G::ExchangeIntegral{A,B}, corb::Conjugate{O}, occ::I=1) where {A,B,O,I}
    corb.orbital ∉ [G.a, G.b] && return 0
    other = G.a == corb.orbital ? G.b : G.a
    inv(occ)*ExchangePotential(other,other)*Ket(corb.orbital)
end

Base.:(==)(a::DirectIntegral, b::ExchangeIntegral) =
    isdiagonal(a) && a.a == b.a && a.b == b.b
Base.:(==)(a::ExchangeIntegral, b::DirectIntegral) =
    b == a

struct TwoBodyIntegral{A,B,C,D} <: AbstractRepulsionIntegral{A,B,C,D}
    a::A
    b::B
    c::C
    d::D
    TwoBodyIntegral{A,B,C,D}(a::A, b::B, c::C=a, d::D=b) where {A,B,C,D} =
        new{A,B,C,D}(a, b, c, d)
end
Base.:(==)(a::TwoBodyIntegral, b::TwoBodyIntegral) =
    a.a == b.a && a.b == b.b && a.c == b.c && a.d == b.d
isdiagonal(FG::TwoBodyIntegral) = FG.a == FG.b == FG.c == FG.d

function Base.diff(FG::TwoBodyIntegral{A,B,C,D}, orb::O, occ::I=1) where {A,B,C,D,O,I}
    orb ∉ [FG.c, FG.d] && return 0
    other,pot_orbs = FG.d == orb ? (FG.b,(FG.a,FG.c)) : (FG.a,(FG.b,FG.d))
    inv(occ)*(DirectPotential(pot_orbs...)*Bra(other)-
              ExchangePotential(pot_orbs...)*Bra(other))
end

function Base.diff(FG::TwoBodyIntegral{A,B,C,D}, orb::O) where {A<:SpinOrbital,B<:SpinOrbital,C<:SpinOrbital,D<:SpinOrbital,O<:SpinOrbital,I}
    orb ∉ [FG.c, FG.d] && return zero(DirectExchangePotentials)
    other,pot_orbs = FG.d == orb ? (FG.b,(FG.a,FG.c)) : (FG.a,(FG.b,FG.d))
    DirectExchangePotentials(pot_orbs..., Bra(other))
end

function Base.diff(FG::TwoBodyIntegral{A,B,C,D}, corb::Conjugate{O}, occ::I=1) where {A,B,C,D,O,I}
    corb.orbital ∉ [FG.a, FG.b] && return 0
    other,pot_orbs = FG.b == corb.orbital ? (FG.d,(FG.a,FG.c)) : (FG.c,(FG.b,FG.d))
    inv(occ)*(DirectPotential(pot_orbs...)*Ket(other)-
              ExchangePotential(pot_orbs...)*Ket(other))
end

function Base.diff(FG::TwoBodyIntegral{A,B,C,D}, corb::Conjugate{O}) where {A<:SpinOrbital,B<:SpinOrbital,C<:SpinOrbital,D<:SpinOrbital,O<:SpinOrbital,I}
    corb.orbital ∉ [FG.a, FG.b] && return zero(DirectExchangePotentials)
    other,pot_orbs = FG.b == corb.orbital ? (FG.d,(FG.a,FG.c)) : (FG.c,(FG.b,FG.d))
    DirectExchangePotentials(pot_orbs..., Ket(other))
end

# ** Multipole expansion

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

# ** Pretty-printing

function Base.show(io::IO, R::GeneralRepulsionIntegral)
    write(io, "⟨")
    show(io, R.a)
    write(io, " ")
    show(io, R.b)
    write(io, "|")
    show(io, R.c)
    write(io, " ")
    show(io, R.d)
    write(io, "⟩")
end

function Base.show(io::IO, F::DirectIntegral)
    write(io, "⟨")
    show(io, F.a)
    write(io, " ")
    show(io, F.b)
    write(io, "|")
    show(io, F.a)
    write(io, " ")
    show(io, F.b)
    write(io, "⟩")
end

function Base.show(io::IO, G::ExchangeIntegral)
    write(io, "⟨")
    show(io, G.a)
    write(io, " ")
    show(io, G.b)
    write(io, "|")
    show(io, G.b)
    write(io, " ")
    show(io, G.a)
    write(io, "⟩")
end

function Base.show(io::IO, FG::TwoBodyIntegral)
    write(io, "[")
    show(io, FG.a)
    write(io, " ")
    show(io, FG.b)
    write(io, "||")
    show(io, FG.c)
    write(io, " ")
    show(io, FG.d)
    write(io, "]")
end

# ** LaTeX

function Symbolics.latex(R::GeneralRepulsionIntegral)
    la,mda = Symbolics.latex(R.a)
    lb,mdb = Symbolics.latex(R.b)
    lc,mdc = Symbolics.latex(R.c)
    ld,mdd = Symbolics.latex(R.d)
    "[$(la),$(lb)|$(lc),$(ld)]", max(mda,mdb,mdc,mdd)
end

function Symbolics.latex(F::DirectIntegral)
    la,mda = Symbolics.latex(F.a)
    lb,mdb = Symbolics.latex(F.b)
    "[$(la),$(lb)|$(la),$(lb)]", max(mda,mdb)
end

function Symbolics.latex(G::ExchangeIntegral)
    la,mda = Symbolics.latex(G.a)
    lb,mdb = Symbolics.latex(G.b)
    "[$(la),$(lb)|$(lb),$(la)]", max(mda,mdb)
end

# * Symbolics registration

@new_number OneBodyHamiltonian
@new_number OneBodyIntegral
@new_number RepulsionPotential
@new_number GeneralRepulsionIntegral
@new_number DirectIntegral
@new_number ExchangeIntegral


export OneBodyHamiltonian, OneBodyIntegral,
    RepulsionPotential, DirectPotential, ExchangePotential, DirectExchangePotentials,
    GeneralRepulsionIntegral, DirectIntegral, ExchangeIntegral, TwoBodyIntegral,
    isdiagonal, multipole_expand
