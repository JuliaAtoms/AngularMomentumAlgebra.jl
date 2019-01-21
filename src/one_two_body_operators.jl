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
    orb ∉ [FG.c, FG.d] && return 0
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
    corb.orbital ∉ [FG.a, FG.b] && return 0
    other,pot_orbs = FG.b == corb.orbital ? (FG.d,(FG.a,FG.c)) : (FG.c,(FG.b,FG.d))
    DirectExchangePotentials(pot_orbs..., Ket(other))
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
GeneralRepulsionIntegral, DirectIntegral, ExchangeIntegral, TwoBodyIntegral
