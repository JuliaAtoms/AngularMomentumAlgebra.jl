using Symbolics
import AngularMomentumAlgebra: @new_number

# * Fock operator

struct Fock <: Symbolic end
Base.:(==)(::Fock, ::Fock) = true
Base.show(io::IO, ::Fock) = write(io, "f̂")

# * OneBodyHamiltonian operator

struct OneBodyHamiltonian <: Symbolic end
Base.:(==)(::OneBodyHamiltonian, ::OneBodyHamiltonian) = true
Base.show(io::IO, ::OneBodyHamiltonian) = write(io, "ĥ")

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

Base.diff(I::OneBodyIntegral{A,B}, orb::O, occ::II=1) where {A,B,O,II} =
    I.b == orb ? OneBodyHamiltonian()*Bra(I.a)/occ : 0

Base.diff(I::OneBodyIntegral{A,B}, corb::Conjugate{O}, occ::II=1) where {A,B,O,II} =
    I.a == corb.orbital ? OneBodyHamiltonian()*Ket(I.b)/occ : 0

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

struct DirectExchangeIntegral{A,B}
    a::A
    b::B
end
Base.:(==)(a::DirectExchangeIntegral, b::DirectExchangeIntegral) =
    a.a == b.a && a.b == b.b
isdiagonal(FG::DirectExchangeIntegral) = FG.a == FG.b

function Base.diff(FG::ExchangeIntegral{A,B}, orb::O, occ::I=1) where {A,B,O,I}
    orb ∉ [FG.a, FG.b] && return 0
    other = FG.a == orb ? FG.b : FG.a
    inv(occ)*(DirectPotential(other,other)*Bra(orb)-
              ExchangePotential(other,other)*Bra(orb))
end


function Base.diff(FG::DirectExchangeIntegral{A,B}, corb::Conjugate{O}, occ::I=1) where {A,B,O,I}
    corb.orbital ∉ [FG.a, FG.b] && return 0
    other = FG.a == corb.orbital ? FG.b : FG.a
    inv(occ)*(DirectPotential(other,other)*Ket(corb.orbital)-
              ExchangePotential(other,other)*Ket(corb.orbital))
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

function Base.show(io::IO, FG::DirectExchangeIntegral)
    write(io, "[")
    show(io, FG.a)
    write(io, " ")
    show(io, FG.b)
    write(io, "||")
    show(io, FG.a)
    write(io, " ")
    show(io, FG.b)
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

@new_number Fock
@new_number OneBodyHamiltonian
@new_number OneBodyIntegral
@new_number RepulsionPotential
@new_number GeneralRepulsionIntegral
@new_number DirectIntegral
@new_number ExchangeIntegral


export Fock, OneBodyHamiltonian, OneBodyIntegral,
    RepulsionPotential, DirectPotential, ExchangePotential
    GeneralRepulsionIntegral, DirectIntegral, ExchangeIntegral, DirectExchangeIntegral
