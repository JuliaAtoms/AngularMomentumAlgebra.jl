# * Repulsion integrals

abstract type AbstractRepulsionIntegral{A,B,C,D} end

# ** General

struct GeneralRepulsionIntegral{A,B,C,D} <: AbstractRepulsionIntegral{A,B,C,D}
    a::A
    b::B
    c::C
    d::D
end
Base.:(==)(a::GeneralRepulsionIntegral, b::GeneralRepulsionIntegral) =
    a.a == b.a && a.b == b.b && a.c == b.c && a.d == b.d

function Base.diff(R::GeneralRepulsionIntegral{A,B,C,D}, orb::O) where {A,B,C,D,O}
    orb ∉ [R.c, R.d] && return zero(RepulsionPotential{:general})
    orb == R.c ? RepulsionPotential{:general,B,D,Conjugate{A}}(R.b,R.d,conj(R.a)) :
        RepulsionPotential{:general,A,C,Conjugate{B}}(R.a,R.c,conj(R.b))
end

function Base.diff(R::GeneralRepulsionIntegral{A,B,C,D}, corb::Conjugate{O}) where {A,B,C,D,O}
    corb.orbital ∉ [R.a, R.b] && return zero(RepulsionPotential{:general})
    corb.orbital == R.a ? RepulsionPotential{:general,B,D,C}(R.b,R.d,R.c) :
        RepulsionPotential{:general,A,C,D}(R.a,R.c,R.d)
end

# ** Direct

struct DirectIntegral{A,B} <: AbstractRepulsionIntegral{A,B,A,B}
    a::A
    b::B
end
Base.:(==)(a::DirectIntegral, b::DirectIntegral) =
    a.a == b.a && a.b == b.b

isdiagonal(F::DirectIntegral) = F.a == F.b

function Base.diff(F::DirectIntegral{A,B}, orb::O) where {A,B,O}
    orb ∉ [F.a, F.b] && return zero(DirectPotential)
    other = F.a == orb ? F.b : F.a
    DirectPotential(other,other,conj(orb))
end

function Base.diff(F::DirectIntegral{A,B}, corb::Conjugate{O}) where {A,B,O}
    corb.orbital ∉ [F.a, F.b] && return zero(DirectPotential)
    other = F.a == corb.orbital ? F.b : F.a
    DirectPotential(other,other,corb.orbital)
end

# ** Exchange

struct ExchangeIntegral{A,B} <: AbstractRepulsionIntegral{A,B,B,A}
    a::A
    b::B
end
Base.:(==)(a::ExchangeIntegral, b::ExchangeIntegral) =
    a.a == b.a && a.b == b.b

function Base.diff(G::ExchangeIntegral{A,B}, orb::O) where {A,B,O}
    orb ∉ [G.a, G.b] && return zero(ExchangePotential)
    other = G.a == orb ? G.b : G.a
    ExchangePotential(other,other,conj(orb))
end

function Base.diff(G::ExchangeIntegral{A,B}, corb::Conjugate{O}) where {A,B,O}
    corb.orbital ∉ [G.a, G.b] && return zero(ExchangePotential)
    other = G.a == corb.orbital ? G.b : G.a
    ExchangePotential(other,other,corb.orbital)
end

Base.:(==)(a::DirectIntegral, b::ExchangeIntegral) =
    isdiagonal(a) && a.a == b.a && a.b == b.b
Base.:(==)(a::ExchangeIntegral, b::DirectIntegral) =
    b == a

# ** Two-body integral

"""
    TwoBodyIntegral(a,b,c,d)

The two-body integral `[ab||cd]` is short-hand notation for
`[ab|cd]-[ab|dc]`, which is useful since that pair is very commonly
seen together.
"""
struct TwoBodyIntegral{A,B,C,D} <: AbstractRepulsionIntegral{A,B,C,D}
    a::A
    b::B
    c::C
    d::D
    TwoBodyIntegral{A,B,C,D}(a::A, b::B, c::C=a, d::D=b) where {A,B,C,D} =
        new{A,B,C,D}(a, b, c, d)
    TwoBodyIntegral(a::A, b::B, c::C=a, d::D=b) where {A,B,C,D} =
        new{A,B,C,D}(a, b, c, d)
end
Base.:(==)(a::TwoBodyIntegral, b::TwoBodyIntegral) =
    a.a == b.a && a.b == b.b && a.c == b.c && a.d == b.d
isdiagonal(FG::TwoBodyIntegral) = FG.a == FG.b == FG.c == FG.d

function Base.diff(FG::TwoBodyIntegral{A,B,C,D}, orb::O) where {A,B,C,D,O,I}
    orb ∉ [FG.c, FG.d] && return zero(DirectExchangePotentials)
    other,pot_orbs = FG.d == orb ? (FG.b,(FG.a,FG.c)) : (FG.a,(FG.b,FG.d))
    DirectExchangePotentials(pot_orbs..., conj(other))
end

function Base.diff(FG::TwoBodyIntegral{A,B,C,D}, corb::Conjugate{O}) where {A,B,C,D,O,I}
    corb.orbital ∉ [FG.a, FG.b] && return zero(DirectExchangePotentials)
    other,pot_orbs = FG.b == corb.orbital ? (FG.d,(FG.a,FG.c)) : (FG.c,(FG.b,FG.d))
    DirectExchangePotentials(pot_orbs..., other)
end

# ** Pretty-printing

function Base.show(io::IO, R::GeneralRepulsionIntegral)
    write(io, "[")
    show(io, R.a)
    write(io, " ")
    show(io, R.b)
    write(io, "|")
    show(io, R.c)
    write(io, " ")
    show(io, R.d)
    write(io, "]")
end

function Base.show(io::IO, F::DirectIntegral)
    write(io, "F(")
    show(io, F.a)
    write(io, ", ")
    show(io, F.b)
    write(io, ")")
end

function Base.show(io::IO, G::ExchangeIntegral)
    write(io, "G(")
    show(io, G.a)
    write(io, ", ")
    show(io, G.b)
    write(io, ")")
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

export GeneralRepulsionIntegral, DirectIntegral, ExchangeIntegral, TwoBodyIntegral,
    isdiagonal
