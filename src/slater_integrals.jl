abstract type AbstractRadialIntegral{O} <: Symbolic end
abstract type AbstractSlaterIntegral{O} <: AbstractRadialIntegral{O} end

# * Integrals

struct OverlapIntegral{O} <: AbstractRadialIntegral{O}
    p::Int
    a::O
    b::O
end
OverlapIntegral(p::Integer, a::O, b::O) where {O} =
    OverlapIntegral{O}(p, a, b)
@new_number OverlapIntegral

Base.diff(::OverlapIntegral, orb, occ) =
    error("Not implemented")

struct DiagonalIntegral{O} <: AbstractRadialIntegral{O}
    o::O
end
@new_number DiagonalIntegral

Base.diff(I::DiagonalIntegral{O}, orb::O, occ::II) where {O,II} =
    I.o == orb ? -Sym(:𝓛)*Ket(orb)/occ : 0

struct GeneralSlaterIntegral{O} <: AbstractSlaterIntegral{O}
    k::Int
    a::O
    b::O
    c::O
    d::O
end
@new_number GeneralSlaterIntegral

Base.diff(::GeneralSlaterIntegral, orb, occ) =
    error("Not implemented")

struct DirectSlaterIntegral{O} <: AbstractSlaterIntegral{O}
    k::Int
    a::O
    b::O
end
@new_number DirectSlaterIntegral

isdiagonal(F::DirectSlaterIntegral) = F.a == F.b

function Base.diff(F::DirectSlaterIntegral{O}, orb::O, occ::I) where {k,O,I}
    F.a != orb && F.b != orb && return 0
    other = F.a == orb ? F.b : F.a
    (2/Sym(:r))*(1+isdiagonal(F))/occ*SlaterPotential(F.k, other, other)*Ket(orb)
end

struct ExchangeSlaterIntegral{O} <: AbstractSlaterIntegral{O}
    k::Int
    a::O
    b::O
end
@new_number ExchangeSlaterIntegral

function Base.diff(G::ExchangeSlaterIntegral{O}, orb::O, occ::I) where {k,O,I}
    G.a != orb && G.b != orb && return 0
    other = G.a == orb ? G.b : G.a
    (2/Sym(:r))/occ*SlaterPotential(G.k, orb, other)*Ket(other)
end

struct SlaterPotential{O} <: AbstractRadialIntegral{O}
    k::Int
    a::O
    b::O
end
@new_number SlaterPotential

isdiagonal(Y::SlaterPotential) = Y.a == Y.b
isdirect(Y::SlaterPotential{O}, o::O) where {O} =
    Y.a != o && Y.b != o
isexchange(Y::SlaterPotential{O}, o::O) where {O} =
    Y.a == o || Y.b == o

Base.log(Y::SlaterPotential) = Y.k

# * Pretty-printing

function Base.show(io::IO, o::OverlapIntegral{O}) where O
    write(io, "⟨")
    show(io, o.a)
    write(io, "|")
    show(io, o.b)
    write(io, "⟩")
    o.p != 1 && write(io, to_superscript(o.p))
end

function Base.show(io::IO, I::DiagonalIntegral{O}) where {O}
    write(io, "I(")
    show(io, I.o)
    write(io, ")")
end

function Base.show(io::IO, R::GeneralSlaterIntegral{O}) where {O}
    write(io, "R", to_superscript(R.k), "(")
    show(io, R.a)
    write(io, ", ")
    show(io, R.b)
    write(io, "; ")
    show(io, R.c)
    write(io, ", ")
    show(io, R.d)
    write(io, ")")
end

function Base.show(io::IO, F::DirectSlaterIntegral{O}) where {O}
    write(io, "F", to_superscript(F.k), "(")
    show(io, F.a)
    if F.a != F.b
        write(io, ", ")
        show(io, F.b)
    end
    write(io, ")")
end

function Base.show(io::IO, G::ExchangeSlaterIntegral{O}) where {O}
    write(io, "G", to_superscript(G.k), "(")
    show(io, G.a)
    write(io, ", ")
    show(io, G.b)
    write(io, ")")
end

function Base.show(io::IO, Y::SlaterPotential{O}) where {O}
    write(io, "Y", to_superscript(Y.k), "(")
    show(io, Y.a)
    write(io, ", ")
    show(io, Y.b)
    write(io, ")")
end

# ** LaTeX

import Symbolics: latex, delimit

function Symbolics.latex(o::OverlapIntegral{O}) where O
    la,mda = latex(o.a)
    lb,mdb = latex(o.b)
    lp,mdp = latex(o.p)
    "\\langle$(la)|$(lb)\\rangle" * (o.p != 1 ? "^{$(lp)}" : ""),max(mda,mdb,mdp)
end

function Symbolics.latex(I::DiagonalIntegral{O}) where {O}
    lo,mdo = delimit(latex(I.o)...)
    "I$(lo)",mdo
end

function Symbolics.latex(R::GeneralSlaterIntegral{O}) where {O}
    la,mda = latex(R.a)
    lb,mdb = latex(R.b)
    lc,mdc = latex(R.c)
    ld,mdd = latex(R.d)
    largs,md = delimit("$(la),$(lb);$(lc),$(ld)",max(mda,mdb,mdc,mdd))
    "R^{$(R.k)}$(largs)",md
end

function Symbolics.latex(F::DirectSlaterIntegral{O}) where {O}
    la,mda = latex(F.a)
    lb,mdb = latex(F.b)
    largs,md = if isdiagonal(F)
        delimit(la,mda)
    else
        delimit("$(la),$(lb)", max(mda,mdb))
    end
    "F^{$(F.k)}$(largs)",md
end

function Symbolics.latex(G::ExchangeSlaterIntegral{O}) where {O}
    la,mda = latex(G.a)
    lb,mdb = latex(G.b)
    largs,md = delimit("$(la),$(lb)",max(mda,mdb))
    "G^{$(G.k)}$(largs)",md
end

function Symbolics.latex(Y::SlaterPotential{O}) where {O}
    la,mda = latex(Y.a)
    lb,mdb = latex(Y.b)
    largs,md = if isdiagonal(Y)
        delimit(la,mda)
    else
        delimit("$(la),$(lb)", max(mda,mdb))
    end
    "Y^{$(Y.k)}$(largs)",md
end

export OverlapIntegral, DiagonalIntegral, GeneralSlaterIntegral, DirectSlaterIntegral, ExchangeSlaterIntegral, SlaterPotential
