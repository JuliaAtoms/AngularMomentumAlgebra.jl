abstract type AbstractRadialIntegral{O} <: Symbolic end
abstract type AbstractSlaterIntegral{k, O} <: AbstractRadialIntegral{O} end

# * Integrals

struct OverlapIntegral{O} <: AbstractRadialIntegral{O}
    p::Int
    a::O
    b::O
end
OverlapIntegral(p::Integer, a::O, b::O) where {O} =
    OverlapIntegral{O}(p, a, b)
Base.:(==)(a::OverlapIntegral, b::OverlapIntegral) =
    a.p == b.p && a.a == b.a && a.b == b.b

Base.diff(::OverlapIntegral, orb, occ) =
    error("Not implemented")

struct DiagonalIntegral{O} <: AbstractRadialIntegral{O}
    o::O
end
Base.:(==)(a::DiagonalIntegral, b::DiagonalIntegral) = a.o == b.o

Base.diff(I::DiagonalIntegral{O}, orb::O, occ::II) where {O,II} =
    I.o == orb ? -Sym(:𝓛)*Ket(orb)/occ : 0

struct GeneralSlaterIntegral{k,O} <: AbstractSlaterIntegral{k,O}
    a::O
    b::O
    c::O
    d::O
end
GeneralSlaterIntegral(k::Integer, a::O, b::O, c::O, d::O) where {O} =
    GeneralSlaterIntegral{k,O}(a, b, c, d)
Base.:(==)(a::GeneralSlaterIntegral{k}, b::GeneralSlaterIntegral{k′}) where {k,k′} =
    k == k′ && a.a == b.a && a.b == b.b && a.c == b.c && a.d == b.d

Base.diff(::GeneralSlaterIntegral, orb, occ) =
    error("Not implemented")

struct DirectSlaterIntegral{k,O} <: AbstractSlaterIntegral{k,O}
    a::O
    b::O
end
DirectSlaterIntegral(k::Integer, a::O, b::O) where {O} =
    DirectSlaterIntegral{k,O}(a, b)
Base.:(==)(a::DirectSlaterIntegral{k}, b::DirectSlaterIntegral{k′}) where {k,k′} =
    k == k′ && a.a == b.a && a.b == b.b

isdiagonal(F::DirectSlaterIntegral) = F.a == F.b

function Base.diff(F::DirectSlaterIntegral{k,O}, orb::O, occ::I) where {k,O,I}
    F.a != orb && F.b != orb && return 0
    other = F.a == orb ? F.b : F.a
    (2/Sym(:r))*(1+isdiagonal(F))/occ*SlaterPotential(k, other, other)*Ket(orb)
end

struct ExchangeSlaterIntegral{k,O} <: AbstractSlaterIntegral{k,O}
    a::O
    b::O
end
ExchangeSlaterIntegral(k::Integer, a::O, b::O) where {O} =
    ExchangeSlaterIntegral{k,O}(a, b)
Base.:(==)(a::ExchangeSlaterIntegral{k}, b::ExchangeSlaterIntegral{k′}) where {k,k′} =
    k == k′ && a.a == b.a && a.b == b.b

function Base.diff(G::ExchangeSlaterIntegral{k,O}, orb::O, occ::I) where {k,O,I}
    G.a != orb && G.b != orb && return 0
    other = G.a == orb ? G.b : G.a
    (2/Sym(:r))/occ*SlaterPotential(k, orb, other)*Ket(other)
end

struct SlaterPotential{k,O} <: AbstractRadialIntegral{O}
    a::O
    b::O
end
SlaterPotential(k::Integer, a::O, b::O) where {O} =
    SlaterPotential{k,O}(a, b)
Base.:(==)(a::SlaterPotential{k}, b::SlaterPotential{k′}) where {k,k′} =
    k == k′ && a.a == b.a && a.b == b.b

isdiagonal(Y::SlaterPotential) = Y.a == Y.b
isdirect(Y::SlaterPotential{k,O}, o::O) where {k,O} =
    Y.a != o && Y.b != o
isexchange(Y::SlaterPotential{k,O}, o::O) where {k,O} =
    Y.a == o || Y.b == o

Base.log(::SlaterPotential{k}) where k = k

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

function Base.show(io::IO, R::GeneralSlaterIntegral{k,O}) where {k,O}
    write(io, "R", to_superscript(k), "(")
    show(io, R.a)
    write(io, ", ")
    show(io, R.b)
    write(io, "; ")
    show(io, R.c)
    write(io, ", ")
    show(io, R.d)
    write(io, ")")
end

function Base.show(io::IO, F::DirectSlaterIntegral{k,O}) where {k,O}
    write(io, "F", to_superscript(k), "(")
    show(io, F.a)
    if F.a != F.b
        write(io, ", ")
        show(io, F.b)
    end
    write(io, ")")
end

function Base.show(io::IO, G::ExchangeSlaterIntegral{k, O}) where {k,O}
    write(io, "G", to_superscript(k), "(")
    show(io, G.a)
    write(io, ", ")
    show(io, G.b)
    write(io, ")")
end

function Base.show(io::IO, Y::SlaterPotential{k, O}) where {k,O}
    write(io, "Y", to_superscript(k), "(")
    show(io, Y.a)
    write(io, ", ")
    show(io, Y.b)
    write(io, ")")
end

# ** LaTeX

function Symbolics.latex(o::Orbital)
    "$(o.n)\\mathrm{$(spectroscopic_label(o.ℓ))}",0
end

function Symbolics.latex(o::OverlapIntegral{O}) where O
    la,mda = Symbolics.latex(o.a)
    lb,mdb = Symbolics.latex(o.b)
    lp,mdp = Symbolics.latex(o.p)
    "\\langle$(la)|$(lb)\\rangle" * (o.p != 1 ? "^{$(lp)}" : ""),max(mda,mdb,mdp)
end

function Symbolics.latex(I::DiagonalIntegral{O}) where {O}
    lo,mdo = Symbolics.latex(I.o)
    "I($(lo))",mdo
end

function Symbolics.latex(R::GeneralSlaterIntegral{k,O}) where {k,O}
    la,mda = Symbolics.latex(R.a)
    lb,mdb = Symbolics.latex(R.b)
    lc,mdc = Symbolics.latex(R.c)
    ld,mdd = Symbolics.latex(R.d)
    args,md = Symbolics.delimit("$(la),$(lb);$(lc),$(ld)", max(mda,mdb,mdc,mdd))
    "R^{$(k)}$(args)",md
end

function Symbolics.latex(F::DirectSlaterIntegral{k,O}) where {k,O}
    la,mda = Symbolics.latex(F.a)
    lb,mdb = Symbolics.latex(F.b)
    args,md = Symbolics.delimit("$(la)"*(!isdiagonal(F) ? ", $(lb)" : ""),max(mda,mdb))
    "F^{$(k)}$(args)",md
end

function Symbolics.latex(G::ExchangeSlaterIntegral{k, O}) where {k,O}
    la,mda = Symbolics.latex(G.a)
    lb,mdb = Symbolics.latex(G.b)
    args,md = Symbolics.delimit("$(la),$(lb)", max(mda,mdb))
    "G^{$(k)}$args",md
end

function Symbolics.latex(Y::SlaterPotential{k,O}) where {k,O}
    la,mda = Symbolics.latex(Y.a)
    lb,mdb = Symbolics.latex(Y.b)
    args,md = Symbolics.delimit("$(la)"*(!isdiagonal(Y) ? ", $(lb)" : ""),max(mda,mdb))
    "Y^{$(k)}$args",md
end
# * Symbolics registration

@new_number OverlapIntegral
@new_number DiagonalIntegral
@new_number GeneralSlaterIntegral
@new_number DirectSlaterIntegral
@new_number ExchangeSlaterIntegral
@new_number SlaterPotential

export OverlapIntegral, DiagonalIntegral, GeneralSlaterIntegral, DirectSlaterIntegral, ExchangeSlaterIntegral, SlaterPotential
