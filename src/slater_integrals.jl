abstract type AbstractRadialIntegral{O<:Orbital} <: Symbolic end
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

struct DiagonalIntegral{O} <: AbstractRadialIntegral{O}
    o::O
end
Base.:(==)(a::DiagonalIntegral, b::DiagonalIntegral) = a.o == b.o

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

struct DirectSlaterIntegral{k,O} <: AbstractSlaterIntegral{k,O}
    a::O
    b::O
end
DirectSlaterIntegral(k::Integer, a::O, b::O) where {O} =
    DirectSlaterIntegral{k,O}(a, b)
Base.:(==)(a::DirectSlaterIntegral{k}, b::DirectSlaterIntegral{k′}) where {k,k′} =
    k == k′ && a.a == b.a && a.b == b.b

struct ExchangeSlaterIntegral{k,O} <: AbstractSlaterIntegral{k,O}
    a::O
    b::O
end
ExchangeSlaterIntegral(k::Integer, a::O, b::O) where {O} =
    ExchangeSlaterIntegral{k,O}(a, b)
Base.:(==)(a::ExchangeSlaterIntegral{k}, b::ExchangeSlaterIntegral{k′}) where {k,k′} =
    k == k′ && a.a == b.a && a.b == b.b

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
    write(io, ", ")
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
    write(io, ", ")
    show(io, F.b)
    write(io, ")")
end

function Base.show(io::IO, G::ExchangeSlaterIntegral{k, O}) where {k,O}
    write(io, "G", to_superscript(k), "(")
    show(io, G.a)
    write(io, ", ")
    show(io, G.b)
    write(io, ")")
end

@new_number OverlapIntegral
@new_number DiagonalIntegral
@new_number GeneralSlaterIntegral
@new_number DirectSlaterIntegral
@new_number ExchangeSlaterIntegral
@gen_compare_false OverlapIntegral DiagonalIntegral GeneralSlaterIntegral DirectSlaterIntegral ExchangeSlaterIntegral

export OverlapIntegral, DiagonalIntegral, GeneralSlaterIntegral, DirectSlaterIntegral, ExchangeSlaterIntegral
