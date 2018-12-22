abstract type AbstractRadialIntegral{O<:Orbital} end
abstract type AbstractSlaterIntegral{k,T,O} <: AbstractRadialIntegral{O} end

struct OverlapIntegral{O} <: AbstractRadialIntegral{O}
    p::Int
    a::O
    b::O
end
OverlapIntegral(p::Integer, a::O, b::O) where {O} =
    OverlapIntegral{O}(p, a, b)

struct DiagonalIntegral{T,O,N} <: AbstractRadialIntegral{O}
    w::T
    o::O
    overlaps::NTuple{N,OverlapIntegral{O}}
end
DiagonalIntegral(w::T, o::O, overlaps::OverlapIntegral{O}...) where {T,O} =
    DiagonalIntegral{T,O,length(overlaps)}(w, o, overlaps)

struct GeneralSlaterIntegral{k,T,O,N} <: AbstractSlaterIntegral{k,T,O}
    w::T
    a::O
    b::O
    c::O
    d::O
    overlaps::NTuple{N,OverlapIntegral{O}}
end
GeneralSlaterIntegral(k::Integer, w::T, a::O, b::O, c::O, d::O, overlaps::OverlapIntegral{O}...) where {T,O} =
    GeneralSlaterIntegral{k,T,O,length(overlaps)}(w, a, b, c, d, overlaps)

struct DirectSlaterIntegral{k,T,O} <: AbstractSlaterIntegral{k,T,O}
    w::T
    a::O
    b::O
end
DirectSlaterIntegral(w::T, a::O, b::O) where {T,O} =
    DirectSlaterIntegral{T,O}(w, a, b)

struct ExchangeSlaterIntegral{k,T,O} <: AbstractSlaterIntegral{k,T,O}
    w::T
    a::O
    b::O
end
ExchangeSlaterIntegral(w::T, a::O, b::O) where {T,O} =
    ExchangeSlaterIntegral{T,O}(w, a, b)



function Base.show(io::IO, o::OverlapIntegral{O}) where O
    write(io, "⟨")
    show(io, o.a)
    write(io, "|")
    show(io, o.b)
    write(io, "⟩")
    o.p != 1 && write(io, to_superscript(o.p))
end

function Base.show(io::IO, I::DiagonalIntegral{T,O}) where {T,O}
    show(io, I.w)
    write(io, " × I(")
    show(io, I.o)
    write(io, ", ")
    show(io, I.o)
    write(io, ")")
    for ov in I.overlaps
        write(io, " × ")
        show(io, ov)
    end
end

function Base.show(io::IO, R::GeneralSlaterIntegral{k,T,O}) where {k,T,O}
    show(io, R.w)
    write(io, " × R", to_superscript(k), "(")
    show(io, R.a)
    write(io, ", ")
    show(io, R.b)
    write(io, "; ")
    show(io, R.c)
    write(io, ", ")
    show(io, R.d)
    write(io, ")")
    for ov in R.overlaps
        write(io, " × ")
        show(io, ov)
    end
end

function Base.show(io::IO, F::DirectSlaterIntegral{k,T,O}) where {k,T,O}
    show(io, F.w)
    write(io, " × F", to_superscript(k), "(")
    show(io, F.a)
    write(io, ", ")
    show(io, F.b)
    write(io, ")")
end

function Base.show(io::IO, G::ExchangeSlaterIntegral{k,T,O}) where {k,T,O}
    show(io, G.w)
    write(io, " × G", to_superscript(k), "(")
    show(io, G.a)
    write(io, ", ")
    show(io, G.b)
    write(io, ")")
end

export OverlapIntegral, DiagonalIntegral, GeneralSlaterIntegral, DirectSlaterIntegral, ExchangeSlaterIntegral
