abstract type AbstractRadialIntegral{T,O<:Orbital} end
abstract type AbstractSlaterIntegral{k,T,O} <: AbstractRadialIntegral{T,O} end

struct DiagonalIntegral{T,O} <: AbstractRadialIntegral{T,O}
    w::T
    o::O
end
DiagonalIntegral(w::T,o::O) where {T,O} =
    DiagonalIntegral{T,O}(w,o)

struct GeneralSlaterIntegral{k,T,O} <: AbstractSlaterIntegral{k,T,O}
    w::T
    a::O
    b::O
    c::O
    d::O
end
GeneralSlaterIntegral(k::Integer, w::T, a::O, b::O, c::O, d::O) where {T,O} =
    GeneralSlaterIntegral{k,T,O}(w, a, b, c, d)

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

function Base.show(io::IO, I::DiagonalIntegral{T,O}) where {T,O}
    show(io, I.w)
    write(io, " × I(")
    show(io, I.o)
    write(io, ", ")
    show(io, I.o)
    write(io, ")")
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

export DiagonalIntegral, GeneralSlaterIntegral, DirectSlaterIntegral, ExchangeSlaterIntegral
