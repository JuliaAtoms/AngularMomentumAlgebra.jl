struct RepulsionPotential{kind,A,B,O}
    a::A
    b::B
    o::O # The orbital o feels the potentials formed by a & b
end
RepulsionPotential{kind}(a::A, b::B, o::O) where {kind,A,B,O} =
    RepulsionPotential{kind,A,B,O}(a,b,o)

Base.zero(::Type{RepulsionPotential{kind}}) where kind =
    RepulsionPotential{kind}(0,0,0)

Base.iszero(JK::RepulsionPotential) =
    JK.a == JK.b == JK.o == 0

Base.:(==)(a::RepulsionPotential{kind}, b::RepulsionPotential{kind′}) where {kind,kind′} =
    kind == kind′ && a.a == b.a && a.b == b.b && a.o == b.o ||
    a.a == a.b == b.a == b.b && a.o == b.o

const DirectPotential{A,B,O} = RepulsionPotential{:direct,A,B,O}
const ExchangePotential{A,B,O} = RepulsionPotential{:exchange,A,B,O}

function Base.show(io::IO, RP::RepulsionPotential{:general})
    write(io, "[")
    show(io, RP.a)
    write(io, "|")
    show(io, RP.b)
    write(io, "]")
    show(io, RP.o)
end

function Base.show(io::IO, J::DirectPotential)
    write(io, "Ĵ{")
    show(io, J.a)
    write(io, ";")
    show(io, J.b)
    write(io, "}")
    show(io, J.o)
end

function Base.show(io::IO, K::ExchangePotential)
    write(io, "K̂{")
    show(io, K.a)
    write(io, ";")
    show(io, K.b)
    write(io, "}")
    show(io, K.o)
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
    write(io, "[")
    show(io, JK.a)
    write(io, "||")
    show(io, JK.b)
    write(io, "]")
    show(io, JK.o)
end

export RepulsionPotential, DirectPotential, ExchangePotential, DirectExchangePotentials
