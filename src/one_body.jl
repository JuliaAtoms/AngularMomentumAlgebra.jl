# * OneBodyHamiltonian operator

struct OneBodyHamiltonian{O}
    orb::O
end
Base.:(==)(a::OneBodyHamiltonian, b::OneBodyHamiltonian) = a.orb == b.orb
function Base.show(io::IO, h::OneBodyHamiltonian)
    h.orb isa Conjugate && show(io, h.orb)
    write(io, "ĥ")
    !(h.orb isa Conjugate) || show(io, h.orb)
end

Base.iszero(h::OneBodyHamiltonian) = h.orb == 0

# * One-body integral

struct OneBodyIntegral{A,B}
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
    OneBodyHamiltonian(I.b == orb ? conj(I.a) : 0)

Base.diff(I::OneBodyIntegral{A,B}, corb::Conjugate{O}) where {A,B,O,II} =
    OneBodyHamiltonian(I.a == corb.orbital ? I.b : 0)

export OneBodyHamiltonian, OneBodyIntegral
