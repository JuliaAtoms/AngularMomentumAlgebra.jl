# * Tensors

abstract type AbstractTensor end

struct BasisTensor{T<:AbstractTensor,q} <: AbstractTensor end

function Base.show(io::IO, ::Type{BasisTensor{T,q}}) where {T,q}
    show(io, T)
    write(io,to_subscript(q))
end

struct TensorContraction{A<:AbstractTensor,B<:AbstractTensor} <: AbstractTensor end
struct TensorProduct{A,B} <: AbstractTensor
    a::A
    b::B
end

LinearAlgebra.dot(::Type{TA}, ::Type{TB}) where {TA<:AbstractTensor,TB<:AbstractTensor} =
    TensorContraction{TA,TB}

function Base.show(io::IO, ::Type{TensorContraction{TA,TB}}) where {TA,TB}
    show(io, TA)
    write(io, "⋅")
    show(io, TB)
end

# ** Spherical tensors

struct SphericalTensor{k} <: AbstractTensor end
order(st::SphericalTensor{k}) where k = k

function Base.show(io::IO, ::Type{SphericalTensor{k}}) where k
    write(io,to_boldface("C"))
    write(io,"⁽",to_superscript(k),"⁾")
end

# * Basis vectors

struct Bra{B} <: AbstractTensor
    v::B
end

struct Ket{K} <: AbstractTensor
    v::K
end

struct Braket{B,K} <: AbstractTensor
    b::B
    k::K
end
Braket(bra::Bra{B}, ket::Ket{K}) where {B,K} =
    Braket{B,K}(bra.v, ket.v)

Base.adjoint(bra::Bra) = Ket(bra.v)
Base.adjoint(ket::Ket) = Bra(ket.v)

Base.:(*)(bra::Bra, ket::Ket) = Braket(bra, ket)

function Base.show(io::IO, bra::Bra)
    write(io, "⟨")
    show(io, bra.v)
    write(io, "|")
end

function Base.show(io::IO, ket::Ket)
    write(io, "|")
    show(io, ket.v)
    write(io, "⟩")
end

function Base.show(io::IO, braket::Braket)
    write(io, "⟨")
    show(io, braket.b)
    write(io, "|")
    show(io, braket.k)
    write(io, "⟩")
end

# ** Angular momenta
struct AngularMomentum
    j::HalfInteger
end

Base.:(*)(bra::Bra{AngularMomentum}, ket::Ket{AngularMomentum}) =
    Kronecker(bra.v, ket.v)

function Base.show(io::IO, a::AngularMomentum)
    r = convert(Rational, a.j)
    write(io, "$(numerator(r))")
    denominator(r) == 2 && write(io, "/2")
end

# * Reduced matrix elements

struct ReducedMatrixElement{T<:AbstractTensor,B,K}
    bra::Bra{B}
    ket::Ket{K}
end

ReducedMatrixElement(::Type{T},bra::Bra{B},ket::Ket{K}) where {T,B,K} =
        ReducedMatrixElement{T,B,K}(bra, ket)

function Base.show(io::IO, rme::ReducedMatrixElement{T,B,K}) where {T,B,K}
    show(io, rme.bra)
    write(io, "|")
    show(io, T)
    write(io, "|")
    show(io, rme.ket)
end

function Base.show(io::IO, rme::ReducedMatrixElement{SphericalTensor{k},A,A}) where {k,A<:AngularMomentum}
    show(io, rme.bra)
    write(io, "|")
    show(io, SphericalTensor{k})
    write(io, "|")
    show(io, rme.ket)

    ℓ = convert(Rational, rme.bra.v.j)
    denominator(ℓ) == 1 && (ℓ = ℓ |> Int)
    ℓ′ = convert(Rational, rme.ket.v.j)
    denominator(ℓ′) == 1 && (ℓ′ = ℓ′ |> Int)

    p = (-1)^ℓ
    write(io, " = $(p > 0 ? "" : "-")")
    write(io, to_root((2ℓ+1)*(2ℓ′+1) |> Int))
    show(io, IIIJ(ℓ, k, ℓ′, 0, 0, 0))
end

Base.:(|)(bra::Bra{B}, ::Type{T}) where {B,T<:AbstractTensor} =
    TensorProduct(bra, T)

Base.:(|)(tp::TensorProduct{A,T}, k::Ket{K}) where {A,T,K} =
    ReducedMatrixElement(tp.b, tp.a, k)

# ** Spherical tensors
function Base.convert(::Type{T}, rme::ReducedMatrixElement{SphericalTensor{k},AngularMomentum,AngularMomentum}) where {T<:Number,k}
    @unpack bra,ket = rme
    j = convert(Rational,bra.v.j)
    j′ = convert(Rational,ket.v.j)
    (-one(T))^j*√((2j+1)*(2j′+1))*wigner3j(j,k,j′,0,0,0)
end

# * Exports

export SphericalTensor, order, dot, ⋅, Bra, Ket, Braket, AngularMomentum
