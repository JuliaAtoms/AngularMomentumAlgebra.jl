struct Kronecker{T}
    a::T
    b::T
end
Base.convert(::Type{N}, k::Kronecker{T}) where {N,T} =
    k.a == k.b ? one(N) : zero(N)

Base.show(io::IO, k::Kronecker) =
    write(io, "δ($(k.a),$(k.b))")
