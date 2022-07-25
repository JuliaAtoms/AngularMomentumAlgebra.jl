import EnergyExpressions: showcoeff

"""
    LinearCombination(Ts::Vector{T}, coeffs::Vector{N})

Represents a general linear combination of objects of type `T`.
"""
struct LinearCombination{T,N<:Number}
    Ts::Vector{T}
    coeffs::Vector{N}
    function LinearCombination(Ts::Vector{T}, coeffs::Vector{N}) where {T,N}
        sel = .!iszero.(coeffs)
        new{T,N}(Ts[sel], coeffs[sel])
    end
end

Base.length(lc::LinearCombination) = length(lc.Ts)
Base.eltype(lc::LinearCombination{T,N}) where {T,N} = (T,N)

Base.:(==)(a::LinearCombination, b::LinearCombination) =
    a.Ts == b.Ts && a.coeffs == b.coeffs

Base.hash(lc::LinearCombination, h::UInt) =
    hash(lc.Ts, hash(lc.coeffs, h))

Base.zero(::LinearCombination{T,N}) where {T,N} =
    LinearCombination(Vector{T}(), Vector{N}())

Base.iszero(lc::LinearCombination) = isempty(lc.Ts)

Base.getindex(lc::LinearCombination, i) =
    (lc.Ts[i],lc.coeffs[i])

Base.iterate(lc::LinearCombination, i=1) =
    i > length(lc) ? nothing : (lc[i], i+=1)

function Base.show(io::IO, lc::LinearCombination)
    n = length(lc)
    for (i,(t,c)) in enumerate(lc)
        showcoeff(io, c, i > 1)
        show(io, t)
        i < n && write(io, " ")
    end
end

Base.:(+)(A::LinearCombination, B::LinearCombination) =
    LinearCombination(vcat(A.Ts, B.Ts), vcat(A.coeffs, B.coeffs))

Base.:(-)(A::LinearCombination, B::LinearCombination) =
    LinearCombination(vcat(A.Ts, B.Ts), vcat(A.coeffs, -1 * B.coeffs))

Base.:(+)(A::LC, B::T) where {T,N,LC<:LinearCombination{<:T,N}} =
    LinearCombination(vcat(A.Ts, B), vcat(A.coeffs, one(N)))

Base.:(+)(A::T, B::LC) where {T,N,LC<:LinearCombination{<:T,N}} =
    LinearCombination(vcat(A, B.Ts), vcat(one(N), B.coeffs))

Base.:(-)(A::LC, B::T) where {T,N,LC<:LinearCombination{<:T,N}} =
    A + (-1)*B

Base.:(-)(A::T, B::LC) where {T,N,LC<:LinearCombination{<:T,N}} =
    A + (-1)*B

Base.:(-)(A::LC) where {LC<:LinearCombination} =
    LinearCombination(A.Ts, -1 * A.coeffs)

Base.:(*)(A::LinearCombination, b::Number) =
    LinearCombination(A.Ts, b * A.coeffs)

Base.:(*)(a::Number, B::LinearCombination) =
    LinearCombination(B.Ts, a * B.coeffs)

Base.:(/)(A::LinearCombination, b::Number) =
    LinearCombination(A.Ts, A.coeffs/b)

"""
    @linearly_combinable TT

Turns the type `TT` into a linearly combinable type, i.e. defines
arithmetic operators.

# Examples

```julia
julia> @linearly_combinable Symbol

julia> 4*(:x) - 5*(:y)
4 :x - 5 :y
```
"""
macro linearly_combinable(TT)
    quote
        Base.:(+)(A::T, B::U) where {T<:$(esc(TT)),U<:$(esc(TT))} =
            LinearCombination(promote_type(T,U)[A,B], [1,1])

        Base.:(-)(A::$(esc(TT))) =
            LinearCombination($(esc(TT))[A], [-1])

        Base.:(*)(a::Number, B::$(esc(TT))) =
            LinearCombination([B], [a])

        Base.:(*)(A::$(esc(TT)), b::Number) =
            LinearCombination([A], [b])
    end
end
