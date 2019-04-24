import EnergyExpressions: showcoeff

"""
    LinearCombination(Ts::Vector{T}, coeffs::Vector{N})

Represents a general linear combination of objects of type `T`.
"""
struct LinearCombination{T,N<:Number}
    Ts::Vector{T}
    coeffs::Vector{N}
end

Base.length(lc::LinearCombination) = length(lc.Ts)
Base.eltype(lc::LinearCombination{T,N}) where {T,N} = (T,N)

Base.iterate(lc::LinearCombination, i=1) =
    i > length(lc) ? nothing : ((lc.Ts[i],lc.coeffs[i]), i+=1)

function Base.show(io::IO, lc::LinearCombination)
    n = length(lc)
    for (i,(t,c)) in enumerate(lc)
        showcoeff(io, c, i > 1)
        write(io, " ")
        show(io, t)
        i < n && write(io, " ")
    end
end

Base.:(+)(A::LC, B::LC) where {LC<:LinearCombination} =
    LinearCombination(vcat(A.Ts, B.Ts), vcat(A.coeffs, B.coeffs))

Base.:(-)(A::LC, B::LC) where {LC<:LinearCombination} =
    LinearCombination(vcat(A.Ts, B.Ts), vcat(A.coeffs, -1 * B.coeffs))

Base.:(+)(A::LC, B::T) where {T,N,LC<:LinearCombination{T,N}} =
    LinearCombination(vcat(A.Ts, B), vcat(A.coeffs, one(N)))

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
        Base.:(+)(A::T, B::T) where {T<:$(esc(TT))} =
            LinearCombination(T[A,B], [1,1])

        Base.:(-)(A::$(esc(TT))) =
            LinearCombination($(esc(TT))[A], [-1])

        Base.:(*)(a::Number, B::$(esc(TT))) =
            LinearCombination([B], [a])

        Base.:(*)(A::$(esc(TT)), b::Number) =
            LinearCombination([A], [b])
    end
end
