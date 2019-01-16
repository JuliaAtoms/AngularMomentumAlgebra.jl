using Symbolics

Base.:(==)(x::A, y::B) where {A<:Symbolic,B<:Symbolic} = false
Base.promote(x::A, y::B) where {A<:Symbolic,B<:Symbolic} =
    (SymExpr(:identity, [x]), SymExpr(:identity, [y]))

macro new_number(T)
    quote
        Base.:(==)(x::$(esc(T)), y::N) where {N<:Union{Real,Complex}} = false
        Base.:(==)(x::N, y::$(esc(T))) where {N<:Union{Real,Complex}} = false
        Base.:(==)(x::$(esc(T)), y::Sym) = false
        Base.:(==)(x::Sym, y::$(esc(T))) = false
        Base.:(==)(x::$(esc(T)), y::SymExpr) = false
        Base.:(==)(x::SymExpr, y::$(esc(T))) = false

        Base.promote(::Type{<:$(esc(T))}) = SymExpr

        Base.promote(x::TT, y::SymExpr) where {TT<:$(esc(T))} =
            (SymExpr(:identity, [x]), y)
        Base.promote(x::SymExpr, y::TT) where {TT<:$(esc(T))} =
            (x, SymExpr(:identity, [y]))
    end
end

Base.diff(::Union{Real,Complex}, orb::O, occ::I) where {O,I} = 0

function Base.diff(s::SymExpr, orb::O, occ::I) where {O,I}
    va = [diff(a, orb, occ) for a in s.args]

    if s.op == :(+)
        sum(va)
    elseif s.op == :(-)
        error("Not implemented")
    elseif s.op == :(*)
        sum([prod(vcat(s.args[1:i-1], va[i], s.args[i+1:end]))
             for i in eachindex(s.args)])
    end
end
