using Symbolics

macro new_number(T)
    quote
        Base.:(==)(x::$T, y) = false
        Base.:(==)(x, y::$T) = false
        Base.:(==)(x::$T, y::Number) = false
        Base.:(==)(x::Number, y::$T) = false
        Base.:(==)(x::$T, y::Sym) = false
        Base.:(==)(x::Sym, y::$T) = false
        Base.:(==)(x::$T, y::SymExpr) = false
        Base.:(==)(x::SymExpr, y::$T) = false
        
        Base.promote(::Type{<:$T}) = SymExpr
        
        Base.promote(x::TT, y::SymExpr) where {TT<:$T} =
            (SymExpr(:identity, [x]), y)
        Base.promote(x::SymExpr, y::TT) where {TT<:$T} =
            (x, SymExpr(:identity, [y]))
        
        # Base.promote(x::A, y::B) where {A<:Symbolic,B<:Symbolic} =
        #     (SymExpr(:identity, [x]), SymExpr(:identity, [y]))
        
        # Base.promote(x::TT, y::S) where {TT<:$T,S<:Symbolic} =
        #     (SymExpr(:identity, [x]), SymExpr(:identity, [y]))
        # Base.promote(x::S, y::TT) where {TT<:$T,S<:Symbolic} =
        #     (SymExpr(:identity, [x]), SymExpr(:identity, [y]))
    end
end

macro gen_compare_false(types...)
    for (i,A) in enumerate(types)
        for (j,B) in enumerate(types)
            j ≤ i && continue
            
            @eval Base.promote(x::$A, y::$B) =
                (SymExpr(:identity, [x]), SymExpr(:identity, [y]))
            
            @eval Base.:(==)(x::$A, y::$B) = false
            @eval Base.:(==)(x::$B, y::$A) = false
        end
    end
end

# macro new_numbers(types...)
#     for T in types
#         eval(@new_number(T))
#     end
#     @gen_compare_false(types...)
# end
