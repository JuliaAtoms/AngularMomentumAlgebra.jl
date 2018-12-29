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
            j < i && continue

            @eval Base.promote(x::$A, y::$B) =
                (SymExpr(:identity, [x]), SymExpr(:identity, [y]))

            j == i && continue

            @eval Base.promote(x::$B, y::$A) = promote(y, x)

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

latex(s::Number) = "$(s)"
function latex(z::Complex{T}) where T
    a,b = real(z),imag(z)
    if a == zero(T)
        latex(b) * "\\mathrm{i}"
    elseif b == zero(T)
        latex(a)
    else
        latex(a) * (b < 0 ? "" : "+") * latex(b) * "\\mathrm{i}"
    end
end
latex(s::Sym) = "$(s)"
latex_op(op::Sym) = op.name == :(*) ? "" : latex(op)
latex_op(op) = latex(op)

function latex(expr::SymExpr)
    op = expr.op
    args = map(expr.args) do arg
        if op.name != :(+) && arg isa Number && !(arg isa SymExpr) &&
            arg isa Real && arg < 0 || arg isa Complex
            "{($(latex(arg)))}"
        else
            "{$(latex(arg))}"
        end
    end
    lop = latex_op(op)
    if op.name == :(+)
        S = args[1]
        for a in args[2:end]
            S *= a[2] == '-' ? "" : "+"
            S *= a
        end
        S
    elseif op.name == :(^)
        arg1 = expr.args[1]
        if arg1 isa Number && !(arg1 isa SymExpr)
            if arg1 isa Real && arg1 == -1
                "(-)^{$(args[2])}"
            else
                "$(args[1])^{$(args[2])}"
            end
        else
            "\\left(args[1]\\right)^{$(args[2])}"
        end
    elseif op.name == :(*)
        if expr.args[1] == -1
            "-" * join(args[2:end],"")
        else
            join(args, "")
        end
    else
        join(args, lop)
    end
end
