using Symbolics

macro new_number(T)
    quote
        Base.:(==)(x::$(esc(T)), y) = false
        Base.:(==)(x, y::$(esc(T))) = false
        Base.:(==)(x::$(esc(T)), y::Number) = false
        Base.:(==)(x::Number, y::$(esc(T))) = false
        Base.:(==)(x::$(esc(T)), y::Sym) = false
        Base.:(==)(x::Sym, y::$(esc(T))) = false
        Base.:(==)(x::$(esc(T)), y::SymExpr) = false
        Base.:(==)(x::SymExpr, y::$(esc(T))) = false

        Base.promote(::Type{<:$(esc(T))}) = SymExpr

        Base.promote(x::TT, y::SymExpr) where {TT<:$(esc(T))} =
            (SymExpr(:identity, [x]), y)
        Base.promote(x::SymExpr, y::TT) where {TT<:$(esc(T))} =
            (x, SymExpr(:identity, [y]))

        # Base.promote(x::A, y::B) where {A<:Symbolic,B<:Symbolic} =
        #     (SymExpr(:identity, [x]), SymExpr(:identity, [y]))

        # Base.promote(x::TT, y::S) where {TT<:$(esc(T)),S<:Symbolic} =
        #     (SymExpr(:identity, [x]), SymExpr(:identity, [y]))
        # Base.promote(x::S, y::TT) where {TT<:$(esc(T)),S<:Symbolic} =
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

# macro new_numbers(types...)
#     for T in types
#         eval(@new_number(T))
#     end
#     @gen_compare_false(types...)
# end

latex(s::Number) = "$(s)"
latex(s::AbstractFloat) = @sprintf("%0.10g", s)

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

latex(o::Orbital) =
    "$(o.n)\\mathrm{$(spectroscopic_label(o.ℓ))}"

function latex(expr::SymExpr)
    op = expr.op
    args = map(enumerate(expr.args)) do (i,arg)
        arg_negative = (arg isa Number && !(arg isa SymExpr) &&
                        arg isa Real && arg < 0)
        if op.name != :(+) && arg_negative && i > 1 || arg isa Complex
            "($(latex(arg)))"
        else
            "$(latex(arg))"
        end
    end
    lop = latex_op(op)
    if op.name == :(+)
        S = args[1]
        for a in args[2:end]
            S *= first(lstrip(a, ['{'])) == '-' ? "" : "+"
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
