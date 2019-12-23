"""
    wrap_show_debug(expr)

Wrap `expr` with a `@show` macro call if `DEBUG_TENSOR_DSL` exists as
an environment flag, otherwise return `expr` as-is. This is applied at
compile time, so no overhead is incurred.
"""
function wrap_show_debug(expr)
    if "DEBUG_TENSOR_DSL" in keys(ENV)
        s = :(@show)
        push!(s.args, expr)
        s
    else
        expr
    end
end

remove_line_numbers(exprs) = filter(n -> !(n isa LineNumberNode), exprs)

generate_signature(f::Function, sym, TensorType::Symbol) =
    Expr(:call, sym, f(Expr(:(::), :tensor, TensorType))...)

generate_signature(f::Function, sym, TensorType::Expr) =
    Expr(:where,
         Expr(:call, sym, f(Expr(:(::), :tensor, TensorType.args[1]))...),
         TensorType.args[2])

"""
    strip′(s::Symbol)

Strip `s` of a trailing `′`, error if it is not present.
"""
function strip′(s::Symbol)
    s = string(s)
    s[end] == '′' || throw(ArgumentError("$s does not end in `′`"))
    Symbol(s[1:prevind(s, end, 1)])
end

"""
    identify_quantum_numbers(selection_rules)

Given the block of selection rules, identify the quantum numbers
pertaining to the tensor under consideration. It is assumed that the
left-hand side consists of a primed quantum number only, which is
specified to be equal to an expression or belonging to a set/interval.
"""
function identify_quantum_numbers(selection_rules::Expr)
    if selection_rules.head == :call
        γj′ = selection_rules.args[2]
        γj′, strip′(γj′)
    elseif selection_rules.head == :block
        γj′ = Symbol[]
        for selection_rule in remove_line_numbers(selection_rules.args)
            push!(γj′, selection_rule.args[2])
        end
        Expr(:tuple, γj′...),Expr(:tuple, map(strip′,γj′)...)
    end
end

recursepm(s) = [s]

"""
    recursepm(e::Expr)

Recursively expand the expression `e` for instances of `±` and `∓`,
generating `+` and `-` minus branches.
"""
function recursepm(e::Expr)
    h = e.head
    if first(e.args) == :(±)
        a = map(recursepm, e.args[2:end])
        ee = map([:(+),:(-)]) do op
            map(Iterators.product(a...)) do a
                Expr(h, op, a...)
            end
        end
        collect(Iterators.flatten(ee))
    elseif first(e.args) == :(∓)
        a = map(recursepm, e.args[2:end])
        ee = map([:(-),:(+)]) do op
            map(Iterators.product(a...)) do a
                Expr(h, op, a...)
            end
        end
        collect(Iterators.flatten(ee))
    else
        a = map(recursepm, e.args)
        ee = map(Iterators.product(a...)) do a
            Expr(h, a...)
        end
        vec(ee)
    end
end

"""
    recursepm(lhs, rhs)

Recursively expand `lhs` and `rhs` for instances of `±` and `∓` and
join all resulting cases into a short-circuiting `||` test.
"""
function recursepm(lhs, rhs)
    ls = recursepm(lhs)
    rs = recursepm(rhs)
    cases = vec(map(((l,r),) -> :($l == $r), Iterators.product(ls, rs)))
    if length(cases) == 1
        first(cases)
    else
        a = pop!(cases)
        b = pop!(cases)
        e = Expr(:(||), a, b)
        while !isempty(cases)
            a = pop!(cases)
            e = Expr(:(||), a, e)
        end
        e
    end
end

"""
    parse_selection_rule(rule)

Parse a selection rule in the DSL of [`@tensor`](@ref); if the
selection rule is an equality, its left- and right-hand sides are
recursively expanded for possible cases of `±` and `∓` with
[`recursepm`](@ref).
"""
function parse_selection_rule(rule)
    rule = if first(rule.args) == :(==)
        recursepm(rule.args[2], rule.args[3])
    else
        rule
    end
    :($rule || return true)
end

function push_selection_rules!(f, args, selection_rules)
    for selection_rule in (selection_rules.head == :call ?
                           [selection_rules] :
                           selection_rules.args)
        if selection_rule isa Expr &&
            selection_rule.head == :call &&
            # ~ implies this particular quantum number does not
            # influence whether this matrix element vanishes or not;
            # this is typically the principal quantum number n.
            first(selection_rule.args) ≠ :(~)
            push!(args, f(selection_rule))
        elseif selection_rule isa LineNumberNode
            push!(args, selection_rule)
        end
    end
end

"""
    generate_iszero(TensorType, selection_rules)

Generate a function for testing if the matrix element of `TensorType`
vanishes, given a set a quantum number deduced from `selection_rules`.
"""
function generate_iszero(TensorType, selection_rules)
    γj′,γj = identify_quantum_numbers(selection_rules)

    signature = generate_signature(t -> (γj′, t, γj), :(Base.iszero), TensorType)
    body = Expr(:block)
    push_selection_rules!(parse_selection_rule, body.args, selection_rules)
    push!(body.args, false)
    # TODO: Should auto-generate docstring for generated iszero
    # function.
    Expr(:function, signature, body)
end

"""
    generate_rme(TensorType, selection_rules, doc, rme)

Generate a function for computing the reduced matrix element of
`TensorType`, given a set a quantum number deduced from
`selection_rules`, along with the docstring and the definition `rme`.
"""
function generate_rme(TensorType, selection_rules, doc, rme)
    γj′,γj = identify_quantum_numbers(selection_rules)
    signature = generate_signature(t -> (γj′, t, γj), :rme, TensorType)
    body = Expr(:block, :(iszero($(γj′), tensor, $(γj)) && return 0), wrap_show_debug(rme))
    fun = Expr(:function, signature, body)
    :(@doc($doc, $fun))
end

"""
    generate_copulings(TensorType, selection_rules)

Generate a function that given the quantum numbers `γj` generates
lists of all permissible `γj′` for which the reduced matrix element
`⟨γj′||::TensorType||γj⟩` does not vanish. This is deduced from the
`selection_rules`.
"""
function generate_couplings(TensorType, selection_rules)
    indeterminates = filter(r -> r.head == :call && first(r.args) == :(~),
                            remove_line_numbers(selection_rules.head == :call ?
                                                [selection_rules] :
                                                selection_rules.args))
    γj′,γj = identify_quantum_numbers(selection_rules)
    signature = generate_signature(t -> (t, γj), :couplings, TensorType)
    body = Expr(:block)
    if !isempty(indeterminates)
        push!(body.args, Expr(:call, :error, "Cannot generate all couplings for "*
                              string(TensorType)*
                              " due to the following indeterminate quantum numbers: "*
                              join(string.(indeterminates), ", ")))
    else
        push_selection_rules!(body.args, selection_rules) do rule
            op = rule.args[1]
            lhs, rhs = rule.args[2], rule.args[3]
            rhs = if op == :(∈)
                rhs
            elseif op == :(==)
                if rhs isa Symbol
                    Expr(:vect, rhs)
                else
                    throw(ArgumentError("Don't know how to handle $(rhs)"))
                end
            else
                throw(ArgumentError("Unknown operator $op"))
            end
            Expr(:(=), lhs, rhs)
        end
        push!(body.args, γj′)
    end
    γj′s = γj′ isa Symbol ? string(γj′) : join(string.(γj′.args), "")
    γjs = γj isa Symbol ? string(γj) : join(string.(γj.args), "")
    TS = string(TensorType isa Expr ? TensorType.args[1] : TensorType)
    doc = """
    $(signature)

Generate all quantum numbers `$(γj′s)` for which
`⟨$(γj′s)||::$(TS)||$(γjs)⟩` does not vanish.
"""
    fun = Expr(:function, signature, body)
    :(@doc($doc, $fun))
end

"""
    @tensor(exprs, TensorType)

Macro to generate `Base.iszero`, [`rme`](@ref) and [`couplings`](@ref)
for `TensorType`, given a set of selection rules and an expression for
the reduced matrix element.
"""
macro tensor(exprs, TensorType)
    args = remove_line_numbers(exprs.args[2].args)
    selection_rules = args[1]
    iszero = generate_iszero(TensorType, selection_rules)
    rme = generate_rme(TensorType, selection_rules, args[2], args[3])
    couplings = generate_couplings(TensorType, selection_rules)
    esc(Expr(:block, iszero, rme, couplings))
end
