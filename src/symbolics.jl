using Symbolics

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

Symbolics.latex(o::Orbital) =
    "$(o.n)\\mathrm{$(spectroscopic_label(o.ℓ))}",0
