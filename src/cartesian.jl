"""
    cartesian_tensor_component(t::Tensor{1}, c)

This returns the Cartesian tensor component `c` (valid choices are
`:x`, `:y`, or `:z`) of the rank-1 tensor `t`, as a linear combination
of its "natural" [`TensorComponent`](@ref)s `1`, `0`, and `-1`. The
transform matrix is `M(+1, 0, -1 ← x, y, z)` given in Table 1.2 on
p. 14 of Varshalovich (1988).
"""
function cartesian_tensor_component(t::Tensor{1}, c::Symbol)
    if c == :x
        (-TensorComponent(t, 1) + TensorComponent(t, -1))/√2
    elseif c == :y
        (im*TensorComponent(t, 1) + im*TensorComponent(t, -1))/√2
    elseif c == :z
        TensorComponent(t, 0)
    else
        throw(ArgumentError("Unknown Cartesian tensor component $(c)"))
    end
end

export cartesian_tensor_component
