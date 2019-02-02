"""
    coupled_states(hamiltonians...[; i₀=1])

Find all states coupled by the `hamiltonians`, starting from the state
with index `i₀`. This can be useful to reduce the necessary basis or
to generate invariant sets for split-operator propagation.
"""
function coupled_states(hamiltonians::Matrix{<:Union{<:OneBodyEnergyExpression,
                                                     <:TwoBodyEnergyExpression}}...;
                        i₀=1)
    couplings = SparseMatrixCSC(sum(.!iszero.(h) for h in hamiltonians))

    m = size(hamiltonians[1],1)
    visited = falses(m)
    visited[i₀] = true

    rows = rowvals(couplings)

    istack = [i₀]
    while !isempty(istack)
        icur = pop!(istack)
        neighbours = rows[nzrange(couplings, icur)]
        for n in neighbours
            if !visited[n]
                push!(istack, n)
                visited[n] = true
            end
        end
    end

    visited
end

function invariant_sets(hamiltonians::Matrix{<:Union{<:OneBodyEnergyExpression,
                                                     <:TwoBodyEnergyExpression}}...)
    m = size(hamiltonians[1],1)
    visited = falses(m)

    sets = Vector{Vector{Bool}}()

    while !all(visited)
        icur = findfirst(.!visited)
        set = coupled_states(hamiltonians...; i₀=icur)
        push!(sets, set)
        visited[:] .|= set
    end
    sets
end

export coupled_states, invariant_sets
