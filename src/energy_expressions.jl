#=

Some references below to

- Cook, D. (2005). Handbook of computational quantum
  chemistry. Mineola, N.Y: Dover Publications, ISBN: 9780486443072.

=#

struct OneBodyEnergyExpression{A<:SpinOrbital,B<:SpinOrbital}
    integrals::Vector{OneBodyIntegral{A,B}}
    signs::Vector{Int}
    function OneBodyEnergyExpression{A,B}(a::Configuration{<:SpinOrbital},
                                          b::Configuration{<:SpinOrbital}) where {A,B}
        m = length(a)
        n = length(b)
        N = m*n
        integrals = Vector{OneBodyIntegral{A,B}}(undef, N)
        signs = Vector{Int}(undef, N)
        ii = 1
        for i = 1:m
            ai = a.orbitals[i]
            for j = 1:n
                # Cook 2005, p. 63
                integrals[ii] = OneBodyIntegral{A,B}(ai,b.orbitals[j])
                signs[ii] = (-1)^(i+j)
                ii += 1
            end
        end
        new{A,B}(integrals, signs)
    end
end

function OneBodyEnergyExpression(a::Configuration{<:SpinOrbital},
                                 b::Configuration{<:SpinOrbital})
    A = promote_type(typeof.(a.orbitals)...)
    B = promote_type(typeof.(b.orbitals)...)
    OneBodyEnergyExpression{A,B}(a, b)
end

function Base.:(==)(a::OneBodyEnergyExpression, b::OneBodyEnergyExpression)
    length(a.integrals) == length(b.integrals) || return false
    # This comparison assumes all integrals are unique.
    for (ai,as) in zip(a.integrals,a.signs)
        i = findfirst(isequal(ai), b.integrals)
        i === nothing && return false
        b.signs[i] == as || return false
    end
    true
end

function Base.adjoint(eng::OneBodyEnergyExpression{A,B}) where {A,B}
    integrals = map(eng.integrals) do integral
        OneBodyIntegral{B,A}(integral.b, integral.a)
    end
    OneBodyEnergyExpression{B,A}(integrals, eng.signs)
end

function Base.show(io::IO, eng::OneBodyEnergyExpression)
    if isempty(eng.integrals)
        write(io, "0")
        return
    end
    skip = if get(io, :compact, false)
        n = length(eng.integrals)
        if n > 2
            2:n-1
        else
            []
        end
    else
        []
    end
    for (i,(integral,sign)) in enumerate(zip(eng.integrals,eng.signs))
        if i ∈ skip
            i == skip[1] && write(io, " + …")
            continue
        end
        (i > 1 || sign != 1) && write(io, " ", sign == 1 ? "+" : "-", " ")
        write(io, string(integral))
    end
end

struct TwoBodyEnergyExpression{A<:SpinOrbital,B<:SpinOrbital}
    integrals::Vector{DirectExchangeIntegral{A,B}}

    function TwoBodyEnergyExpression{A,B}(a::Configuration{<:SpinOrbital},
                                          b::Configuration{<:SpinOrbital}) where {A,B}
        m = length(a)
        n = length(b)

        integrals = Vector{DirectExchangeIntegral{A,B}}()
        map(1:m) do i
            ai = a.orbitals[i]
            # Cook 2005, p. 69
            append!(integrals,
                    map(j -> DirectExchangeIntegral{A,B}(ai,b.orbitals[j]), 1:i-1))
        end |> oo -> vcat(oo...)
        new{A,B}(integrals)
    end
end

function TwoBodyEnergyExpression(a::Configuration{<:SpinOrbital},
                                 b::Configuration{<:SpinOrbital})
    A = promote_type(typeof.(a.orbitals)...)
    B = promote_type(typeof.(b.orbitals)...)
    TwoBodyEnergyExpression{A,B}(a,b)
end

Base.:(==)(a::TwoBodyEnergyExpression, b::TwoBodyEnergyExpression) =
    all([ai ∈ b.integrals for ai ∈ a.integrals])

function Base.adjoint(eng::TwoBodyEnergyExpression{A,B}) where {A,B}
    integrals = map(eng.integrals) do integral
        DirectExchangeIntegral{B,A}(integral.b, integral.a)
    end
    TwoBodyEnergyExpression{B,A}(integrals)
end

function Base.show(io::IO, eng::TwoBodyEnergyExpression)
    if isempty(eng.integrals)
        write(io, "0")
        return
    end
    skip = if get(io, :compact, false)
        n = length(eng.integrals)
        if n > 2
            2:n-1
        else
            []
        end
    else
        []
    end
    for (i,integral) in enumerate(eng.integrals)
        if i ∈ skip
            i == skip[1] && write(io, " + …")
            continue
        end
        i > 1 && write(io, " + ")
        write(io, string(integral))
    end
end

function hamiltonian_matrix(::Type{Eng}, ::Type{O},
                            spcs::Vector{<:Configuration{<:SpinOrbital}};
                            selector::Function = peel) where {Eng,O<:SpinOrbital}
    spcs = selector.(spcs)
    m = length(spcs)
    H = Matrix{Eng{O,O}}(undef, m, m)
    for i = 1:m
        for j = 1:m
            H[i,j] = Eng{O,O}(spcs[i],spcs[j])
        end
    end
    H
end

one_body_hamiltonian_matrix(::Type{O}, spcs::Vector{<:Configuration{<:SpinOrbital}}; kwargs...) where {O<:SpinOrbital} =
    hamiltonian_matrix(OneBodyEnergyExpression, O, spcs; kwargs...)

two_body_hamiltonian_matrix(::Type{O}, spcs::Vector{<:Configuration{<:SpinOrbital}}; kwargs...) where {O<:SpinOrbital} =
    hamiltonian_matrix(TwoBodyEnergyExpression, O, spcs; kwargs...)


export OneBodyEnergyExpression, TwoBodyEnergyExpression,
    hamiltonian_matrix, one_body_hamiltonian_matrix, two_body_hamiltonian_matrix
