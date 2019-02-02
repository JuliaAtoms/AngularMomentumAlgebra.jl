#=

Some references below to

- Cook, D. (2005). Handbook of computational quantum
  chemistry. Mineola, N.Y: Dover Publications, ISBN: 9780486443072.

=#

# * One-body energy expression

mutable struct OneBodyEnergyExpression{A<:SpinOrbital,B<:SpinOrbital}
    integrals::Vector{OneBodyIntegral{A,B}}
    signs::Vector{Int}
    OneBodyEnergyExpression{A,B}(integrals::Vector{OneBodyIntegral{A,B}},
                                 signs::Vector{Int}) where {A,B} =
                                     new{A,B}(integrals,signs)
    function OneBodyEnergyExpression{A,B}(a::Configuration{<:SpinOrbital},
                                          b::Configuration{<:SpinOrbital};
                                          orthogonal=true) where {A,B}
        m = length(a)
        n = length(b)
        integrals = Vector{OneBodyIntegral{A,B}}()
        signs = Vector{Int}()

        if orthogonal # Cook 2005, p. 72--73
            sab = substitutions(a,b)
            if isempty(sab) # Case 1, a == b
                for i = 1:m
                    ai = a.orbitals[i]
                    push!(integrals, OneBodyIntegral{A,B}(ai,ai))
                    push!(signs, 1)
                end
            elseif length(sab) == 1 # Case 2, one substitution
                push!(integrals, OneBodyIntegral{A,B}(sab[1]...))
                push!(signs, 1)
            end
        else # Cook 2005, p. 63
            for i = 1:m
                ai = a.orbitals[i]
                for j = 1:n
                    push!(integrals, OneBodyIntegral{A,B}(ai,b.orbitals[j]))
                    push!(signs, (-1)^(i+j))
                end
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

Base.iszero(eng::OneBodyEnergyExpression) = isempty(eng.integrals)

Base.copy(eng::OneBodyEnergyExpression{A,B}) where {A,B} =
    OneBodyEnergyExpression{A,B}(copy(eng.integrals), copy(eng.signs))

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

function adjacent_in_coord!(comp::Function, eng::OneBodyEnergyExpression, coord::Symbol)
    pred = i -> comp(getproperty(i.a, coord), getproperty(i.b, coord))
    sel = map(pred, eng.integrals)
    eng.integrals = eng.integrals[sel]
    eng.signs = eng.signs[sel]
    eng
end

adjacent_in_coord(comp::Function, eng::OneBodyEnergyExpression, coord::Symbol) =
    adjacent_in_coord!(comp, copy(eng), coord)

diagonal_in_coord!(eng::OneBodyEnergyExpression, coord::Symbol) =
    adjacent_in_coord!(isequal, eng, coord)

diagonal_in_coord(eng::OneBodyEnergyExpression, coord::Symbol) =
    diagonal_in_coord!(copy(eng), coord)

# * Two-body energy expression

mutable struct TwoBodyEnergyExpression{A<:SpinOrbital,B<:SpinOrbital,
                                       C<:SpinOrbital,D<:SpinOrbital}
    integrals::Vector{TwoBodyIntegral{A,B,C,D}}
    # signs::Vector{Int} # Needed for the non-orthogonal case
    TwoBodyEnergyExpression{A,B,C,D}(integrals::Vector{TwoBodyIntegral{A,B,C,D}}) where {A,B,C,D} =
                                     new{A,B,C,D}(integrals)

    function TwoBodyEnergyExpression{A,B,C,D}(a::Configuration{<:SpinOrbital},
                                              b::Configuration{<:SpinOrbital};
                                              orthogonal=true) where {A,B,C,D}
        m = length(a)

        integrals = Vector{TwoBodyIntegral{A,B,C,D}}()
        if orthogonal
            sab = substitutions(a,b)
            if isempty(sab) # Case 1, a == b
                for i = 1:m
                    ai = a.orbitals[i]
                    for j = 1:i-1
                        push!(integrals, TwoBodyIntegral{A,B,C,D}(ai,b.orbitals[j]))
                        # push!(signs, 1)
                    end
                end
            elseif length(sab) == 1 # Case 2, one substitution
                ai,bi′ = sab[1]
                for k = 1:m
                    ak = a.orbitals[k]
                    ak == ai && continue
                    push!(integrals, TwoBodyIntegral{A,B,C,D}(ai,ak,bi′))
                    # push!(signs, 1)
                end
            elseif length(sab) == 2 # Case 3, two substitutions
                push!(integrals, TwoBodyIntegral{A,B,C,D}(first.(sab)...,last.(sab)...))
                # push!(signs, 1)
            end
        else # Cook 2005, p. 72
            error("Not yet implemented")
            # for i = 1:m
            #     ai = a.orbitals[i]

            #     append!(integrals,
            #             map(j -> TwoBodyIntegral{A,B}(ai,b.orbitals[j]), 1:i-1))
            # end
        end
        new{A,B,C,D}(integrals)
    end
end

function TwoBodyEnergyExpression(a::Configuration{<:SpinOrbital},
                                 b::Configuration{<:SpinOrbital})
    A = promote_type(typeof.(a.orbitals)...)
    B = promote_type(typeof.(b.orbitals)...)
    TwoBodyEnergyExpression{A,A,B,B}(a,b)
end

Base.:(==)(a::TwoBodyEnergyExpression, b::TwoBodyEnergyExpression) =
    all([ai ∈ b.integrals for ai ∈ a.integrals])

Base.iszero(eng::TwoBodyEnergyExpression) = isempty(eng.integrals)

function Base.adjoint(eng::TwoBodyEnergyExpression{A,B,C,D}) where {A,B,C,D}
    integrals = map(eng.integrals) do integral
        TwoBodyIntegral{D,C,B,A}(integral.b, integral.a)
    end
    TwoBodyEnergyExpression{D,C,B,A}(integrals)
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

function diagonal_in_coord!(eng::TwoBodyEnergyExpression, coord::Symbol)
    pred = i -> getproperty(i.a, coord) == getproperty(i.b, coord) == getproperty(i.c, coord) == getproperty(i.d, coord)
    filter!(pred, eng.integrals)
end

diagonal_in_coord(eng::TwoBodyEnergyExpression, coord::Symbol) =
    diagonal_in_coord!(copy(eng), coord)

# * Hamiltonian matrix generation

function hamiltonian_matrix(::Type{Eng},
                            spcs::Vector{<:Configuration{<:SpinOrbital}};
                            selector::Function = peel) where {Eng<:Union{<:OneBodyEnergyExpression,<:TwoBodyEnergyExpression}}
    spcs = selector.(spcs)
    m = length(spcs)
    H = Matrix{Eng}(undef, m, m)
    for i = 1:m
        for j = 1:m
            H[i,j] = Eng(spcs[i],spcs[j])
        end
    end
    H
end

one_body_hamiltonian_matrix(::Type{O}, spcs::Vector{<:Configuration{<:SpinOrbital}}; kwargs...) where {O<:SpinOrbital} =
    hamiltonian_matrix(OneBodyEnergyExpression{O,O}, spcs; kwargs...)

two_body_hamiltonian_matrix(::Type{O}, spcs::Vector{<:Configuration{<:SpinOrbital}}; kwargs...) where {O<:SpinOrbital} =
    hamiltonian_matrix(TwoBodyEnergyExpression{O,O,O,O}, spcs; kwargs...)

export OneBodyEnergyExpression, TwoBodyEnergyExpression,
    adjacent_in_coord!, diagonal_in_coord!,
    adjacent_in_coord, diagonal_in_coord,
    hamiltonian_matrix, one_body_hamiltonian_matrix, two_body_hamiltonian_matrix
