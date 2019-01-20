struct SlaterDeterminant
    c::Configuration{<:SpinOrbital}
end

function Base.show(io::IO, sd::SlaterDeterminant)
    orbitals = sd.c.orbitals
    n = length(orbitals)
    for (j,p) in enumerate(permutations(1:n))
        s = Combinatorics.parity(p)
        (j > 1 || s == 1) && write(io, " ", s == 0 ? "+" : "-", " ")
        for (i,c) in enumerate(p)
            show(io, orbitals[i])
            write(io, "($(c))")
        end
    end
end

export SlaterDeterminant
