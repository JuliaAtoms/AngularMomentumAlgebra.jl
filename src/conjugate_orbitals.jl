import AtomicLevels: AbstractOrbital

struct Conjugate{O<:AbstractOrbital}
    orbital::O
end
function Base.show(io::IO, co::Conjugate{O}) where O
    show(io, co.orbital)
    write(io, "†")
end
Base.conj(o::O) where {O<:AbstractOrbital} = Conjugate(o)

export Conjugate
