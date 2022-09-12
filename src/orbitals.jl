"""
    System

The abstract base type for a (quantum) system.
"""
abstract type System end

"""
    SubSystem{which}

A subsystem, i.e. limited to a set of coordinates.
"""
struct SubSystem{which} <: System end

"""
    quantum_numbers(system, a, b)

Return the quantum numbers characterizing `system` for the orbitals `a` and `b`.
"""
quantum_numbers(system, a, b) = quantum_numbers(system, a), quantum_numbers(system, b)

"""
    other_quantum_numbers(system, a, b)

Return the quantum numbers characterizing the orthogonal complement to
`system` for the orbitals `a` and `b`.
"""
other_quantum_numbers(system, a, b) = other_quantum_numbers(system, a), other_quantum_numbers(system, b)

# * Full system

"""
    FullSystem

The entire system, i.e. all coordinates.
"""
struct FullSystem <: System end

"""
    quantum_numbers(::FullSystem, o::SpinOrbital{<:Orbital})

The full system of an uncoupled spin-orbital is ``n\\ell m_\\ell; s
m_s``, where ``;`` denotes that the spatial and spin subsystems are
separable.
"""
quantum_numbers(::FullSystem, o::SpinOrbital{<:Orbital}) =
    quantum_numbers(SpatialSubSystem(), o), quantum_numbers(SpinSubSystem(), o)

"""
    quantum_numbers(::FullSystem, o::SpinOrbital{<:RelativisticOrbital})

The full system of a coupled spin-orbital is ``n\\ell s j m_j``.
"""
quantum_numbers(::FullSystem, o::SpinOrbital{<:RelativisticOrbital}) =
    ((o.orb.n, o.orb.ℓ, half(1), o.orb.j), o.m[1])

"""
    other_quantum_numbers(::FullSystem, ::SpinOrbital)

No quantum numbers characterize the orthogonal complement to
[`FullSystem`](@ref).
"""
other_quantum_numbers(::FullSystem, ::SpinOrbital) = ((),missing)

# * Spatial subsystem

"""
    SpatialSubSystem

The spatial subsystem, i.e. the coordinates `r`, `θ`, and `ϕ`.
"""
const SpatialSubSystem = SubSystem{:spatial}

"""
    quantum_numbers(::SpatialSubSystem, o::SpinOrbital{<:Orbital})

The spatial subsystem of an uncoupled spin-orbital is ``n\\ell m_\\ell``.
"""
quantum_numbers(::SpatialSubSystem, o::SpinOrbital{<:Orbital}) =
    ((o.orb.n, o.orb.ℓ), o.m[1])

"""
    quantum_numbers(::SpatialSubSystem, o::SpinOrbital{<:RelativisticOrbital})

The spatial subsystem of a coupled spin-orbital is just ``n\\ell
m_\\ell``; ``m_\\ell`` is not a good quantum number.
"""
quantum_numbers(::SpatialSubSystem, o::SpinOrbital{<:RelativisticOrbital}) =
    ((o.orb.n, o.orb.ℓ), missing)

"""
    other_quantum_numbers(::SpatialSubSystem, o::SpinOrbital{<:Orbital})

The orthogonal complement to [`SpatialSubSystem`](@ref) is
[`SpinSubSystem`](@ref), which is characterized by ``sm_s``.
"""
other_quantum_numbers(::SpatialSubSystem, o::SpinOrbital{<:Orbital}) =
    (half(1), o.m[2])

"""
    other_quantum_numbers(::SpatialSubSystem, o::SpinOrbital{<:RelativisticOrbital})

The orthogonal complement to [`SpatialSubSystem`](@ref) is
[`SpinSubSystem`](@ref), which is characterized by ``s``; its
projection is not a good quantum number in the coupled basis.
"""
other_quantum_numbers(::SpatialSubSystem, o::SpinOrbital{<:RelativisticOrbital}) =
    (half(1), missing)

# * Orbital angular momentum subsystem

"""
    OrbitalAngularMomentaSubSystem

The orbital angular momentum subsystem, i.e. the coordinates, `θ` and
`ϕ`.
"""
const OrbitalAngularMomentumSubSystem = SubSystem{:orbital_angular_momentum}

"""
    quantum_numbers(::OrbitalAngularMomentumSubSystem, o::SpinOrbital{<:Orbital})

The orbital angular momentum subsystem of an uncoupled spin-orbital is ``\\ell m_\\ell``.
"""
quantum_numbers(::OrbitalAngularMomentumSubSystem, o::SpinOrbital{<:Orbital}) =
    (o.orb.ℓ, o.m[1])

"""
    quantum_numbers(::OrbitalAngularMomentumSubSystem, o::SpinOrbital{<:RelativisticOrbital})

The orbital angular momentum subsystem of a coupled spin-orbital is
just ``\\ell``; ``m_\\ell`` is not a good quantum number.
"""
quantum_numbers(::OrbitalAngularMomentumSubSystem, o::SpinOrbital{<:RelativisticOrbital}) =
    (o.orb.ℓ, missing)

"""
    other_quantum_numbers(::OrbitalAngularMomentumSubSystem, o::SpinOrbital{<:Orbital})

The orthogonal complement to [`OrbitalAngularMomentumSubSystem`](@ref)
is characterized by the quantum numbers ``n s m_s``; however, ``n``
does not affect the matrix elements of ``𝐋``, so it is silently
ignored.
"""
other_quantum_numbers(::OrbitalAngularMomentumSubSystem, o::SpinOrbital{<:Orbital}) =
    ((half(1),), o.m[2])

"""
    other_quantum_numbers(::OrbitalAngularMomentumSubSystem, o::SpinOrbital{<:RelativisticOrbital})

The orthogonal complement to [`OrbitalAngularMomentumSubSystem`](@ref)
is characterized by the quantum numbers ``n s``; however, ``n`` does
not affect the matrix elements of ``𝐋``, so it is silently
ignored. Additoinally, the projection of ``s`` is not a good quantum
number in the coupled basis.
"""
other_quantum_numbers(::OrbitalAngularMomentumSubSystem, o::SpinOrbital{<:RelativisticOrbital}) =
    ((half(1),), missing)

const SpatialSubSystems = Union{SpatialSubSystem,OrbitalAngularMomentumSubSystem}

# * Spin subsystem

"""
    SpinSubSystem

The spin subsystem, i.e. the coordinate `s`.
"""
const SpinSubSystem = SubSystem{:spin}

"""
    quantum_numbers(::SpinSubSystem, o::SpinOrbital{<:Orbital})

The spin subsystem of an uncoupled spin-orbital is ``s m_s``.
"""
quantum_numbers(::SpinSubSystem, o::SpinOrbital{<:Orbital}) =
    (half(1), o.m[2])

"""
    quantum_numbers(::SpinSubSystem, o::SpinOrbital{<:RelativisticOrbital})

The spin subsystem of a coupled spin-orbital is just ``s``; ``m_s`` is
not a good quantum number.
"""
quantum_numbers(::SpinSubSystem, o::SpinOrbital{<:RelativisticOrbital}) =
    (half(1), missing)

"""
    other_quantum_numbers(::SpinSubSystem, o::SpinOrbital{<:Orbital})

The orthogonal complement to [`SpinSubSystem`](@ref) is
[`SpatialSubSystem`](@ref) which is characterized by the quantum
numbers ``n \\ell m_\\ell``; however, ``n`` does not affect the matrix
elements of ``𝐒``, so it is silently ignored.
"""
other_quantum_numbers(::SpinSubSystem, o::SpinOrbital{<:Orbital}) =
    ((o.orb.ℓ,), o.m[1])

"""
    other_quantum_numbers(::SpinSubSystem, o::SpinOrbital{<:RelativisticOrbital})

The orthogonal complement to [`SpinSubSystem`](@ref) is
[`SpatialSubSystem`](@ref) which is characterized by the quantum
numbers ``n \\ell``; however, ``n`` does not affect the matrix
elements of ``𝐒``, so it is silently ignored. Additionally, the
projection of ``\\ell`` is not a good quantum number in the coupled
basis.
"""
other_quantum_numbers(::SpinSubSystem, o::SpinOrbital{<:RelativisticOrbital}) =
    ((o.orb.ℓ,), missing)

# * Total angular momentum subsystem

"""
    TotalAngularMomentumSubSystem

The total angular momentum subsystem, i.e. the coordinates, `θ`, `ϕ`, and
`s`.
"""
const TotalAngularMomentumSubSystem = SubSystem{:total_angular_momentum}

"""
    quantum_numbers(::TotalAngularMomentumSubSystem, o::SpinOrbital{<:Orbital})

The total angular momentum of an uncoupled spin-orbital is not a good
quantum number; only its projection is known. The system is specified
by ``\\ell m_\\ell; s m_s``, where ``;`` denotes that the spatial and
spin subsystems are separable.
"""
quantum_numbers(::TotalAngularMomentumSubSystem, o::SpinOrbital{<:Orbital}) =
    quantum_numbers(OrbitalAngularMomentumSubSystem(), o), quantum_numbers(SpinSubSystem(), o)

"""
    quantum_numbers(::TotalAngularMomentumSubSystem, o::SpinOrbital{<:RelativisticOrbital})

The total angular momentum subsystem of a coupled spin-orbital just
``\\ell s j m_j``.
"""
quantum_numbers(::TotalAngularMomentumSubSystem, o::SpinOrbital{<:RelativisticOrbital}) =
    ((o.orb.ℓ, half(1), o.orb.j), o.m[1])

"""
    other_quantum_numbers(::TotalAngularMomentumSubSystem, o::SpinOrbital{<:Orbital})

The orthogonal complement to [`TotalAngularMomentumSubSystem`](@ref)
is characterized by the principal quantum number ``n``; however, this
does not affect the matrix elements of ``𝐉``, so it is silently
ignored.
"""
other_quantum_numbers(::TotalAngularMomentumSubSystem, o::SpinOrbital{<:Orbital}) =
    ((), missing)

"""
    other_quantum_numbers(::TotalAngularMomentumSubSystem, o::SpinOrbital{<:RelativisticOrbital})

The orthogonal complement to [`TotalAngularMomentumSubSystem`](@ref)
is characterized by the principal quantum number ``n``; however, this
does not affect the matrix elements of ``𝐉``, so it is silently
ignored.
"""
other_quantum_numbers(::TotalAngularMomentumSubSystem, o::SpinOrbital{<:RelativisticOrbital}) =
    ((), missing)


export quantum_numbers, other_quantum_numbers,
    SubSystem, FullSystem,
    SpatialSubSystem,
    OrbitalAngularMomentumSubSystem, SpinSubSystem,
    TotalAngularMomentumSubSystem
