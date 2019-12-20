# Orbitals and quantum systems

When calculating matrix elements of tensor components (using the [The
Wigner–Eckart theorem](@ref wigner_eckart)), it is important to know
which part of the quantum system the tensor acts on,
e.g. [``\tensor{C}``](@ref tensors_spherical_tensors) and
[``\tensor{\nabla}``](@ref gradients) act on the spatial part of a
spin-orbital (the coordinates ``r``, ``\theta``, and ``\phi``, or
equivalently, the quantum numbers ``n``, ``\ell``, ``m_\ell``),
whereas [``\tensor{S}``](@ref spin_angular_momentum) acts on the spin
part (the coordinate ``s``, or equivalently, the quantum numbers ``s``
and ``m_s``).


```@docs
AngularMomentumAlgebra.System
AngularMomentumAlgebra.SubSystem
AngularMomentumAlgebra.FullSystem
AngularMomentumAlgebra.SpatialSubSystem
AngularMomentumAlgebra.OrbitalAngularMomentumSubSystem
AngularMomentumAlgebra.SpinSubSystem
AngularMomentumAlgebra.TotalAngularMomentumSubSystem
```

Given the different systems and subsystems listed above, it is
interesting to access the quantum numbers of an orbital pertaining to
these. For this, `quantum_numbers(::System, ::SpinOrbital)` is provided,
which returns `Pair`s of `(magnitudes,) => projection`, where
`projection` is either a `Number` or `missing` if it is not a good
quantum number, e.g. ``m_\ell`` being the projection quantum number
for both [`AngularMomentumAlgebra.SpatialSubSystem`](@ref) and
[`AngularMomentumAlgebra.OrbitalAngularMomentumSubSystem`](@ref) is
not a good quantum number for coupled spin-orbitals (``\ket{n \ell j
m_j}``).

```jldoctest
julia> using AngularMomentumAlgebra, AtomicLevels, HalfIntegers

julia> o = SpinOrbital(o"3d", 1, -half(1))
3d₁β

julia> ro = SpinOrbital(ro"3d", half(1))
3d(1/2)

julia> quantum_numbers(FullSystem(), o)
((3, 2) => 1, 1/2 => -1/2)

julia> quantum_numbers(FullSystem(), ro)
(3, 2, 1/2, 5/2) => 1/2

julia> quantum_numbers(SpatialSubSystem(), o)
(3, 2) => 1

julia> quantum_numbers(SpatialSubSystem(), ro)
(3, 2) => missing

julia> quantum_numbers(OrbitalAngularMomentumSubSystem(), o)
2 => 1

julia> quantum_numbers(OrbitalAngularMomentumSubSystem(), ro)
2 => missing

julia> quantum_numbers(SpinSubSystem(), o)
1/2 => -1/2

julia> quantum_numbers(SpinSubSystem(), ro)
1/2 => missing

julia> quantum_numbers(TotalAngularMomentumSubSystem(), o)
(2, 1/2, missing) => 1/2

julia> quantum_numbers(TotalAngularMomentumSubSystem(), ro)
(2, 1/2, 5/2) => 1/2
```

```@docs
AngularMomentumAlgebra.quantum_numbers
```
