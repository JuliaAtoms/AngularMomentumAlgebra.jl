"""
    Tensor{k,label}

Abstract base for any tensor of rank `k`.
"""
abstract type Tensor{k,label} end

(::Type{T})(k::Int) where {T<:Tensor} = T{k}()

LinearAlgebra.rank(::Tensor{k}) where k = k

"""
    system(::Tensor)

A general tensor acts on the full system, i.e. all coordinates.
"""
system(::Type{<:Tensor}) = FullSystem()
system(::T) where {T<:Tensor} = system(T)

couples(a::SpinOrbital{<:Orbital}, ::Type{T}, b::SpinOrbital{<:Orbital}) where {T<:Tensor} =
    isequal(other_quantum_numbers(system(T), a, b)...)

couples(a::SpinOrbital{<:RelativisticOrbital}, ::Type{T}, b::SpinOrbital{<:RelativisticOrbital}) where {T<:Tensor} =
    isequal(first.(other_quantum_numbers(system(T), a, b))...)

# This is only true for SO(3)
components(::Tensor{k}) where k = -k:k

function Base.show(io::IO, ::Tensor{k,label}) where {k,label}
    write(io,to_boldface(label)*"̂")
    write(io,"⁽",to_superscript(k),"⁾")
end

# * Tensor components

"""
    TensorComponent(tensor, q)

Represents the `q`th component of a `tensor`; `abs(q) ≤ rank(tensor)`.
"""
struct TensorComponent{T<:Tensor}
    tensor::T
    q::Int
    function TensorComponent(tensor::T, q::Int) where {T<:Tensor}
        k = rank(tensor)
        q ∈ components(tensor) ||
            throw(ArgumentError("Tensor component $(q) not possible for tensor $(tensor) of rank $k"))
        new{T}(tensor, q)
    end
end

function Base.show(io::IO, Tq::TensorComponent)
    show(io, Tq.tensor)
    write(io, to_subscript(Tq.q))
end

Base.parent(Tᵏq::TensorComponent) = Tᵏq.tensor
component(Tᵏq::TensorComponent) = Tᵏq.q

system(::TensorComponent{T}) where T = system(T)

Base.:(==)(a::TensorComponent, b::TensorComponent) =
    a.tensor == b.tensor && a.q == b.q
Base.hash(Tq::TensorComponent, h::UInt) = hash(Tq.tensor, hash(Tq.q, h))

# * Tensor products

@doc raw"""
    TensorProduct{K}(T, U)

A tensor of rank `K` formed from the product of the tensors `T` and
`U`, according to

```math
\begin{equation}
\tag{V3.1.20}
\tensor{X}^{(K)}_Q \equiv
\{\tensor{T}^{(k_1)}\tensor{U}^{(k_2)}\}^{(K)}_Q \defd
\tensor{T}^{(k_1)}_{q_1}
\tensor{U}^{(k_2)}_{q_2}
C_{k_1q_1k_2q_2}^{KQ}
\end{equation}
```
"""
struct TensorProduct{K,A<:Tensor,B<:Tensor} <: Tensor{K,'X'}
    T::A
    U::B
end
TensorProduct(K::Int, T::A, U::B) where {A<:Tensor, B<:Tensor} =
    TensorProduct{K,A,B}(T, U)

function Base.show(io::IO, X::TensorProduct{K}) where K
    write(io,"{")
    show(io, X.T)
    show(io, X.U)
    write(io,"}⁽",to_superscript(K),"⁾")
end

system(::TensorProduct{<:Any,A,B}) where {A,B} = (system(A), system(B))

# ** Scalar product

@doc raw"""
    TensorScalarProduct(T, U)

A tensor of rank 0 (and thus implicitly 0 projection) formed from the
product of the tensors `T` and `U` (which have to have the same rank),
according to

```math
\begin{equation}
\tag{V3.1.30,35}
(\tensor{T}^{(k)} \cdot
\tensor{U}^{(k)}) \defd
(-)^k\angroot{k}
\{\tensor{T}^{(k)}
\tensor{U}^{(k)}\}^{(0)}_0 \equiv
(-)^q
\tensor{T}^{(k)}_{q}
\tensor{U}^{(k)}_{-q}
\end{equation}
```
"""
struct TensorScalarProduct{A<:Tensor,B<:Tensor} <: Tensor{0,'X'}
    T::A
    U::B
    function TensorScalarProduct(T::A, U::B) where {A<:Tensor,B<:Tensor}
        rank(T) == rank(U) ||
            throw(ArgumentError("Cannot form scalar product of tensors of differing ranks"))
        new{A,B}(T,U)
    end
end

function Base.show(io::IO, X::TensorScalarProduct)
    write(io,"(")
    show(io, X.T)
    write(io,"⋅")
    show(io, X.U)
    write(io,")")
end

"""
    dot(T::Tensor, U::Tensor)

Form the scalar product of the two tensors `T` and `U`, which need to
have the same rank.

# Examples

```jldoctest
julia> SphericalTensor(4)⋅SphericalTensor(4)
(𝐂̂⁽⁴⁾⋅𝐂̂⁽⁴⁾)
```
"""
LinearAlgebra.dot(T::Tensor, U::Tensor) =
    TensorScalarProduct(T, U)

system(::TensorScalarProduct{A,B}) where {A,B} = (system(A), system(B))

# * Linear combination of tensors

"""
    LinearCombinationTensor

Represents a linear combination of tensor components.
"""
const LinearCombinationTensor{T<:Tensor,N<:Number} = LinearCombination{<:TensorComponent{<:T},N}

@linearly_combinable TensorComponent

# * OrbitalRadialOverlap

@doc raw"""
    OrbitalRadialOverlap(a,b)

Represents the radial overlap between the orbitals `a` and `b` in a
N-body matrix element expansion. This is different from
`EnergyExpressions.OrbitalOverlap`, which represents integration over
_all_ coordinates. An `OrbitalRadialOverlap` might result when
integrating over the spin–angular degrees of freedom of a matrix
element. As an example, with spin-orbitals on the form
```math
\phi_{n\ell m_\ell s m_s}(\spatialspin) =
\frac{P_{n\ell m_\ell s m_s}(r)}{r}
Y^\ell_m(\theta,\phi)
\chi_{m_s}(s),
\quad
\spatialspin = \{\vec{r},s\},
```
in the spin-restricted case, we have
```math
\int\diff{r}
\conj{P_{\textrm{1s}_0\alpha}}(r)
P_{\textrm{1s}_0\beta}(r) = 1
```
(provided ``P`` is normalized) whereas integration over all
coordinates (spatial and spin)
```math
\int\diff{\spatialspin}
\conj{\phi_{\textrm{1s}_0\alpha}}(\spatialspin)
\phi_{\textrm{1s}_0\beta}(\spatialspin) = 0
```
by symmetry.

# Examples

```jldoctest
julia> AngularMomentumAlgebra.OrbitalRadialOverlap(so"1s₀α", so"1s₀β")
⟨1s₀α|1s₀β⟩ᵣ
```
"""
struct OrbitalRadialOverlap{A,B} <: NBodyTermFactor
    a::A
    b::B
end

# By default, we assume that all orbital radial overlaps are non-zero,
# i.e. the orbitals are non-orthogonal.
Base.iszero(::OrbitalRadialOverlap) = false

Base.:(==)(a::OrbitalRadialOverlap, b::OrbitalRadialOverlap) =
    a.a == b.a && a.b == b.b

Base.hash(o::OrbitalRadialOverlap, h::UInt) =
    hash(hash(hash(o.a), hash(o.b)), h)

Base.adjoint(o::OrbitalRadialOverlap{A,B}) where {A,B} = OrbitalRadialOverlap{B,A}(o.b, o.a)

radial_integral(oo::OrbitalOverlap) = OrbitalRadialOverlap(oo.a, oo.b)

"""
    numbodies(::OrbitalRadialOverlap)

Returns the number of bodies coupled by the zero-body operator in the
orbital overlap, i.e. `0`.
"""
EnergyExpressions.numbodies(::OrbitalRadialOverlap) = 0

"""
    isdependent(o::OrbitalRadialOverlap, orbital)

Returns `true` if the [`OrbitalRadialOverlap`](@ref) `o` depends on `orbital`.

# Examples

```jldoctest
julia> isdependent(OrbitalRadialOverlap(:a,:b), :a)
false

julia> isdependent(OrbitalRadialOverlap(:a,:b), Conjugate(:a))
true

julia> isdependent(OrbitalRadialOverlap(:a,:b), :b)
true
```
"""
EnergyExpressions.isdependent(o::OrbitalRadialOverlap, corb::Conjugate{O}) where O = o.a == corb.orbital
EnergyExpressions.isdependent(o::OrbitalRadialOverlap, orb::O) where O = o.b == orb

Base.show(io::IO, o::OrbitalRadialOverlap) =
    write(io, "⟨$(o.a)|$(o.b)⟩ᵣ")

export Tensor, TensorComponent,
    TensorProduct, TensorScalarProduct,
    system
