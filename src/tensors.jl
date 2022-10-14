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

# * Radial integrals

@doc raw"""
    OrbitalRadialMatrixElement(a,b)

Represents the radial overlap between the orbitals `a` and `b` in a
N-body matrix element expansion. This is different from
`EnergyExpressions.OrbitalMatrixElement`, which represents integration
over _all_ coordinates. An `OrbitalRadialMatrixElement` might result when
integrating over the spin–angular degrees of freedom of a matrix
element.

# Examples

```jldoctest
julia> AngularMomentumAlgebra.OrbitalRadialMatrixElement(so"1s₀α", RadialOperator(), so"2p₀α")
⟨1s₀α|r|2p₀α⟩ᵣ
```
"""
struct OrbitalRadialMatrixElement{A,O,B} <: NBodyTermFactor
    a::A
    o::O
    b::B
end

# By default, we assume that all orbital radial matrix elements are
# non-zero.
Base.iszero(::OrbitalRadialMatrixElement) = false

Base.:(==)(a::OrbitalRadialMatrixElement, b::OrbitalRadialMatrixElement) =
    a.a == b.a && a.o == b.o && a.b == b.b

Base.hash(me::OrbitalRadialMatrixElement, h::UInt) =
    hash(hash(hash(me.a), hash(me.o, hash(me.b))), h)

Base.adjoint(me::OrbitalRadialMatrixElement) = OrbitalRadialMatrixElement(me.b, me.o, me.a)

"""
    numbodies(::OrbitalRadialMatrixElement)

Returns the number of bodies coupled by the operator in the
radial matrix element.
"""
EnergyExpressions.numbodies(me::OrbitalRadialMatrixElement) = length(me.a)

"""
    isdependent(o::OrbitalRadialMatrixElement, orbital)

Returns `true` if the [`OrbitalRadialMatrixElement`](@ref) `o` depends on `orbital`.

# Examples

```jldoctest
julia> isdependent(OrbitalRadialMatrixElement(:a,I,:b), :a)
false

julia> isdependent(OrbitalRadialMatrixElement(:a,I,:b), Conjugate(:a))
true

julia> isdependent(OrbitalRadialMatrixElement(:a,I,:b), :b)
true
```
"""
EnergyExpressions.isdependent(me::OrbitalRadialMatrixElement, corb::Conjugate{O}) where O = corb.orbital ∈ me.a
EnergyExpressions.isdependent(me::OrbitalRadialMatrixElement, orb::O) where O = orb ∈ me.b

function Base.show(io::IO, me::OrbitalRadialMatrixElement)
    write(io, "⟨", join(string.(me.a), " "))
    write(io, "|")
    show(io, me.o)
    write(io, "|")
    write(io, join(string.(me.b), " "), "⟩ᵣ")
end

function radial_integral(a, op, b)
    na = length(a)
    nb = length(b)
    na == nb ||
        throw(DimensionMismatch("Trying to form radial integral from $(na) and $(nb) orbitals"))
    OrbitalRadialMatrixElement(a, op, b)
end

radial_integral(a, op::EnergyExpressions.LinearCombinationOperator, b) =
    sum(c*radial_integral(a, o, b) for (o,c) in op.operators)

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

function radial_integral(a, n::Number, b)
    @assert length(a) == length(b)
    NBodyTerm([OrbitalRadialOverlap(aa, bb) for (aa,bb) in zip(a,b)], n)
end

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

# * Tensor operators

"""
    TensorOperator{N}(T)

Create an `NBodyOperator` out of a tensor component (either a
[`TensorComponent`](@ref) or a [`TensorScalarProduct`](@ref), which
implicitly only has one component, `0`).

# Example

```julia-repl
julia> 𝐉 = TotalAngularMomentum()
𝐉̂⁽¹⁾

julia> 𝐉₀ = TensorComponent(𝐉, 0)
𝐉̂⁽¹⁾₀

julia> A = TensorOperator{1}(𝐉₀)
[𝐉̂⁽¹⁾₀]

julia> cfgs = scs"1s2"
1-element Vector{SpinConfiguration{SpinOrbital{Orbital{Int64}, Tuple{Int64, HalfIntegers.Half{Int64}}}}}:
 1s₀α 1s₀β

julia> Matrix(A, cfgs)
1×1 SparseArrays.SparseMatrixCSC{NBodyMatrixElement, Int64} with 1 stored entry:
 0.5⟨1s₀α|1s₀α⟩ - 0.5⟨1s₀β|1s₀β⟩
```
"""
struct TensorOperator{N,Tt} <: NBodyOperator{N}
    T::Tt
    TensorOperator{N}(T::Tt) where {N,Tt} =
        new{N,Tt}(T)
end

Base.hash(t::TensorOperator{N}, h::UInt) where N = hash(N, hash(t.T, h))

spin_ang_coeff(me::OrbitalMatrixElement{1,<:SpinOrbital,<:TensorOperator{1},<:SpinOrbital}) =
    dot(me.a[1], me.o.T, me.b[1])

spin_ang_coeff(me::OrbitalMatrixElement{<:Any,<:SpinOrbital,<:TensorOperator,<:SpinOrbital}) =
    dot(Tuple(me.a), me.o.T, Tuple(me.b))

Base.iszero(me::OrbitalMatrixElement{<:Any,<:SpinOrbital,<:TensorOperator,<:SpinOrbital}) =
    iszero(spin_ang_coeff(me))

function integrate_spinor(me::OrbitalMatrixElement{N,<:SpinOrbital,<:TensorOperator{N},<:SpinOrbital}) where N
    coeff = spin_ang_coeff(me)
    iszero(coeff) && return zero(NBodyMatrixElement)
    NBodyMatrixElement(radial_integral(me.a, coeff, me.b))
end

function Base.show(io::IO, to::TensorOperator)
    write(io, "[")
    show(io, to.T)
    write(io, "]")
end

function Base.show(io::IO, to::TensorOperator{1,<:TensorScalarProduct})
    T = to.T.T
    U = to.T.U
    write(io, "[")
    show(io, T)
    if T ≠ U
        write(io, "⋅")
        show(io, U)
    end
    write(io, "]")
    T == U && write(io, "²")
end

function Base.show(io::IO, to::TensorOperator{2,<:TensorScalarProduct})
    T = to.T.T
    U = to.T.U
    write(io, "[")
    show(io, T)
    write(io, "(1)⋅")
    show(io, U)
    write(io, "(2)]")
end

@doc raw"""
    many_electron_scalar_product(𝐓::Tensor{k}, 𝐔::Tensor{k}=𝐓) where k

Create the total tensor acting on a many-electron state according to

```math
\begin{equation}
\begin{aligned}
(\tensor{T}^{(k)} \cdot
\tensor{U}^{(k)})
&=
\sum_{i,j}
[\tensor{t}^{(k)}(i)
 \cdot
\tensor{u}^{(k)}(j)] \\
&\equiv
\sum_i
[\tensor{t}^{(k)}(i)
 \cdot
\tensor{u}^{(k)}(i)] +
\sum_{i\ne j}
2[\tensor{t}^{(k)}(i)
 \cdot
\tensor{u}^{(k)}(j)],
\end{aligned}
\tag{H11-32*}
\end{equation}
```
where ``\tensor{t}^{(k)}(i)`` and ``\tensor{u}^{(k)}(j)`` act only on
electron ``i`` and ``j``, respectively [cf. John E. Harriman:
_Theoretical Foundations of Electron Spin Resonance_ (1978); note that
Harriman uses another normalization of the ladder operators compared
to ours: [Eq. (V3.1.1)](@ref angular_momenta), which explains why his
Eq. (H11-32) is missing a factor of ``2``].

# Examples

```julia-repl
julia> 𝐉 = TotalAngularMomentum()
𝐉̂⁽¹⁾

julia> A = many_electron_scalar_product(𝐉)
[𝐉̂⁽¹⁾]² + 2.0[𝐉̂⁽¹⁾(1)⋅𝐉̂⁽¹⁾(2)]

julia> 𝐋 = OrbitalAngularMomentum()
𝐋̂⁽¹⁾

julia> 𝐒 = SpinAngularMomentum()
𝐒̂⁽¹⁾

julia> B = many_electron_scalar_product(𝐋, 𝐒)
[𝐋̂⁽¹⁾⋅𝐒̂⁽¹⁾] + 2.0[𝐋̂⁽¹⁾(1)⋅𝐒̂⁽¹⁾(2)]

julia> cfgs = rscs"1s2"
1-element Vector{SpinConfiguration{SpinOrbital{RelativisticOrbital{Int64}, Tuple{HalfIntegers.Half{Int64}}}}}:
 1s(-1/2) 1s(1/2)

julia> Matrix(A, cfgs)
1×1 SparseArrays.SparseMatrixCSC{NBodyMatrixElement, Int64} with 1 stored entry:
 0.75⟨1s(-1/2)|1s(-1/2)⟩ + 0.75⟨1s(1/2)|1s(1/2)⟩ - ⟨1s(-1/2)|1s(1/2)⟩⟨1s(1/2)|1s(-1/2)⟩ - 0.5⟨1s(-1/2)|1s(-1/2)⟩⟨1s(1/2)|1s(1/2)⟩

julia> Matrix(B, cfgs)
1×1 SparseArrays.SparseMatrixCSC{NBodyMatrixElement, Int64} with 0 stored entries:
 ⋅
```
"""
function many_electron_scalar_product(𝐓::Tensor{k}, 𝐔::Tensor{k}=𝐓) where k
    𝐗 = 𝐓⋅𝐔
    TensorOperator{1}(𝐗) + 2TensorOperator{2}(𝐗)
end

export Tensor, TensorComponent,
    TensorProduct, TensorScalarProduct,
    system,
    TensorOperator, many_electron_scalar_product
