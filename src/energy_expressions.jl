Base.Matrix(::QuantumOperator, csfs::Vector{<:CSF},
            overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[]) =
    throw(ArgumentError("Derivation of energy expression for atomic CSFs not yet implemented"))
# * Spin-configurations

"""
    Matrix(op::QuantumOperator,
           spin_cfgs::Vector{<:Configuration{<:SpinOrbital}}[, overlaps])

Generate the energy-expression associated with the quantum operator
`op`, in the basis of the spin-configurations `spin_cfgs`, with an
optional set of orbital `overlaps`, specifying any desired
non-orthogonalities. The energy expression is generated in a
basis-agnostic way by EnergyExpressions.jl and then
multipole-expanded.
"""
function Base.Matrix(op::QuantumOperator, spin_cfgs::VSC,
                     overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[]) where {VSC<:AbstractVector{<:Configuration{<:SpinOrbital}}}
    M = Matrix(op, SlaterDeterminant.(spin_cfgs), overlaps)
    m,n = size(M)
    Mm = spzeros(eltype(M), m, n)
    for i = 1:m
        for j = 1:n
            me = transform(multipole_expand, M[i,j])
            iszero(me) && continue
            Mm[i,j] = me
        end
    end
    Mm
end

# * Configurational averages

function NBodyMatrixElement(a::Configuration, op::OneBodyOperator, b::Configuration, ::Nothing)
    terms = NBodyTerm[]
    if a == b
        for (orb,q,_) in a
            me = q*OrbitalMatrixElement([orb], op, [orb])
            !iszero(me) && push!(terms, me)
        end
    else
        @warn "Configurational averages only make sense for a single configuration"
    end

    NBodyMatrixElement(terms)
end

"""
    f_av(ℓ)

Average interactions between equivalent electrons of orbital angular
momentum `ℓ` (belonging to the same subshell). Returns multipole
orders and associated expansion coefficients.

# Examples

```jldoctest
julia> collect(AngularMomentumAlgebra.f_av(3))
3-element Array{Tuple{Int64,Float64},1}:
 (2, -0.02051282051282051)
 (4, -0.013986013986013984)
 (6, -0.017930787161556393)
```
"""
function f_av(ℓ)
    # The coefficients are those in front of Fᵏ(ii) in Eq. (6.39) of
    # Cowan (1981).
    μ = (2ℓ+1)/(4ℓ+1)
    ((k,-μ*wigner3j(ℓ,k,ℓ,0,0,0)^2)
     for k ∈ 2:2:2ℓ)
end

"""
    g_av(ℓ,ℓ′)

Average interactions between non-equivalent electrons of orbital
angular momenta `ℓ` and `ℓ′` (not belonging to the same
subshell). Returns multipole orders and associated expansion
coefficients.

# Examples

```jldoctest
julia> collect(AngularMomentumAlgebra.g_av(1,2))
2-element Array{Tuple{Int64,Float64},1}:
 (1, -0.06666666666666665)
 (3, -0.042857142857142864)
```
"""
function g_av(ℓ,ℓ′)
    # The coefficients are those in front of Gᵏ(ij) in Eq. (6.38) of
    # Cowan (1981).
    ((k,-(wigner3j(ℓ,k,ℓ′,0,0,0)^2)/2)
     for k ∈ triangle_range(ℓ,ℓ′))
end

function NBodyMatrixElement(a::Configuration, op::CoulombInteraction, b::Configuration, ::Nothing)
    g = k -> CoulombInteractionMultipole(k,op)
    R = (k,a,b,c,d) -> OrbitalMatrixElement([a,b], g(k), [c,d])
    F = (k,a,b) -> R(k,a,b,a,b)
    G = (k,a,b) -> R(k,a,b,b,a)
    terms = NBodyTerm[]
    if a == b
        # Fischer (1977), Eq. (2-2)
        #
        # Keep occupancies in energy expression to yield correct total
        # energy; divide equations of motion by degeneracy of orbital.
        for (i,(orbᵢ,qᵢ,_)) in enumerate(a)
            # Equivalent electrons
            μ = qᵢ*(qᵢ-1)/2
            iszero(μ) || push!(terms, μ*F(0,orbᵢ,orbᵢ))
            ℓᵢ = orbᵢ.ℓ
            for (k,fₖ) in f_av(ℓᵢ)
                me = μ*fₖ*F(k,orbᵢ,orbᵢ)
                iszero(me) || push!(terms, me)
            end

            # Non-equivalent electrons
            for (j,(orbⱼ,qⱼ,_)) in enumerate(b)
                j == i && break
                ν = qᵢ*qⱼ
                push!(terms, ν*F(0,orbᵢ,orbⱼ))
                ℓⱼ = orbⱼ.ℓ
                for (k,gₖ) in g_av(ℓᵢ,ℓⱼ)
                    me = ν*gₖ*G(k,orbᵢ,orbⱼ)
                    iszero(me) || push!(terms, me)
                end
            end
        end
    else
        @warn "Configurational averages only make sense for a single configuration"
    end

    NBodyMatrixElement(terms)
end



function overlap_matrix(a::Wfn, b::Wfn, overlaps) where {Wfn<:Union{Configuration,CSF}}
    isempty(overlaps) ||
        throw(ArgumentError("Non-orthogonal orbitals not yet supported for $(Wfn)s"))
    nothing
end
