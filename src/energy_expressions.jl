# * Integration of spinors

"""
    integrate_spinor(me)

Perform the spin-angular integration of the matrix element `me`,
leaving only a radial integral multiplied by a spin-angular
coefficient. The spin-angular integral is dependent on the specific
combination of spin-orbitals and the operator (expressed as a tensor);
the default implementation is to leave `me` as-is, corresponding to a
spin-angular integral of unity.
"""
integrate_spinor(me) = me

"""
    integrate_spinor(me::OrbitalMatrixElement{2,<:Any,<:CoulombInteraction,<:Any})

Perform the spin-angular integration of the two-electron matrix
element `me`, by first multipole-expanding the Coulomb interaction and
then integrating all the resulting terms over the spin-angular
coordinates (see [`multipole_expand`](@ref)).
"""
integrate_spinor(me::OrbitalMatrixElement{2,<:Any,<:CoulombInteraction,<:Any}) =
    multipole_expand(me)


"""
    integrate_spinor(integral::NBodyTermFactor)

Dummy method that returns `integral` unchanged, used for all
`NBodyTermFactor`s that are _not_ to be multipole-expanded.
"""
integrate_spinor(integral::NBodyTermFactor) = NBodyMatrixElement([integral])

# * Energy Expressions

# ** Spin-configurations

"""
    Matrix(op::QuantumOperator,
           spin_cfgs::Vector{<:Configuration{<:SpinOrbital}}[, overlaps])

Generate the energy-expression associated with the quantum operator
`op`, in the basis of the spin-configurations `spin_cfgs`, with an
optional set of orbital `overlaps`, specifying any desired
non-orthogonalities. The energy expression is generated in a
basis-agnostic way by EnergyExpressions.jl and each term is then
integrated over the spin-angular coordinates using
[`integrate_spinor`](@ref).
"""
function Base.Matrix(op::QuantumOperator, spin_cfgs::VSC,
                     overlaps::Vector{<:OrbitalOverlap}=OrbitalOverlap[];
                     kwargs...) where {VSC<:AbstractVector{<:Configuration{<:SpinOrbital}}}
    bcs = BitConfigurations(spin_cfgs, overlaps)
    E = Matrix(bcs, op; kwargs...)
    transform(integrate_spinor, E; kwargs...)
end

# ** Average-of-configuration

function NBodyMatrixElement(a::Configuration{O}, op::OneBodyOperator, b::Configuration{O}, overlap) where {O<:Union{Orbital,RelativisticOrbital}}
    a == b || return 0
    overlap == I || @error "We don't know what to do about non-orthogonal orbitals" overlap
    terms = NBodyTerm[]
    for (o,w,_) in b
        me = integrate_spinor(OrbitalMatrixElement([o], op, [o]))
        nbme = convert(NBodyMatrixElement, w*me)
        append!(terms, nbme.terms)
    end
    NBodyMatrixElement(terms)
end

function NBodyMatrixElement(ac::Configuration{O}, op::CoulombInteraction, bc::Configuration{O}, overlap) where {O<:Union{Orbital,RelativisticOrbital}}
    ac == bc || return 0
    overlap == I || @error "We don't know what to do about non-orthogonal orbitals" overlap
    # Eq. (2-2) of
    #
    # - Froese Fischer, C. (1977). The Hartree–Fock Method for Atoms : A
    #   Numerical Approach. New York: Wiley.
    #
    # and Eq. (2.40) of
    #
    # - Froese Fischer, C., Brage, T., & Jönsson, P. (1997). Computational
    #   Atomic Structure : An MCHF Approach. Bristol, UK Philadelphia, Penn:
    #   Institute of Physics Publ.

    terms = NBodyTerm[]
    add!(v) = append!(terms, convert(NBodyMatrixElement, v).terms)
    for (a,w,_) in bc
        w == 1 && continue
        μ = w*(w-1)/2
        ℓa = a.ℓ
        for k = 0:2ℓa
            me = radial_integral([a,a], (k, op), [a,a])
            f = k == 0 ? 1 : (-(2ℓa+1)/(4ℓa+1))*wigner3j(ℓa, k, ℓa,
                                                         0,  0, 0)^2
            add!(μ*f*me)
        end
    end
    m = length(bc)
    for i=2:m
        (a,wa,_) = bc[i]
        ℓa = a.ℓ
        for j = 1:i-1
            (b,wb,_) = bc[j]
            ℓb = b.ℓ
            μ = wa*wb
            add!(μ*radial_integral([a,b], (0, op), [a,b]))
            for k = abs(ℓa-ℓb):(ℓa+ℓb)
                me = radial_integral([a,b], (k, op), [b,a])
                g = -1/2*wigner3j(ℓa, k, ℓb,
                                  0,  0, 0)^2
                add!(μ*g*me)
            end
        end
    end
    NBodyMatrixElement(terms)
end
