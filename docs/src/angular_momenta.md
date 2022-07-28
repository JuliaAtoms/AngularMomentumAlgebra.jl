# [Angular Momenta](@id angular_momenta)

An irreducible tensor ``\tensor{T}^{(J)}_M`` obeys the following
commutation relations with the angular momentum ``\tensor{J}``:
```math
\begin{equation}
\begin{aligned}
{}[J_{\pm},\tensor{T}^{(J)}_M]
&=
\sqrt{\frac{J(J+1)-M(M\pm1)}{2}}\tensor{T}^{(J)}_{M\pm1}, \\
[J_0,\tensor{T}^{(J)}_M]
&=
M\tensor{T}^{(J)}_M.
\end{aligned}
\label{eqn:ladder-operators}
\tag{V3.1.1}
\end{equation}
```

## [Orbital angular momentum](@id orbital_angular_momentum)

```@docs
OrbitalAngularMomentum
rme(ℓ′, ::OrbitalAngularMomentum, ℓ)
AngularMomentumAlgebra.couplings(tensor::OrbitalAngularMomentum, ℓ)
```

## [Spin angular momentum](@id spin_angular_momentum)

```@docs
SpinAngularMomentum
rme(s′, ::SpinAngularMomentum, s)
AngularMomentumAlgebra.couplings(tensor::SpinAngularMomentum, s)
```

## [Total angular momentum](@id total_angular_momentum)

```@docs
TotalAngularMomentum
rme((ℓ′,s′,J′)::Tuple{<:Number, <:Number, <:Number}, ::TotalAngularMomentum, (ℓ,s,J)::Tuple{<:Number, <:Number, <:Number})
AngularMomentumAlgebra.couplings(tensor::TotalAngularMomentum, (ℓ,s,J)::Tuple{<:Number, <:Number, <:Number})
```

