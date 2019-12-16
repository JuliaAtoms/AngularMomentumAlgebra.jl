# Angular Momenta

## [Orbital angular momentum](@id orbital_angular_momentum)

```@docs
OrbitalAngularMomentum
rme(ℓ′, ::OrbitalAngularMomentum, ℓ)
```

## [Spin angular momentum](@id spin_angular_momentum)

```@docs
SpinAngularMomentum
rme(s′, ::SpinAngularMomentum, s)
```

## [Total angular momentum](@id total_angular_momentum)

```@docs
TotalAngularMomentum
rme((ℓ′,s′,J′)::Tuple{<:Number, <:Number, <:Number}, ::TotalAngularMomentum, (ℓ,s,J)::Tuple{<:Number, <:Number, <:Number})
```

