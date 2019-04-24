# Definitions

```@meta
CurrentModule = AngularMomentumAlgebra
```

## Notation

The abbreviation $(-)^k\defd(-1)^k$ is used for the roots of negative
unity (see also [`powneg1`](@ref)).

## Spherical harmonics

We assume the following definition of the [spherical harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics)

```math
Y_{m}^{\ell}(\theta,\varphi) = \sqrt{\frac{2\ell+1}{4\pi}\frac{(\ell-m)!}{(\ell+m)!}} P_{\ell}^m(\cos\theta) \mathrm{e}^{\im m \varphi}
```

where ``\theta`` and ``\varphi`` are the usual [spherical
coordinates](https://en.wikipedia.org/wiki/Spherical_coordinate_system) and ``P_\ell^m(z)``
are the [associated Legendre
polynomials](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Definition_for_non-negative_integer_parameters_%E2%84%93_and_m).
The Condon-Shortley phase ``(-1)^m`` is included in the definition of the Legendre
polynomials, consistent with Varshalovich (1988) and the [ISO 80000-2:2009
standard](https://en.wikipedia.org/wiki/ISO_80000-2).

!!! note "Index placement convention"
    Unlike the ISO standard, we put the ``\ell`` index on top and ``m`` on the bottom, to be
    consistent with the way the ``J`` and ``M`` indices are normally written for tensor
    operators.


```@meta
CurrentModule = nothing
```
