# Definitions

```@meta
CurrentModule = AngularMomentumAlgebra
```

This page defines much of the general notation and conventions used in the code. Where possible,
we are consistent with Varshalovich (1988).

## Shorthands

The abbreviation ``(-)^k\defd(-1)^k`` is used for the roots of negative unity (see also [`powneg1`](@ref)).

A commonly occurring factor in angular momentum algebra is
```math
\angroot{j_1j_2...j_n}
\defd[(2j_1+1)(2j_2+1)...(2j_n+1)]^{1/2}.
\tag{V13.1.3½}
```
It can be calculated with the unexported [`AngularMomentumAlgebra.∏`](@ref) function.

Indices appearing in pairs on only one side of an equation are
[implicitly summed
over](https://en.wikipedia.org/wiki/Einstein_notation).

## Spherical harmonics

We assume the following definition of the [spherical harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics)

```math
Y_{m}^{\ell}(\theta,\varphi) = \sqrt{\frac{2\ell+1}{4\pi}\frac{(\ell-m)!}{(\ell+m)!}} P_{\ell}^m(\cos\theta) \mathrm{e}^{\im m \varphi}
\tag{V5.2.1}
```

where ``\theta`` and ``\varphi`` are the usual [spherical
coordinates](https://en.wikipedia.org/wiki/Spherical_coordinate_system) and ``P_\ell^m(z)``
are the [associated Legendre
polynomials](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Definition_for_non-negative_integer_parameters_%E2%84%93_and_m).
The Condon-Shortley phase ``(-)^m`` is included in the definition of the Legendre
polynomials, consistent with Varshalovich (1988).

An explicit expression for the Legendre polynomials is given by the [Rodrigues formula](https://en.wikipedia.org/wiki/Rodrigues%27_formula):

```math
P_{\ell}^{m}(x) = \frac{(-)^m}{2^{\ell} \ell !} (1 - x^2)^{m/2} \frac{\mathrm{d}^{\ell+m}}{\mathrm{d}x^{\ell+m}} (x^2 - 1)^\ell
```

**Properties.**
The complex conjugate of a spherical harmonic can be expressed in terms of spherical harmonics:
```math
\bar{Y}^{\ell}_{m}(\theta,\varphi) = (-)^m Y^{\ell}_{-m}(\theta,\varphi)
```

The spherical harmonics are normalized

```math
\int_{0}^{2\pi} \int_{0}^{\pi}
\bar{Y}^{\ell_1}_{m_1}(\theta,\varphi)
Y^{\ell_2}_{m_2}(\theta,\varphi)
\sin\theta \diff{\theta} \diff{\varphi}
= \delta_{\ell_1 \ell_2} \delta_{m_1 m_2}
\tag{V5.1.6}
```

and the integral of three spherical harmonics is given by

```math
\int_{0}^{2\pi} \int_{0}^{\pi}
Y^{\ell_1}_{m_1}(\theta,\varphi)
Y^{\ell_2}_{m_2}(\theta,\varphi)
Y^{\ell_3}_{m_3}(\theta,\varphi)
\sin\theta \diff{\theta} \diff{\varphi}
\\= \frac{1}{\sqrt{4\pi}} \angroot{\ell_1,\ell_2,\ell_3}
\begin{pmatrix}
\ell_1 & \ell_2 & \ell_3 \\
0      & 0      & 0
\end{pmatrix}
\begin{pmatrix}
\ell_1 & \ell_2 & \ell_3 \\
m_1    & m_2    & m_3
\end{pmatrix}
\tag{V5.9.5}
```

!!! note "Differences from ISO 80000-2:2009"
    The [ISO 80000-2:2009 standard](https://en.wikipedia.org/wiki/ISO_80000-2) standardizes
    some mathematical notations and conventions for definitions of some special functions.

    **Index placement convention**
    Unlike the ISO standard, we put the ``\ell`` index on top and ``m`` on the bottom, to be
    consistent with the way the ``k`` and ``q`` indices are normally written for tensor
    operators.

    **Spherical harmonics with negative ``m``.**
    The Condon-Shortley phase in the Legendre polynomials is consistent with the ISO standard.
    However, the definition of spherical harmonics differs slightly. Namely, the standard defines
    the spherical harmonics as follows

    ```math
    Y_{m}^{\ell}(\theta,\varphi) =
    \sqrt{\frac{2\ell+1}{4\pi}\frac{(\ell-|m|)!}{(\ell+|m|)!}}
    P_{\ell}^{|m|}(\cos\theta)
    \mathrm{e}^{\im m \varphi}
    \tag{ISO19.17}
    ```

    which leads to the following relationship for the complex conjugate of ``Y_{m}^{\ell}``

    ```math
    \bar{Y}_{m}^{\ell}(\theta,\varphi) = Y_{-m}^{\ell}(\theta,\varphi)
    ```

    We opt to follow the (arguably, more common) non-ISO definition to stay consistent with the
    primary reference of Varshalovich (1988).

## Clebsch–Gordan coefficients

The Clebsch–Gordan coefficients are related to the 3j symbols as

```math
C_{j_1m_1j_2m_2}^{j_3m_3} \equiv
\braket{j_1m_1j_2m_2}{j_3m_3} =
(-)^{j_1-j_2+m_3}\angroot{j_3}
\begin{pmatrix}
j_1&j_2&j_3\\
m_1&m_2&-m_3
\end{pmatrix}.
\tag{V8.1.12}
```

They can be calculated with the [`clebsch_gordan`](@ref) function.

## Wigner–Eckart theorem

For the Wigner–Eckart theorem, which defines the **reduced matrix elements** (RMEs)
``\redmatrixel{n' j'}{\tensor{T}^{(k)}}{n j}``
of a tensor operator of rank ``k``, the convention is the following

```math
\begin{align}
\matrixel{n' j' m'}{\tensor{T}^{(k)}_q}{n j m}
&\defd
(-)^{2k} \frac{1}{\angroot{j'}}
C_{jm;kq}^{j'm'}
\redmatrixel{n' j'}{\tensor{T}^{(k)}}{n j} \nonumber \\
&=
(-)^{j'-m'}
\begin{pmatrix}
 j' & k & j \\
-m' & q & m
\end{pmatrix}
\redmatrixel{n' j'}{\tensor{T}^{(k)}}{n j}
\tag{V13.1.2}
\end{align}
```

The second form can be derived by using the relationship between the Clebsch–Gordan coefficients and the Wigner 3j symbols, and the permutation symmetries of the 3j symbol.
The ``n`` and ``n'`` labels represent all non angular momentum quantum numbers.

!!! note "Other conventions for RMEs"
    A simpler convention used by some books, that also generalizes to other symmetry groups, is
    ```math
    \matrixel{n' j' m'}{\tensor{T}^{(k)}_q}{n j m} =
    \braket{j_1 m_1 j_2 m_2}{j_3 m_3}
    \redmatrixel{n' j'}{\tensor{T}^{(k)}}{n j}
    ```
    However, again to stay consistent with Varshalovich (1988), we shall not use it.
    But it must be noted that, as the Wigner–Eckart theorem functions as a definition
    for the reduced matrix elements, this choice will change the values of the RMEs.

```@meta
CurrentModule = nothing
```
