# [Gradients](@id gradients)

The gradient is a rank-1 tensor, the reduced matrix element of which
is given by

```math
\begin{equation}
\label{eqn:gradient-rme}
\tag{V13.2.22}
\redmatrixel{n'\ell'}{\tensor{\nabla}^{(1)}}{n\ell} =
\sqrt{\ell+1}
A_{n'\ell'n\ell}\delta_{\ell',\ell+1} -
\sqrt{\ell}
B_{n'\ell'n\ell}\delta_{\ell',\ell-1},
\end{equation}
```

where

```math
\begin{equation}
\begin{aligned}
A_{n'\ell'n\ell} &\defd
\int_0^\infty\diff{r}r^2
\conj{\Psi}_{n'\ell'}(r)
\left(\partial_r-\frac{\ell}{r}\right)
\Psi_{n\ell}(r),\\
B_{n'\ell'n\ell} &\defd
\int_0^\infty\diff{r}r^2
\conj{\Psi}_{n'\ell'}(r)
\left(\partial_r+\frac{\ell+1}{r}\right)
\Psi_{n\ell}(r).
\end{aligned}
\tag{V13.2.23}
\label{eqn:gradient-rme-terms}
\end{equation}
```

!!! note
    When working in spherical coordinates, it is common to use the
    _reduced wavefunction_, since that simplifies the Laplacian and
    enforces vanishing Dirichlet boundary conditions at ``r=0``:

    ```math
    \Phi(\vec{r}) \defd r\Psi(\vec{r}).
    ```

    This must be taken into account when computing matrix elements of
    differential operators with respect to the partial waves of the
    reduced wavefunction:

    ```math
    \nabla \Psi = \nabla \frac{\Phi}{r} =
    \left(\nabla\frac{1}{r}\right)\Phi + \frac{1}{r}(\nabla\Phi) =
    -\frac{\tensor{n}_1}{r^2}\Phi + \frac{1}{r}\nabla\Phi =
    \frac{1}{r}\left(\nabla - \frac{\tensor{n}_1}{r}\right)\Phi,
    ```

    which means that if ``\ket{n'\ell'}`` and ``\ket{n\ell}`` are
    partial waves of ``\Phi``, instead of
    ``\eqref{eqn:gradient-rme}``, we need to use

    ```math
    \begin{equation}
    \tag{\ref{eqn:gradient-rme}*}
    \redmatrixel{n'\ell'}{\tensor{\nabla}^{(1)}-\frac{\tensor{n}_1}{r}}{n\ell}
    =
    \sqrt{\ell+1}
    \tilde{A}_{n'\ell'n\ell}\delta_{\ell',\ell+1} -
    \sqrt{\ell}
    \tilde{B}_{n'\ell'n\ell}\delta_{\ell',\ell-1},
    \end{equation}
    ```

    where

    ```math
    \begin{equation}
    \begin{aligned}
    \tilde{A}_{n'\ell'n\ell}\delta_{\ell',\ell+1}
    &\defd
    \int_0^\infty\diff{r}r^2
    \conj{\Psi}_{n'\ell'}(r)
    \left(\partial_r-\frac{\ell+1}{r}\right)
    \Psi_{n\ell}(r),\\
    \tilde{B}_{n'\ell'n\ell}\delta_{\ell',\ell-1}
    &\defd
    \int_0^\infty\diff{r}r^2
    \conj{\Psi}_{n'\ell'}(r)
    \left(\partial_r+\frac{\ell}{r}\right)
    \Psi_{n\ell}(r),
    \end{aligned}
    \tag{\ref{eqn:gradient-rme-terms}*}
    \end{equation}
    ```

    since

    ```math
    \begin{equation}
    \tag{V13.2.11}
    \redmatrixel{\ell'}{\tensor{n}_1}{\ell} =
    \angroot{\ell}
    C_{\ell 0;10}^{\ell'0}
    \end{equation}
    ```

    (cf [Spherical tensors](@ref tensors_spherical_tensors)
    since ``\tensor{n}_1\equiv\tensor{C}^{(1)}``) and

    ```math
    \begin{equation}
    \tag{V8.5.33,34}
    C_{\ell 0;10}^{\ell'0} =
    % (\ell+1)
    % \sqrt{\frac{(2\ell)!2!}{(2\ell+2)!}}
    % \delta_{\ell',\ell+1}
    % - \ell
    % \sqrt{\frac{2!(2\ell-1)!}{(2\ell+1)!}}
    % \delta_{\ell',\ell-1}
    \sqrt{\frac{(\ell+1)}{(2\ell+1)}}
    \delta_{\ell',\ell+1}
    - \sqrt{\frac{\ell}{(2\ell+1)}}
    \delta_{\ell',\ell-1}
    \end{equation}
    ```

    ```math
    \begin{equation}
    \tag{V13.2.11*}
    \implies\redmatrixel{\ell'}{\tensor{n}_1}{\ell} =
    \sqrt{\ell+1}
    \delta_{\ell',\ell+1} -
    \sqrt{\ell}
    \delta_{\ell',\ell-1},
    \end{equation}
    ```

    i.e. ``r^{-1}`` is subtracted from the centrifugal terms of ``A``
    and ``B`` in ``\eqref{eqn:gradient-rme-terms}``.

## Example

```jldoctest
julia> using AngularMomentumAlgebra, AtomicLevels

julia> orbitals = sos"k[s-d]"
18-element Array{SpinOrbital{Orbital{Symbol},Tuple{Int64,HalfIntegers.Half{Int64}}},1}:
 ks₀α
 ks₀β
 kp₋₁α
 kp₋₁β
 kp₀α
 kp₀β
 kp₁α
 kp₁β
 kd₋₂α
 kd₋₂β
 kd₋₁α
 kd₋₁β
 kd₀α
 kd₀β
 kd₁α
 kd₁β
 kd₂α
 kd₂β

julia> a,b,c = orbitals[[3,9,13]]
3-element Array{SpinOrbital{Orbital{Symbol},Tuple{Int64,HalfIntegers.Half{Int64}}},1}:
 kp₋₁α
 kd₋₂α
 kd₀α

julia> ∂x = cartesian_tensor_component(Gradient(), :x)
- 0.707107 𝛁̂⁽¹⁾₁ + 0.707107 𝛁̂⁽¹⁾₋₁

julia> dot(a, ∂x, b)
0.447214(∂ᵣ + 3/r)

julia> dot(b, ∂x, a)
0.447214(∂ᵣ - 1/r)

julia> dot(a, ∂x, c)
- 0.182574(∂ᵣ + 3/r)

julia> dot(c, ∂x, a)
- 0.182574(∂ᵣ - 1/r)
```

## Reference

```@docs
Gradient
system(::Type{Gradient})
AngularMomentumAlgebra.RadialGradientMatrixElement
rme((n′,ℓ′)::Tuple{<:Number, <:Number}, ::Gradient, (n,ℓ)::Tuple{<:Number, <:Number})
AngularMomentumAlgebra.couplings(tensor::Gradient, (n,ℓ)::Tuple{<:Number, <:Number})
```
