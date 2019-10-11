var documenterSearchIndex = {"docs":
[{"location":"tensors/#Tensors-1","page":"Tensors","title":"Tensors","text":"","category":"section"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"DocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\n    using LinearAlgebra\nend","category":"page"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"Tensor\nTensorComponent\nAngularMomentumAlgebra.LinearCombinationTensor","category":"page"},{"location":"tensors/#AngularMomentumAlgebra.Tensor","page":"Tensors","title":"AngularMomentumAlgebra.Tensor","text":"Tensor{k,label}\n\nAbstract base for any tensor of rank k.\n\n\n\n\n\n","category":"type"},{"location":"tensors/#AngularMomentumAlgebra.TensorComponent","page":"Tensors","title":"AngularMomentumAlgebra.TensorComponent","text":"TensorComponent(tensor, q)\n\nRepresents the qth component of a tensor; abs(q) ≤ rank(tensor).\n\n\n\n\n\n","category":"type"},{"location":"tensors/#AngularMomentumAlgebra.LinearCombinationTensor","page":"Tensors","title":"AngularMomentumAlgebra.LinearCombinationTensor","text":"LinearCombinationTensor\n\nRepresents a linear combination of tensor components.\n\n\n\n\n\n","category":"type"},{"location":"tensors/#Product-tensors-1","page":"Tensors","title":"Product tensors","text":"","category":"section"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"TensorProduct","category":"page"},{"location":"tensors/#AngularMomentumAlgebra.TensorProduct","page":"Tensors","title":"AngularMomentumAlgebra.TensorProduct","text":"TensorProduct{K}(T, U)\n\nA tensor of rank K formed from the product of the tensors T and U, according to\n\nbeginequation\ntagV3120\ntensorX^(K)_Q equiv\ntensorT^(k_1)tensorU^(k_2)^(K)_Q defd\ntensorT^(k_1)_q_1\ntensorU^(k_2)_q_2\nC_k_1q_1k_2q_2^KQ\nendequation\n\n\n\n\n\n","category":"type"},{"location":"tensors/#Scalar-product-1","page":"Tensors","title":"Scalar product","text":"","category":"section"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"The matrix element of scalar product of two tensors acting on different coordinates is given by","category":"page"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"beginequation\nbeginaligned\nmatrixeln_aj_am_an_bj_bm_btensorP^(k)(1)cdottensorQ^(k)(2)n_cj_cm_cn_dj_dm_d\n=\nfrac1angrootj_aj_b\nsum_alpha(-)^-alpha\nC_j_cm_ckalpha^j_am_a\nC_j_dm_dk-alpha^j_bm_b\ntimes\nredmatrixeln_aj_atensorP^(k)(1)n_cj_c\nredmatrixeln_bj_btensorQ^(k)(2)n_dj_d\nendaligned\ntagV13126\nendequation","category":"page"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"Since the Clebsch–Gordan coefficients can be rewritten using 3j symbols and the 3j symbols vanish unless m_c + alpha - m_3 = m_d - alpha - m_b = 0, we have","category":"page"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"beginequation\nalpha = m_a - m_c = m_d-m_b\nimplies\n-alpha + m_a + m_b = m_b + m_c\nendequation","category":"page"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"TensorScalarProduct\ndot(::Tensor, ::Tensor)\nAngularMomentumAlgebra.integrate_spinors((a,b), X::TensorScalarProduct, (c,d))","category":"page"},{"location":"tensors/#AngularMomentumAlgebra.TensorScalarProduct","page":"Tensors","title":"AngularMomentumAlgebra.TensorScalarProduct","text":"TensorScalarProduct(T, U)\n\nA tensor of rank 0 formed from the product of the tensors T and U (which have to have the same rank), according to\n\nbeginequation\ntagV313035\n(tensorT^(k) cdot\ntensorU^(k)) defd\n(-)^kangrootk\ntensorT^(k)\ntensorU^(k)^(0)_0 equiv\n(-)^q\ntensorT^(k)_q\ntensorU^(k)_-q\nendequation\n\n\n\n\n\n","category":"type"},{"location":"tensors/#LinearAlgebra.dot-Tuple{Tensor,Tensor}","page":"Tensors","title":"LinearAlgebra.dot","text":"dot(T::Tensor, U::Tensor)\n\nForm the scalar product of the two tensors T and U, which need to have the same rank.\n\nExamples\n\njulia> SphericalTensor(4)⋅SphericalTensor(4)\n(𝐂⁽⁴⁾⋅𝐂⁽⁴⁾)\n\n\n\n\n\n","category":"method"},{"location":"tensors/#AngularMomentumAlgebra.integrate_spinors-Tuple{Any,TensorScalarProduct,Any}","page":"Tensors","title":"AngularMomentumAlgebra.integrate_spinors","text":"integrate_spinors((a,b), X, (c,d))\n\nPerform the spin-angular integration of the scalar-product tensor X≡(T⁽ᵏ⁾⋅U⁽ᵏ⁾), where T acts on the coordinate of orbitals a & c and similarly, U acts on the coordinate of orbitals b & d, according to Eq. (13.1.26) of Varshalovich (1988).\n\n\n\n\n\n","category":"method"},{"location":"tensors/#The-Wigner–Eckart-theorem-1","page":"Tensors","title":"The Wigner–Eckart theorem","text":"","category":"section"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"The Wigner–Eckart theorem states that the matrix element of a tensor component tensorT^(k)_q can be evaluated as","category":"page"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"beginequation\nmatrixelnjmtensorT^(k)_qnjm=\n(-)^j-m\nbeginpmatrix\njkj\n-mqm\nendpmatrix\nredmatrixelnjtensorT^(k)nj\ntagV1312\nendequation","category":"page"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"where the reduced matrix element redmatrixelnjtensorT^(k)nj does not depend on mm.","category":"page"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"wigner_eckart","category":"page"},{"location":"tensors/#AngularMomentumAlgebra.wigner_eckart","page":"Tensors","title":"AngularMomentumAlgebra.wigner_eckart","text":"wigner_eckart(j′, m′, T⁽ᵏ⁾q, j, m)\n\nComputes the (spin-angular part of the) matrix element ⟨n′j′m′|Tᵏq|njm⟩, where T⁽ᵏ⁾q is the qth component of a tensor of rank k, using the definition of Eq. (13.1.2) in Varshalovich (1988).\n\n\n\n\n\nwigner_eckart(o′, T⁽ᵏ⁾q, o)\n\nComputes the (spin-angular part of the) matrix element ⟨o′|Tᵏq|o⟩, where T⁽ᵏ⁾q is the qth component of a tensor of rank k, using the definition of Eq. (13.1.2) in Varshalovich (1988).\n\n\n\n\n\n","category":"function"},{"location":"tensors/#","page":"Tensors","title":"Tensors","text":"CurrentModule = nothing\nDocTestSetup = nothing","category":"page"},{"location":"common/#Common-routines-1","page":"Common routines","title":"Common routines","text":"","category":"section"},{"location":"common/#","page":"Common routines","title":"Common routines","text":"CurrentModule = AngularMomentumAlgebra\nDocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\nend","category":"page"},{"location":"common/#","page":"Common routines","title":"Common routines","text":"∏\ntriangle_range\npowneg1\njmⱼ\nspin\nLinearCombination\n@linearly_combinable\nclebsch_gordan_condon_shortley","category":"page"},{"location":"common/#AngularMomentumAlgebra.∏","page":"Common routines","title":"AngularMomentumAlgebra.∏","text":"∏(ℓs...)\n\nCalculates √((2ℓ₁+1)(2ℓ₂+1)...(2ℓₙ+1)), which is a common factor in angular momentum algebra.\n\n\n\n\n\n","category":"function"},{"location":"common/#AngularMomentumAlgebra.triangle_range","page":"Common routines","title":"AngularMomentumAlgebra.triangle_range","text":"triangle_range(a,b)\n\nFind all (even) k such that |a-b| ≤ k ≤ a + b. This is useful when expanding matrix elements of tensors between angular momenta a and b in multipoles k; triangle_range can then be used to decided which multipole terms are required.\n\n\n\n\n\n","category":"function"},{"location":"common/#AngularMomentumAlgebra.powneg1","page":"Common routines","title":"AngularMomentumAlgebra.powneg1","text":"powneg1(k) = (-)ᵏ\n\nCalculates powers of negative unity for integer k.\n\n\n\n\n\n","category":"function"},{"location":"common/#AngularMomentumAlgebra.jmⱼ","page":"Common routines","title":"AngularMomentumAlgebra.jmⱼ","text":"jmⱼ(o::SpinOrbital)\n\nReturn the angular momentum and its projection on the z axis of the spin-orbital o.\n\n\n\n\n\n","category":"function"},{"location":"common/#AngularMomentumAlgebra.spin","page":"Common routines","title":"AngularMomentumAlgebra.spin","text":"spin(o::SpinOrbital)\n\nReturn the spin of the spin-orbital o.\n\n\n\n\n\n","category":"function"},{"location":"common/#AngularMomentumAlgebra.LinearCombination","page":"Common routines","title":"AngularMomentumAlgebra.LinearCombination","text":"LinearCombination(Ts::Vector{T}, coeffs::Vector{N})\n\nRepresents a general linear combination of objects of type T.\n\n\n\n\n\n","category":"type"},{"location":"common/#AngularMomentumAlgebra.@linearly_combinable","page":"Common routines","title":"AngularMomentumAlgebra.@linearly_combinable","text":"@linearly_combinable TT\n\nTurns the type TT into a linearly combinable type, i.e. defines arithmetic operators.\n\nExamples\n\njulia> @linearly_combinable Symbol\n\njulia> 4*(:x) - 5*(:y)\n4 :x - 5 :y\n\n\n\n\n\n","category":"macro"},{"location":"common/#AngularMomentumAlgebra.clebsch_gordan_condon_shortley","page":"Common routines","title":"AngularMomentumAlgebra.clebsch_gordan_condon_shortley","text":"clebsch_gordan_condon_shortley(j₁, m₁, j₂, m₂, j₃, m₃=m₁+m₂)\n\nCalculate the vector coupling coefficient ⟨j₁ m₁, j₂ m₂|j₃ m₃⟩ according to the Condon–Shortley phase convention and using the definition of Eq. (8.1.12) in Varshalovich (1988).\n\n\n\n\n\n","category":"function"},{"location":"common/#","page":"Common routines","title":"Common routines","text":"CurrentModule = nothing\nDocTestSetup = nothing","category":"page"},{"location":"spherical_tensors/#tensors_spherical_tensors-1","page":"Spherical Tensors","title":"Spherical tensors","text":"","category":"section"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"DocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\n    using LinearAlgebra\nend","category":"page"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"The spherical tensors are related to the spherical harmonics as","category":"page"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"tensorC^(k)_q defd\nsqrtfrac4pi2k+1\nY^k_q\ntagV517","category":"page"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"SphericalTensor","category":"page"},{"location":"spherical_tensors/#AngularMomentumAlgebra.SphericalTensor","page":"Spherical Tensors","title":"AngularMomentumAlgebra.SphericalTensor","text":"SphericalTensor(k)\n\nConstruct a spherical tensor of rank k.\n\n\n\n\n\n","category":"type"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"The reduced matrix element of the spherical tensor is given by","category":"page"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"beginaligned\nredmatrixelelltensorC^(k)ell\n=\nangrootell\nC_ell 0k0^ell0 =\n(-)^ell-k\nangrootellell\nwignerthreejellkell000\nendaligned\ntagV132107","category":"page"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"For the reduced matrix element of the same spherical tensor in the ell s j basis, we must use the uncoupling formula together with the above formula in the uncoupled basis:","category":"page"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"beginaligned\nredmatrixelell s jtensorC^(k)ell s j\n=\ndelta_ss(-)^j+ell+s+k\nangrootjj\nwignersixjellsjjkell\nredmatrixelelltensorC^(k)ell\n=\n(-)^ell+ell+s+j\nangrootellelljj\nwignersixjellsjjkell\nwignerthreejellkell000\nendaligned\ntagV1314025","category":"page"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"rme\nAngularMomentumAlgebra.ranks","category":"page"},{"location":"spherical_tensors/#AngularMomentumAlgebra.rme","page":"Spherical Tensors","title":"AngularMomentumAlgebra.rme","text":"rme(ℓ′,Cᵏ,ℓ)\n\nCalculate the reduced matrix element ⟨ℓ′||C⁽ᵏ⁾||ℓ⟩ of the spherical tensor of rank k. Condon–Shortley phase convention and using the definition of Eq. (13.2.107) in Varshalovich (1988).\n\n\n\n\n\nrme(o′, Cᵏ, o)\n\nCalculate the reduced matrix element ⟨o′||C⁽ᵏ⁾||o⟩ of the spherical tensor of rank k, between the relativistic spin-orbitals o′ and o. Since the spherical tensors only act on the quantum numbers ℓ′ and ℓ, the reduced matrix has to be evaluated via the uncoupling formula given in Eqs. (13.1.40) & (13.2.5) in Varshalovich (1988), together with reduced matrix element for the uncoupled angular momenta ℓ′ and ℓ given by Eq. (13.2.107) (ibid).\n\n\n\n\n\n","category":"function"},{"location":"spherical_tensors/#AngularMomentumAlgebra.ranks","page":"Spherical Tensors","title":"AngularMomentumAlgebra.ranks","text":"ranks(a, ::Type{SphericalTensor}, b)\n\nReturn which tensor ranks for spherical tensors that fulfill the triangle condition between spin-orbitals a and b.\n\n\n\n\n\n","category":"function"},{"location":"spherical_tensors/#Dipole-operator-1","page":"Spherical Tensors","title":"Dipole operator","text":"","category":"section"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"The dipole operator is a rank-1 Cartesian tensor that may be expressed using the rank-1 spherical tensor:","category":"page"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"hatvecr equiv\nbeginbmatrixhatxhatyhatzendbmatrix\nequiv\nbeginbmatrix\nfrac1sqrt2-tensorC^(1)_1 + tensorC^(1)_-1\nfracimsqrt2tensorC^(1)_1 + tensorC^(1)_-1\ntensorC^(1)_0\nendbmatrix","category":"page"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"CurrentModule = AngularMomentumAlgebra","category":"page"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"Dipoles.𝐫̂","category":"page"},{"location":"spherical_tensors/#AngularMomentumAlgebra.Dipoles.𝐫̂","page":"Spherical Tensors","title":"AngularMomentumAlgebra.Dipoles.𝐫̂","text":"𝐫̂\n\nThe angular part of the dipole operator; the elements correspond to [x,y,z], i.e. the Cartesian tensor components. Can be entered as \\bfr\\hat.\n\nExamples\n\njulia> using AngularMomentumAlgebra.Dipoles\n\njulia> z = 𝐫̂[3]\n𝐂⁽¹⁾₀\n\njulia> wigner_eckart(0, 0, z, 1, 0)\n0.5773502691896256\n\n\n\n\n\n","category":"constant"},{"location":"spherical_tensors/#","page":"Spherical Tensors","title":"Spherical Tensors","text":"CurrentModule = nothing\nDocTestSetup = nothing","category":"page"},{"location":"definitions/#Definitions-1","page":"Definitions","title":"Definitions","text":"","category":"section"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"CurrentModule = AngularMomentumAlgebra","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"This page defines much of the general notation and conventions used in the code. Where possible, we are consistent with Varshalovich (1988).","category":"page"},{"location":"definitions/#Shorthands-1","page":"Definitions","title":"Shorthands","text":"","category":"section"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"The abbreviation (-)^kdefd(-1)^k is used for the roots of negative unity (see also powneg1).","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"A commonly occurring factor in angular momentum algebra is","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"angrootj_1j_2j_n\ndefd(2j_1+1)(2j_2+1)(2j_n+1)^12\ntagV1313½","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"It can be calculated with the unexported AngularMomentumAlgebra.∏ function.","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"Indices appearing in pairs on only one side of an equation are implicitly summed over.","category":"page"},{"location":"definitions/#Spherical-harmonics-1","page":"Definitions","title":"Spherical harmonics","text":"","category":"section"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"We assume the following definition of the spherical harmonics","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"Y_m^ell(thetavarphi) = sqrtfrac2ell+14pifrac(ell-m)(ell+m) P_ell^m(costheta) mathrme^im m varphi\ntagV521","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"where theta and varphi are the usual spherical coordinates and P_ell^m(z) are the associated Legendre polynomials. The Condon-Shortley phase (-)^m is included in the definition of the Legendre polynomials, consistent with Varshalovich (1988).","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"An explicit expression for the Legendre polynomials is given by the Rodrigues formula:","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"P_ell^m(x) = frac(-)^m2^ell ell  (1 - x^2)^m2 fracmathrmd^ell+mmathrmdx^ell+m (x^2 - 1)^ell","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"Properties. The complex conjugate of a spherical harmonic can be expressed in terms of spherical harmonics:","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"barY^ell_m(thetavarphi) = (-)^m Y^ell_-m(thetavarphi)","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"The spherical harmonics are normalized","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"int_0^2pi int_0^pi\nbarY^ell_1_m_1(thetavarphi)\nY^ell_2_m_2(thetavarphi)\nsintheta difftheta diffvarphi\n= delta_ell_1 ell_2 delta_m_1 m_2\ntagV516","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"and the integral of three spherical harmonics is given by","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"int_0^2pi int_0^pi\nY^ell_1_m_1(thetavarphi)\nY^ell_2_m_2(thetavarphi)\nY^ell_3_m_3(thetavarphi)\nsintheta difftheta diffvarphi\n= frac1sqrt4pi angrootell_1ell_2ell_3\nbeginpmatrix\nell_1  ell_2  ell_3 \n0       0       0\nendpmatrix\nbeginpmatrix\nell_1  ell_2  ell_3 \nm_1     m_2     m_3\nendpmatrix\ntagV595","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"note: Differences from ISO 80000-2:2009\nThe ISO 80000-2:2009 standard standardizes some mathematical notations and conventions for definitions of some special functions.Index placement convention Unlike the ISO standard, we put the ell index on top and m on the bottom, to be consistent with the way the k and q indices are normally written for tensor operators.Spherical harmonics with negative m. The Condon-Shortley phase in the Legendre polynomials is consistent with the ISO standard. However, the definition of spherical harmonics differs slightly. Namely, the standard defines the spherical harmonics as followsY_m^ell(thetavarphi) = sqrtfrac2ell+14pifrac(ell-m)(ell+m) P_ell^m(costheta) mathrme^i m varphi\ntagISO1917which leads to the following relationship for the complex conjugate of Y_m^ellbarY_m^ell(thetavarphi) = Y_-m^ell(thetavarphi)We opt to follow the (arguably, more common) non-ISO definition to stay consistent with the primary reference of Varshalovich (1988).","category":"page"},{"location":"definitions/#Clebsch–Gordan-coefficients-1","page":"Definitions","title":"Clebsch–Gordan coefficients","text":"","category":"section"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"The Clebsch–Gordan coefficients are related to the 3j symbols as","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"C_j_1m_1j_2m_2^j_3m_3 equiv\nbraketj_1m_1j_2m_2j_3m_3 =\n(-)^j_1-j_2+m_3angrootj_3\nbeginpmatrix\nj_1j_2j_3\nm_1m_2-m_3\nendpmatrix\ntagV8112","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"They can be calculated with the clebsch_gordan_condon_shortley function.","category":"page"},{"location":"definitions/#Wigner–Eckart-theorem-1","page":"Definitions","title":"Wigner–Eckart theorem","text":"","category":"section"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"For the Wigner–Eckart theorem, which defines the reduced matrix elements (RMEs) redmatrixeln jtensorT^(k)n j of a tensor operator of rank k, the convention is the following","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"beginalign\nmatrixeln j mtensorT^(k)_qn j m\ndefd\n(-)^2k fracbraketjmkqj mangrootj\nredmatrixeln jtensorT^(k)n j nonumber \n=\n(-)^j-m\nbeginpmatrix\n j  k  j \n-m  q  m\nendpmatrix\nredmatrixeln jtensorT^(k)n j\ntagV1312\nendalign","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"The second form can be derived by using the relationship between the Clebsch–Gordan coefficients and the Wigner 3j symbols, and the permutation symmetries of the 3j symbol. The n and n labels represent all non angular momentum quantum numbers.","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"note: Other conventions for RMEs\nA simpler convention used by some books, that also generalizes to other symmetry groups, ismatrixeln j mtensorT^(k)_qn j m =\nbraketj_1 m_1 j_2 m_2j_3 m_3\nredmatrixeln jtensorT^(k)n jHowever, again to stay consistent with Varshalovich (1988), we shall not use it. But it must be noted that, as the Wigner–Eckart theorem functions as a definition for the reduced matrix elements, this choice will change the values of the RMEs.","category":"page"},{"location":"definitions/#","page":"Definitions","title":"Definitions","text":"CurrentModule = nothing","category":"page"},{"location":"energy_expressions/#Energy-expressions-1","page":"Energy expressions","title":"Energy expressions","text":"","category":"section"},{"location":"energy_expressions/#","page":"Energy expressions","title":"Energy expressions","text":"DocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\nend","category":"page"},{"location":"energy_expressions/#","page":"Energy expressions","title":"Energy expressions","text":"AngularMomentumAlgebra.integrate_spinor\nMatrix","category":"page"},{"location":"energy_expressions/#AngularMomentumAlgebra.integrate_spinor","page":"Energy expressions","title":"AngularMomentumAlgebra.integrate_spinor","text":"integrate_spinor(me)\n\nPerform the spin-angular integration of the matrix element me, leaving only a radial integral multiplied by a spin-angular coefficient. The spin-angular integral is dependent on the specific combination of spin-orbitals and the operator (expressed as a tensor); the default implementation is to leave me as-is, corresponding to a spin-angular integral of unity.\n\n\n\n\n\nintegrate_spinor(me::OrbitalMatrixElement{2,<:Any,<:CoulombInteraction,<:Any})\n\nPerform the spin-angular integration of the two-electron matrix element me, by first multipole-expanding the Coulomb interaction and then integrating all the resulting terms over the spin-angular coordinates (see multipole_expand).\n\n\n\n\n\nintegrate_spinor(integral::NBodyTermFactor)\n\nDummy method that returns integral unchanged, used for all NBodyTermFactors that are not to be multipole-expanded.\n\n\n\n\n\n","category":"function"},{"location":"energy_expressions/#Base.Matrix","page":"Energy expressions","title":"Base.Matrix","text":"Matrix(a::Vector{<:Configuration{<:SpinOrbital}},\n       op::QuantumOperator,\n       b::Vector{<:Configuration{<:SpinOrbital}}[, overlaps])\n\nGenerate the energy-expression associated with the quantum operator op, in the basis of the spin-configurations a and b, with an optional set of orbital overlaps, specifying any desired non-orthogonalities. The energy expression is generated in a basis-agnostic way by EnergyExpressions.jl and each term is then integrated over the spin-angular coordinates using integrate_spinor.\n\n\n\n\n\nMatrix(op::QuantumOperator,\n       spin_cfgs::Vector{<:Configuration{<:SpinOrbital}}[, overlaps])\n\nGenerate the energy-expression associated with the quantum operator op, in the basis of the spin-configurations spin_cfgs, with an optional set of orbital overlaps, specifying any desired non-orthogonalities. The energy expression is generated in a basis-agnostic way by EnergyExpressions.jl and each term is then integrated over the spin-angular coordinates using integrate_spinor.\n\n\n\n\n\n","category":"type"},{"location":"energy_expressions/#","page":"Energy expressions","title":"Energy expressions","text":" DocTestSetup = nothing","category":"page"},{"location":"multipole_expansions/#Multipole-expansions-1","page":"Multipole expansions","title":"Multipole expansions","text":"","category":"section"},{"location":"multipole_expansions/#","page":"Multipole expansions","title":"Multipole expansions","text":"CurrentModule = AngularMomentumAlgebra\nDocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\nend","category":"page"},{"location":"multipole_expansions/#","page":"Multipole expansions","title":"Multipole expansions","text":"multipole_expand_scalar_product\nmultipole_expand","category":"page"},{"location":"multipole_expansions/#AngularMomentumAlgebra.multipole_expand_scalar_product","page":"Multipole expansions","title":"AngularMomentumAlgebra.multipole_expand_scalar_product","text":"multipole_expand_scalar_product(a, b, P, Q, c, d)\n\nMultipole-expand the matrix element ⟨ab|P⋅Q|cd⟩, where the tensor P acts on orbitals a & c, and the tensor Q acts on orbitals b & d. The definition is taken from Eq. (13.1.26) of Varshalovich (1988).\n\n\n\n\n\n","category":"function"},{"location":"multipole_expansions/#AngularMomentumAlgebra.multipole_expand","page":"Multipole expansions","title":"AngularMomentumAlgebra.multipole_expand","text":"multipole_expand(integral::OrbitalMatrixElement{2,A,<:CoulombInteraction,B})\n\nMultipole-expand the two-body integral resulting from the Coulomb repulsion between two electrons.\n\n\n\n\n\n","category":"function"},{"location":"multipole_expansions/#","page":"Multipole expansions","title":"Multipole expansions","text":" DocTestSetup = nothing","category":"page"},{"location":"#AngularMomentumAlgebra.jl-1","page":"Home","title":"AngularMomentumAlgebra.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Documentation for AngularMomentumAlgebra.jl","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Equation labelled (Vx.y.z) references equation (z) from section x.y of","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Varshalovich, D. A. (1988). Quantum Theory of Angular Momentum: Irreducible Tensors, Spherical Harmonics, Vector Coupling Coefficients, 3nj Symbols. Singapore Teaneck, NJ, USA: World Scientific Pub.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The Condon–Shortley convention is employed throughout (if you spot a mistake, please tell me!), and formulas are taken from Varshalovich (1988) only, to ensure consistency.","category":"page"},{"location":"coulomb/#Coulomb-interaction-1","page":"Coulomb interaction","title":"Coulomb interaction","text":"","category":"section"},{"location":"coulomb/#","page":"Coulomb interaction","title":"Coulomb interaction","text":"DocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\nend","category":"page"},{"location":"coulomb/#","page":"Coulomb interaction","title":"Coulomb interaction","text":"The Coulomb interaction between two electrons at coordinates 1 and 2, respectively, can in spherical coordinates be multipole-expanded as:","category":"page"},{"location":"coulomb/#","page":"Coulomb interaction","title":"Coulomb interaction","text":"beginequation\nfrac1r_12 =\nsum_k=0^infty\nfracr_^kr_^k+1\nP_k(costheta)\nendequation","category":"page"},{"location":"coulomb/#","page":"Coulomb interaction","title":"Coulomb interaction","text":"which by the addition theorem for the spherical harmonics can be further expanded as","category":"page"},{"location":"coulomb/#","page":"Coulomb interaction","title":"Coulomb interaction","text":"beginequation\nfrac1r_12 =\nsum_kq\nfrac4pi2k+1\nfracr_^kr_^k+1\nconjY^k_q(1)Y^k_q(2)equiv\nsum_kq\nfracr_^kr_^k+1\ntensorC^k_q(1)cdottensorC^k_q(2)\ntagV5179\nendequation","category":"page"},{"location":"coulomb/#","page":"Coulomb interaction","title":"Coulomb interaction","text":"where we in the last step have used the definition of the Spherical tensors.","category":"page"},{"location":"coulomb/#","page":"Coulomb interaction","title":"Coulomb interaction","text":"CoulombInteractionMultipole\nCoulombPotentialMultipole","category":"page"},{"location":"coulomb/#AngularMomentumAlgebra.CoulombInteractionMultipole","page":"Coulomb interaction","title":"AngularMomentumAlgebra.CoulombInteractionMultipole","text":"CoulombInteractionMultipole(k, g)\n\nRepresents the kth multipole of the multipole expansion of the Coulomb interaction g.\n\n\n\n\n\n","category":"type"},{"location":"coulomb/#AngularMomentumAlgebra.CoulombPotentialMultipole","page":"Coulomb interaction","title":"AngularMomentumAlgebra.CoulombPotentialMultipole","text":"CoulombPotentialMultipole\n\nType alias for contraction of the CoulombInteractionMultipole over one coordinate, thereby forming a potential in the other coordinate.\n\n\n\n\n\n","category":"type"},{"location":"coulomb/#","page":"Coulomb interaction","title":"Coulomb interaction","text":" DocTestSetup = nothing","category":"page"}]
}
