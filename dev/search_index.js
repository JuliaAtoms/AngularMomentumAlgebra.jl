var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#AngularMomentumAlgebra.jl-1",
    "page": "Home",
    "title": "AngularMomentumAlgebra.jl",
    "category": "section",
    "text": "Documentation for AngularMomentumAlgebra.jlEquations labelled (Vx.y.z) referenceVarshalovich, D. A. (1988). Quantum Theory of Angular Momentum: Irreducible Tensors, Spherical Harmonics, Vector Coupling Coefficients, 3nj Symbols. Singapore Teaneck, NJ, USA: World Scientific Pub.The Condon–Shortley convention is employed throughout (if you spot a mistake, please tell me!), and formulas are taken from Varshalovich (1988) only, to ensure consistency."
},

{
    "location": "common/#",
    "page": "Common routines",
    "title": "Common routines",
    "category": "page",
    "text": ""
},

{
    "location": "common/#AngularMomentumAlgebra.∏",
    "page": "Common routines",
    "title": "AngularMomentumAlgebra.∏",
    "category": "function",
    "text": "∏(ℓs...)\n\nCalculates √((2ℓ₁+1)(2ℓ₂+1)...(2ℓₙ+1)), which is a common factor in angular momentum algebra.\n\n\n\n\n\n"
},

{
    "location": "common/#AngularMomentumAlgebra.triangle_range",
    "page": "Common routines",
    "title": "AngularMomentumAlgebra.triangle_range",
    "category": "function",
    "text": "triangle_range(a,b)\n\nFind all (even) k such that |a-b| ≤ k ≤ a + b. This is useful when expanding matrix elements of tensors between angular momenta a and b in multipoles k; triangl_range can then be used to decided which multipole terms are required.\n\n\n\n\n\n"
},

{
    "location": "common/#AngularMomentumAlgebra.powneg1",
    "page": "Common routines",
    "title": "AngularMomentumAlgebra.powneg1",
    "category": "function",
    "text": "powneg1(k) = (-)ᵏ\n\nCalculates powers of negative unity for integer k.\n\n\n\n\n\n"
},

{
    "location": "common/#AngularMomentumAlgebra.jmⱼ",
    "page": "Common routines",
    "title": "AngularMomentumAlgebra.jmⱼ",
    "category": "function",
    "text": "jmⱼ(o::SpinOrbital)\n\nReturn the angular momentum and its projection on the z axis of the spin-orbital o.\n\n\n\n\n\n"
},

{
    "location": "common/#AngularMomentumAlgebra.spin",
    "page": "Common routines",
    "title": "AngularMomentumAlgebra.spin",
    "category": "function",
    "text": "spin(o::SpinOrbital)\n\nReturn the spin of the spin-orbital o.\n\n\n\n\n\n"
},

{
    "location": "common/#Common-routines-1",
    "page": "Common routines",
    "title": "Common routines",
    "category": "section",
    "text": "CurrentModule = AngularMomentumAlgebra\nDocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\nendA commonly occurring factor in angular momentum algebra isbeginequation\nangrootj_1j_2j_n\ndefd(2j_1+1)(2j_2+1)(2j_n+1)^12\ntagV1313½\nendequation∏\ntriangle_range\npowneg1\njmⱼ\nspin DocTestSetup = nothing"
},

{
    "location": "clebsch_gordan/#",
    "page": "Clebsch–Gordan coefficients",
    "title": "Clebsch–Gordan coefficients",
    "category": "page",
    "text": ""
},

{
    "location": "clebsch_gordan/#AngularMomentumAlgebra.clebsch_gordan_condon_shortley",
    "page": "Clebsch–Gordan coefficients",
    "title": "AngularMomentumAlgebra.clebsch_gordan_condon_shortley",
    "category": "function",
    "text": "clebsch_gordan_condon_shortley(j₁, m₁, j₂, m₂, j₃, m₃=m₁+m₂)\n\nCalculate the vector coupling coefficient ⟨j₁ m₁, j₂ m₂|j₃ m₃⟩ according to the Condon–Shortley phase convention and using the definition of Eq. (8.1.12) in Varshalovich (1988).\n\n\n\n\n\n"
},

{
    "location": "clebsch_gordan/#Clebsch–Gordan-coefficients-1",
    "page": "Clebsch–Gordan coefficients",
    "title": "Clebsch–Gordan coefficients",
    "category": "section",
    "text": "The Clebsch–Gordan coefficients are related to the 3j symbols asbeginequation\nC_j_1m_1j_2m_2^j_3m_3 equiv\nbraketj_1m_1j_2m_2j_3m_3 =\n(-)^j_1-j_2+m_3angrootj_3\nbeginpmatrix\nj_1j_2j_3\nm_1m_2-m_3\nendpmatrix\ntagV8112\nendequationDocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\nendclebsch_gordan_condon_shortley DocTestSetup = nothing"
},

{
    "location": "tensors/#",
    "page": "Tensors",
    "title": "Tensors",
    "category": "page",
    "text": ""
},

{
    "location": "tensors/#AngularMomentumAlgebra.Tensor",
    "page": "Tensors",
    "title": "AngularMomentumAlgebra.Tensor",
    "category": "type",
    "text": "Tensor\n\nAbstract base for any tensor.\n\n\n\n\n\n"
},

{
    "location": "tensors/#Tensors-1",
    "page": "Tensors",
    "title": "Tensors",
    "category": "section",
    "text": "DocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\nendTensor"
},

{
    "location": "tensors/#AngularMomentumAlgebra.SphericalTensor",
    "page": "Tensors",
    "title": "AngularMomentumAlgebra.SphericalTensor",
    "category": "type",
    "text": "SphericalTensor(k)\n\nConstruct a spherical tensor of rank k.\n\n\n\n\n\n"
},

{
    "location": "tensors/#AngularMomentumAlgebra.rme",
    "page": "Tensors",
    "title": "AngularMomentumAlgebra.rme",
    "category": "function",
    "text": "rme(ℓ′,Cᵏ,ℓ)\n\nCalculate the reduced matrix element ⟨ℓ′||C⁽ᵏ⁾||ℓ⟩ of the spherical tensor of rank k. Condon–Shortley phase convention and using the definition of Eq. (13.2.107) in Varshalovich (1988).\n\n\n\n\n\n"
},

{
    "location": "tensors/#tensors_spherical_tensors-1",
    "page": "Tensors",
    "title": "Spherical tensors",
    "category": "section",
    "text": "The spherical tensors are related to the spherical harmonics asbeginequation\ntensorC^(k)_q defd\nsqrtfrac4pi2k+1\nY^k_q\ntagV517\nendequationSphericalTensorThe reduced matrix element of the spherical tensor is given bybeginequation\nbeginaligned\nredmatrixelelltensorC^(k)ell\n=\nangrootell\nC_ell 0k0^ell0 =\n(-)^ell-k\nangrootellell\nbeginpmatrix\nellkell000\nendpmatrix\nendaligned\ntagV132107\nendequationrme DocTestSetup = nothing"
},

{
    "location": "coulomb/#",
    "page": "Coulomb interaction",
    "title": "Coulomb interaction",
    "category": "page",
    "text": ""
},

{
    "location": "coulomb/#AngularMomentumAlgebra.CoulombInteractionMultipole",
    "page": "Coulomb interaction",
    "title": "AngularMomentumAlgebra.CoulombInteractionMultipole",
    "category": "type",
    "text": "CoulombInteractionMultipole(k)\n\nRepresents the kth multipole of the multipole expansion of the Coulomb interaction.\n\n\n\n\n\n"
},

{
    "location": "coulomb/#AngularMomentumAlgebra.CoulombPotentialMultipole",
    "page": "Coulomb interaction",
    "title": "AngularMomentumAlgebra.CoulombPotentialMultipole",
    "category": "type",
    "text": "CoulombPotentialMultipole\n\nType alias for contraction of the CoulombInteractionMultipole over one coordinate, thereby forming a potential in the other coordinate.\n\n\n\n\n\n"
},

{
    "location": "coulomb/#Coulomb-interaction-1",
    "page": "Coulomb interaction",
    "title": "Coulomb interaction",
    "category": "section",
    "text": "DocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\nendThe Coulomb interaction between two electrons at coordinates 1 and 2, respectively, can in spherical coordinates be multipole-expanded as:beginequation\nfrac1r_12 = ()\nsum_k=0^infty\nfracr_^kr_^k+1\nP_k(costheta)\nendequationwhich by the addition theorem for the spherical harmonics can be further expanded asbeginequation\nfrac1r_12 =\nsum_kq\nfrac4pi2k+1\nfracr_^kr_^k+1\nconjY^k_q(1)Y^k_q(2)equiv\nsum_kq\nfracr_^kr_^k+1\ntensorC^k_q(1)cdottensorC^k_q(2)\ntagV5179\nendequationwhere we in the last step have used the definition of the Spherical tensors.CoulombInteractionMultipole\nCoulombPotentialMultipole DocTestSetup = nothing"
},

{
    "location": "multipole_expansions/#",
    "page": "Multipole expansions",
    "title": "Multipole expansions",
    "category": "page",
    "text": ""
},

{
    "location": "multipole_expansions/#AngularMomentumAlgebra.multipole_expand_scalar_product",
    "page": "Multipole expansions",
    "title": "AngularMomentumAlgebra.multipole_expand_scalar_product",
    "category": "function",
    "text": "multipole_expand_scalar_product(a, b, P, Q, c, d)\n\nMultipole expand the matrix element ⟨ab|P⋅Q|cd⟩, where the tensor P acts on orbitals a & c, and the tensor Q acts on orbitals b & d. The definition is taken from Eq. (13.1.26) of Varshalovich (1988).\n\n\n\n\n\n"
},

{
    "location": "multipole_expansions/#AngularMomentumAlgebra.multipole_expand",
    "page": "Multipole expansions",
    "title": "AngularMomentumAlgebra.multipole_expand",
    "category": "function",
    "text": "multipole_expand(integral::OrbitalMatrixElement{2,A,CoulombInteraction,B})\n\nMultipole-expand the two-body integral resulting from the Coulomb repulsion between two electrons.\n\n\n\n\n\nmultipole_expand(integral::NBodyTermFactor)\n\nDummy method that returns integral unchanged, used for all NBodyTermFactors that are not to be multipole-expanded.\n\n\n\n\n\n"
},

{
    "location": "multipole_expansions/#Multipole-expansions-1",
    "page": "Multipole expansions",
    "title": "Multipole expansions",
    "category": "section",
    "text": "CurrentModule = AngularMomentumAlgebra\nDocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\nendThe matrix element of scalar product of two tensors acting on different coordinates is given bybeginequation\nbeginaligned\nmatrixeln_aj_am_an_bj_bm_btensorP^(k)(1)cdottensorQ^(k)(2)n_cj_cm_cn_dj_dm_d\n=\nfrac1angrootj_aj_b\nsum_alpha(-)^-alpha\nC_j_cm_ckalpha^j_am_a\nC_j_dm_dk-alpha^j_bm_b\ntimes\nredmatrixeln_aj_atensorP^(k)(1)n_cj_c\nredmatrixeln_bj_btensorQ^(k)(2)n_dj_d\nendaligned\ntagV13126\nendequationSince the Clebsch–Gordan coefficients can be rewritten using 3j symbols and the 3j symbols vanish unless m_c + alpha - m_3 = m_d - alpha - m_b = 0, we havebeginequation\nalpha = m_a - m_c = m_d-m_b\nimplies\n-alpha + m_a + m_b = m_b + m_c\nendequationmultipole_expand_scalar_product\nmultipole_expand"
},

{
    "location": "multipole_expansions/#Spherical-tensors-1",
    "page": "Multipole expansions",
    "title": "Spherical tensors",
    "category": "section",
    "text": "For spherical tensors, we can insert the expression for the reduced matrix elements (V13.2.107) into (V13.1.26) above, to getbeginequation\nbeginaligned\nimplies\nmatrixelell_am_aell_bm_btensorC^(k)(1)cdottensorC^(k)(2)ell_cm_cell_dm_d\n=\nfrac1angrootell_aell_b\nsum_alpha(-)^-alpha\nC_ell_cm_ckalpha^ell_am_a\nC_ell_dm_dk-alpha^ell_bm_b\nredmatrixelell_atensorC^(k)(1)ell_c\nredmatrixelell_btensorC^(k)(2)ell_d\n=\nfracangrootell_cell_dangrootell_aell_b\nC_ell_c0k0^ell_a0\nC_ell_d0k0^ell_b0\nsum_alpha(-)^-alpha\nC_ell_cm_ckalpha^ell_am_a\nC_ell_dm_dk-alpha^ell_bm_b\n=beginalignedt\nsum_alpha\n(-)^-alpha + m_a + m_b\nangrootell_aell_bell_cell_d\ntimes\nbeginpmatrix\nell_ckell_a\nm_calpha-m_a\nendpmatrix\nbeginpmatrix\nell_dkell_b\nm_d-alpha-m_b\nendpmatrix\ntimes\nbeginpmatrix\nell_ckell_a000\nendpmatrix\nbeginpmatrix\nell_dkell_b000\nendpmatrix\nendaligned\nendaligned\nendequation DocTestSetup = nothing"
},

{
    "location": "energy_expressions/#",
    "page": "Energy expressions",
    "title": "Energy expressions",
    "category": "page",
    "text": ""
},

{
    "location": "energy_expressions/#Base.Matrix",
    "page": "Energy expressions",
    "title": "Base.Matrix",
    "category": "type",
    "text": "Matrix(op::QuantumOperator,\n       spin_cfgs::Vector{<:Configuration{<:SpinOrbital}}[, overlaps])\n\nGenerate the energy-expression associated with the quantum operator op, in the basis of the spin-configurations spin_cfgs, with an optional set of orbital overlaps, specifying any desired non-orthogonalities. The energy expression is generated in a basis-agnostic way by EnergyExpressions.jl and then multipole-expanded.\n\n\n\n\n\n"
},

{
    "location": "energy_expressions/#Energy-expressions-1",
    "page": "Energy expressions",
    "title": "Energy expressions",
    "category": "section",
    "text": "DocTestSetup = quote\n    using AngularMomentumAlgebra\n    using AtomicLevels\nendMatrix DocTestSetup = nothing"
},

]}
