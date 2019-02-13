"""
    clebsch_gordan_condon_shortley(j₁, m₁, j₂, m₂, j₃, m₃=m₁+m₂)

Calculate the vector coupling coefficient `⟨j₁ m₁, j₂ m₂|j₃ m₃⟩`
according to the Condon–Shortley phase convention and using the
definition of Eq. (8.1.12) in Varshalovich (1988).
"""
clebsch_gordan_condon_shortley(j₁, m₁, j₂, m₂, j₃, m₃=m₁+m₂) =
    powneg1(j₁ - j₂ + m₃)*∏(j₃)*wigner3j(j₁, j₂, j₃,
                                         m₁, m₂, -m₃)

export clebsch_gordan_condon_shortley
