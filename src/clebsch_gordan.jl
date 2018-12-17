struct ClebschGordan{A,B,C,j₁T,j₂T}
    j₁::A
    m₁::A
    j₂::B
    m₂::B
    J::C
    M::C
end

Base.zero(::Type{ClebschGordan{A,B,C,j₁T,j₂T}}) where {A,B,C,j₁T,j₂T} =
    ClebschGordan{A,B,C,j₁T,j₂T}(zero(A),zero(A),zero(B),zero(B),zero(C),zero(C))

ClebschGordan(j₁::A, m₁::A, j₂::B, m₂::B, J::C, M::C) where {A,B,C} =
    ClebschGordan{A,B,C,:J,:M}(j₁, m₁, j₂, m₂, J, M)

const ClebschGordanℓs{I<:Integer} = ClebschGordan{I,Rational{I},Rational{I},:ℓ,:s}
ClebschGordanℓs(j₁::I, m₁::I, j₂::R, m₂::R, J::R, M::R) where {I<:Integer,R<:Rational{I}} =
    ClebschGordanℓs{I}(j₁, m₁, j₂, m₂, J, M)

Base.convert(::Type{T}, cg::ClebschGordan) where {T<:Real} =
    clebschgordan(cg.j₁,cg.m₁,cg.j₂,cg.m₂,cg.J,cg.M)

Base.show(io::IO, cg::ClebschGordan) =
    write(io, "⟨$(rs(cg.j₁)),$(rs(cg.j₂));$(rs(cg.m₁)),$(rs(cg.m₂))|$(rs(cg.J)),$(rs(cg.M))⟩")

function Base.show(io::IO, cg::ClebschGordanℓs)
    spin = cg.m₂ > 0 ? "↑" : "↓"
    write(io, "⟨$(spectroscopic_label(cg.j₁));$(cg.m₁),$(spin)|$(rs(cg.J)),$(rs(cg.M))⟩")
end

export ClebschGordan, ClebschGordanℓs
