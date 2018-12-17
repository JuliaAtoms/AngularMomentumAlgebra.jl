const spectroscopic = "spdfghiklmnoqrtuv"
spectroscopic_label(ℓ) =
    ℓ + 1 ≤ length(spectroscopic) ? spectroscopic[ℓ+1] : "[$(ℓ)]"

# Nicer string representation for rationals
rs(r::Number) = "$(r)"
rs(r::Rational) = denominator(r) == 1 ? "$(numerator(r))" : "$(numerator(r))/$(denominator(r))"
