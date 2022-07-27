@testset "Gradients" begin
    import AngularMomentumAlgebra: RadialGradientMatrixElement

    @test system(Gradient) == SpatialSubSystem()
    @test system(ReducedGradient) == SpatialSubSystem()

    @test string(RadialGradientMatrixElement(0)) == "∂ᵣ"
    @test string(RadialGradientMatrixElement(2)) == "(∂ᵣ + 2/r)"
    @test string(RadialGradientMatrixElement(-3)) == "(∂ᵣ - 3/r)"

    ∇ = Gradient()

    @test iszero(rme((1,0), ∇, (1,0)))
    @test rme((2,1), ∇, (3,0)) == 1RadialGradientMatrixElement(0)
    @test rme((1,0), ∇, (2,1)) == -1RadialGradientMatrixElement(2)
    # Risky? Until JuliaAtoms/EnergyExpressions.jl#18 is implemented,
    # we have no better option.
    @test rme((2,1), ∇, (3,2)) == -√2*RadialGradientMatrixElement(3)

    ∇̃ = ReducedGradient()

    @test rme((2,1), ∇̃, (3,0)) == 1RadialGradientMatrixElement(-1)
    @test rme((1,0), ∇̃, (2,1)) == -1RadialGradientMatrixElement(1)
    @test rme((2,1), ∇̃, (3,2)) == -√2*RadialGradientMatrixElement(2)
end
