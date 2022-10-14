@testset "Gradients" begin
    import AngularMomentumAlgebra: RadialGradientOperator

    @test system(Gradient) == SpatialSubSystem()
    @test system(ReducedGradient) == SpatialSubSystem()

    @test string(RadialGradientOperator(0)) == "∂ᵣ"
    @test string(RadialGradientOperator(2)) == "(∂ᵣ + 2/r)"
    @test string(RadialGradientOperator(-3)) == "(∂ᵣ - 3/r)"

    ∇ = Gradient()

    @test iszero(rme((1,0), ∇, (1,0)))
    @test rme((2,1), ∇, (3,0)) == 1RadialGradientOperator(0)
    @test rme((1,0), ∇, (2,1)) == -1RadialGradientOperator(2)
    # Risky? Until JuliaAtoms/EnergyExpressions.jl#18 is implemented,
    # we have no better option.
    @test rme((2,1), ∇, (3,2)) == -√2*RadialGradientOperator(3)

    ∇̃ = ReducedGradient()

    @test rme((2,1), ∇̃, (3,0)) == 1RadialGradientOperator(-1)
    @test rme((1,0), ∇̃, (2,1)) == -1RadialGradientOperator(1)
    @test rme((2,1), ∇̃, (3,2)) == -√2*RadialGradientOperator(2)
end
