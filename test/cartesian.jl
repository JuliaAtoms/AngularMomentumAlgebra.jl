@testset "Cartesian tensor components" begin
    𝐉 = TotalAngularMomentum()
    J₁ = TensorComponent(𝐉, 1)
    J₀ = TensorComponent(𝐉, 0)
    J₋₁ = TensorComponent(𝐉, -1)

    Jx = cartesian_tensor_component(𝐉, :x)
    @test Jx ≈ (-J₁ + J₋₁)/√2
    Jy = cartesian_tensor_component(𝐉, :y)
    @test Jy ≈ im*(J₁ + J₋₁)/√2
    Jz = cartesian_tensor_component(𝐉, :z)
    @test Jz == J₀
    @test_throws ArgumentError cartesian_tensor_component(𝐉, :w)
end
