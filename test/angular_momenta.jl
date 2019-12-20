@testset "Angular momenta" begin
    @testset "Total angular momentum" begin
        𝐉 = TotalAngularMomentum()
        for ℓ = 0:10
            for j = ℓ .+ (-half(1):half(1))
                for m = -j:j
                    b = SpinOrbital(RelativisticOrbital(:m, ℓ, j), m)
                    for (±,∓) in [(+,-),(-,+)]
                        a = SpinOrbital(RelativisticOrbital(:k, ℓ, j), m)
                        # (V13.2.41,2)
                        @test dot(a, TensorComponent(𝐉, 0), b) ≈ m
                        # (V13.2.42)
                        @test dot(a, cartesian_tensor_component(𝐉, :z), b) ≈ m

                        m′ = m ± 1
                        abs(m′) > j && continue
                        a = SpinOrbital(RelativisticOrbital(:k, ℓ, j), m′)

                        # (V13.2.41)
                        @test dot(a, TensorComponent(𝐉, 0 ± 1), b) ≈ 0 ∓ √((j ± m + 1)*(j ∓ m)/2)
                        # (V13.2.42)
                        @test dot(a, cartesian_tensor_component(𝐉, :x), b) ≈ 0.5*√((j ± m + 1)*(j ∓ m))
                        @test dot(a, cartesian_tensor_component(𝐉, :y), b) ≈ 0 ∓ 0.5im*√((j ± m + 1)*(j ∓ m))
                    end
                end
            end
        end
    end
end
