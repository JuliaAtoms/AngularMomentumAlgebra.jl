@testset "Angular momenta" begin
    function test_angular_momenta(𝐉, js, args...)
        # The matrix_element method for TotalAngularMomentum expects
        # (ℓ, s, j), whereas the other only expect j.
        makearg = isempty(args) ? j -> j : j -> (args..., j)
        for j ∈ js
            for mⱼ ∈ -j:j
                b = (makearg(j), mⱼ)
                for (±,∓) ∈ [(+,-),(-,+)]
                    a = (makearg(j), mⱼ)
                    # (V13.2.41,2)
                    @test matrix_element(a, TensorComponent(𝐉, 0), b) ≈ mⱼ
                    # (V13.2.42)
                    @test matrix_element(a, cartesian_tensor_component(𝐉, :z), b) ≈ mⱼ

                    mⱼ′ = mⱼ ± 1
                    abs(mⱼ′) > j && continue
                    a = (makearg(j), mⱼ′)

                    # (V13.2.41)
                    @test matrix_element(a, TensorComponent(𝐉, 0 ± 1), b) ≈ 0 ∓ √((j ± mⱼ + 1)*(j ∓ mⱼ)/2)
                    # (V13.2.42)
                    @test matrix_element(a, cartesian_tensor_component(𝐉, :x), b) ≈ 0.5*√((j ± mⱼ + 1)*(j ∓ mⱼ))
                    @test matrix_element(a, cartesian_tensor_component(𝐉, :y), b) ≈ 0 ∓ 0.5im*√((j ± mⱼ + 1)*(j ∓ mⱼ))
                end
            end
        end
    end
    @testset "Orbital angular momentum" begin
        𝐋 = OrbitalAngularMomentum()
        @test string(𝐋) == "𝐋̂⁽¹⁾"
        @test system(OrbitalAngularMomentum) == OrbitalAngularMomentumSubSystem()

        test_angular_momenta(𝐋, 0:20)
    end

    @testset "Spin angular momentum" begin
        𝐒 = SpinAngularMomentum()
        @test string(𝐒) == "𝐒̂⁽¹⁾"
        @test system(SpinAngularMomentum) == SpinSubSystem()

        test_angular_momenta(𝐒, half.(1:2:9))
    end

    @testset "Total angular momentum" begin
        𝐉 = TotalAngularMomentum()
        @test string(𝐉) == "𝐉̂⁽¹⁾"
        @test system(TotalAngularMomentum) == TotalAngularMomentumSubSystem()

        for ℓ = 0:0
            test_angular_momenta(𝐉, ℓ .+ (-half(1):half(1)), ℓ, half(1))
        end
    end
end
