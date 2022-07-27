@testset "Spherical tensors" begin
    import AngularMomentumAlgebra: ranks, RadialMatrixElement

    approx_comp(a, b) = a ≈ b
    function approx_comp(a::EnergyExpressions.LinearCombinationOperator,
                         b::EnergyExpressions.LinearCombinationOperator)
        length(a.operators) == length(b.operators) || return false
        for ((ao,ac),(bo,bc)) in zip(a.operators,b.operators)
            ao == bo && ac ≈ bc || return fals
        end
        true
    end

    function test_components(𝐓¹, ℓs, μ=1)
        for ℓ = ℓs
            for m = -ℓ:ℓ
                b = SpinOrbital(Orbital(:k, ℓ), m, half(1))
                for (±,∓) in [(+,-),(-,+)]
                    m′ = m ± 1
                    a = SpinOrbital(Orbital(:k, ℓ+1), m′, half(1))
                    @test approx_comp(dot(a, TensorComponent(𝐓¹, 0 ± 1), b),
                                      √((ℓ ± m + 1)*(ℓ ± m + 2)/
                                          (2*(2ℓ+1)*(2ℓ+3)))*μ)

                    m′ = m
                    a = SpinOrbital(Orbital(:k, ℓ+1), m′, half(1))
                                      @test approx_comp(dot(a, TensorComponent(𝐓¹, 0), b),
                                                        √((ℓ - m + 1)*(ℓ + m + 1)/
                                                            ((2ℓ+1)*(2ℓ+3)))*μ)

                    m′ = m ± 1
                    if abs(m′) ≤ ℓ-1
                        a = SpinOrbital(Orbital(:k, ℓ-1), m′, half(1))
                                      @test approx_comp(dot(a, TensorComponent(𝐓¹, 0 ± 1), b),
                                                        -√((ℓ ∓ m - 1)*(ℓ ∓ m)/
                                                            (2*(2ℓ+1)*(2ℓ-1)))*μ)
                    end

                    m′ = m
                    if abs(m′) ≤ ℓ-1
                        a = SpinOrbital(Orbital(:k, ℓ-1), m′, half(1))
                                          @test approx_comp(dot(a, TensorComponent(𝐓¹, 0), b),
                                                            √((ℓ - m)*(ℓ + m)/
                                                                ((2ℓ+1)*(2ℓ-1)))*μ)
                    end
                end
            end
        end
    end

    @test system(TensorComponent(SphericalTensor(1), 1)) == OrbitalAngularMomentumSubSystem()
    @test ranks(so"1s₀α", SphericalTensor, so"3d₂α") == 2:2:2
    @test ranks(so"3d₀α", SphericalTensor, so"4f₂α") == 1:2:5

    @testset "Reduced matrix elements" begin
        for k = 0:6
            𝐂ᵏ = SphericalTensor(k)
            for (a,b) in [(0,k), (k,0), (1,k+1), (k+1,1)]
                @test !iszero(a, 𝐂ᵏ, b)
                # Since this is verbatim the definition of the reduced
                # matrix element, this is more a test of the tensor
                # DSL producing the correct code.
                @test rme(a, 𝐂ᵏ, b) ≈ AngularMomentumAlgebra.∏(b)*clebschgordan(b, 0, k, 0, a)
            end
        end
    end
    @testset "Components" begin
        𝐧 = SphericalTensor(1)
        test_components(𝐧, 0:10)

        @testset "Cartesian components" begin
            𝐂¹ = SphericalTensor(1)
            C₁ = TensorComponent(𝐂¹, 1)
            C₀ = TensorComponent(𝐂¹, 0)
            C₋₁ = TensorComponent(𝐂¹, -1)

            x̂,ŷ,ẑ = AngularMomentumAlgebra.Dipoles.𝐫̂
            @test x̂ ≈ (-C₁ + C₋₁)/√2
            @test ŷ ≈ im*(C₁ + C₋₁)/√2
            @test ẑ == C₀
        end
    end

    @testset "Dipoles" begin
        𝐃 = Dipole()
        D₁ = TensorComponent(𝐃, 1)
        D₀ = TensorComponent(𝐃, 0)
        D₋₁ = TensorComponent(𝐃, -1)
        @test system(𝐃) == SpatialSubSystem()

        r = RadialMatrixElement()
        @test string(r) == "r"

        @test iszero(dot(so"1s₀α", D₀, so"3d₀α"))

        test_components(𝐃, 0:10, r)

        @testset "Cartesian components" begin
            x,y,z = AngularMomentumAlgebra.Dipoles.𝐫
            @test x ≈ (-D₁ + D₋₁)/√2
            @test y ≈ im*(D₁ + D₋₁)/√2
            @test z == D₀
        end
    end
end
