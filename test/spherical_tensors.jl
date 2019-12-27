@testset "Spherical tensors" begin
    @test system(TensorComponent(SphericalTensor(1), 1)) == OrbitalAngularMomentumSubSystem()
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
        for ℓ = 0:10
            for m = -ℓ:ℓ
                b = SpinOrbital(Orbital(:k, ℓ), m, half(1))
                for (±,∓) in [(+,-),(-,+)]
                    m′ = m ± 1
                    a = SpinOrbital(Orbital(:k, ℓ+1), m′, half(1))
                    @test dot(a, TensorComponent(𝐧, 0 ± 1), b) ≈ √((ℓ ± m + 1)*(ℓ ± m + 2)/
                                                                   (2*(2ℓ+1)*(2ℓ+3)))

                    m′ = m
                    a = SpinOrbital(Orbital(:k, ℓ+1), m′, half(1))
                    @test dot(a, TensorComponent(𝐧, 0), b) ≈ √((ℓ - m + 1)*(ℓ + m + 1)/
                                                               ((2ℓ+1)*(2ℓ+3)))

                    m′ = m ± 1
                    if abs(m′) ≤ ℓ-1
                        a = SpinOrbital(Orbital(:k, ℓ-1), m′, half(1))
                        @test dot(a, TensorComponent(𝐧, 0 ± 1), b) ≈ -√((ℓ ∓ m - 1)*(ℓ ∓ m)/
                                                                        (2*(2ℓ+1)*(2ℓ-1)))
                    end

                    m′ = m
                    if abs(m′) ≤ ℓ-1
                        a = SpinOrbital(Orbital(:k, ℓ-1), m′, half(1))
                        @test dot(a, TensorComponent(𝐧, 0), b) ≈ √((ℓ - m)*(ℓ + m)/
                                                                   ((2ℓ+1)*(2ℓ-1)))
                    end
                end
            end
        end
    end
end
