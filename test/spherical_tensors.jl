@testset "Spherical tensors" begin
    𝐧 = SphericalTensor(1)
    for ℓ = 0:10
        for m = -ℓ:ℓ
            b = SpinOrbital(Orbital(:m, ℓ), m, half(1))
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
