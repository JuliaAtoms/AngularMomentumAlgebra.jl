@testset "Tensor matrix elements" begin
    @testset "Wigner–Eckart" begin
        @test matrix_element((2, 1), TensorComponent(OrbitalAngularMomentum(), 1), (2, 0)) ≈ -√3
        @test matrix_element((0, 0), TensorComponent(SphericalTensor(1), 0), (1, 0)) ≈ 1/√3
        @test matrix_element(((1,half(1),half(1)), -half(1)),
                             TensorComponent(TotalAngularMomentum(), -1),
                             ((1,half(1),half(1)), half(1))) ≈ 1/√2
    end

    @testset "Wigner–Eckart for orbitals" begin
        @test dot(so"3d₁α", TensorComponent(OrbitalAngularMomentum(), 1), so"3d₀α") ≈ -√3
        @test dot(so"2s₀α", TensorComponent(SphericalTensor(1), 0), so"2p₀α") ≈ 1/√3
        # The spherical tensor only couples the angular dimension, not r
        @test iszero(dot(so"1s₀α", TensorComponent(SphericalTensor(1), 0), so"2p₀α"))
        @test dot(rso"2p-(-1/2)", TensorComponent(TotalAngularMomentum(), -1), rso"2p-(1/2)") ≈ 1/√2
        @test dot(rso"2p(-1/2)", TensorComponent(TotalAngularMomentum(), -1), rso"2p(1/2)") ≈ √2

        # Orthogonality in spin
        @test iszero(dot(so"3d₁α", TensorComponent(OrbitalAngularMomentum(), 1), so"3d₀β"))
        @test iszero(dot(so"2s₀α", TensorComponent(SphericalTensor(1), 0), so"2p₀β"))
    end

    𝐋 = OrbitalAngularMomentum()
    𝐋² = 𝐋⋅𝐋
    𝐋₀ = TensorComponent(𝐋, 0)

    𝐒 = SpinAngularMomentum()
    𝐒² = 𝐒⋅𝐒
    𝐒₀ = TensorComponent(𝐒, 0)

    𝐉 = TotalAngularMomentum()
    𝐉² = 𝐉⋅𝐉
    𝐉₀ = TensorComponent(𝐉, 0)
    𝐉₁ = TensorComponent(𝐉, 1)

    𝐂⁰ = SphericalTensor(0)

    @testset "Scalar product tensor" begin
        @test matrix_element((half(1), half(1)),
                             𝐒², (half(1), half(1))) ≈ 1/2*(1/2 + 1)
        @test matrix_element(((1, half(1), half(3)), half(3)),
                             𝐉², ((1, half(1), half(3)), half(3))) ≈ 3/2*(3/2+1)

        @test matrix_element2((0, 0), (0, 0), 𝐂⁰⋅𝐂⁰, (0,0), (0, 0)) ≈ 1
        @test matrix_element2((1, 1), (half(1), half(1)),
                              𝐋⋅𝐒, (1,1), (half(1), half(1))) ≈ 1/2
        @test matrix_element2((1, half(1), half(3), half(3)),
                              𝐋⋅𝐒, (1, half(1), half(3), half(3))) ≈ 1/2
        @test iszero(matrix_element2((1, half(1), half(3), half(3)),
                                     𝐋⋅𝐒, (1, half(1), half(5), half(3))))
        @test iszero(matrix_element2((1, half(1), half(3), half(3)),
                                     𝐋⋅𝐒, (1, half(1), half(3), half(1))))

        @test dot(so"2p₁α", 𝐋², so"2p₁α") ≈ 1*(1+1)
        @test dot(so"2p₁α", 𝐒², so"2p₁α") ≈ 1/2*(1/2+1)
        @test iszero(dot(so"2p₁α", 𝐋², so"2p₁β"))
        @test iszero(dot(so"2p₁α", 𝐒², so"2p₀α"))
    end

    @testset "Tensor acts on entire system" begin
        @test matrix_element((1, 1), 𝐋₀, (1, 1)) ≈ 1
        @test matrix_element((half(1),half(1)), 𝐒₀, (half(1),half(1))) ≈ 1/2

        @testset "Uncoupled basis functions" begin
            @test matrix_element((1,1), (half(1),half(1)),
                                 𝐉₀, (1,1), (half(1), half(1))) ≈ 3/2

            @test matrix_element((1,-1), (half(1),half(1)),
                                 𝐉₀, (1,-1), (half(1), half(1))) ≈ -1/2

            @test matrix_element((1,1), (half(1),half(1)),
                                 𝐉₁, (1,0), (half(1), half(1))) ≈ -1

            @test matrix_element((1,1), (half(1),half(1)),
                                 𝐉², (1,1), (half(1), half(1))) ≈ 3/2*(3/2+1)

            @test dot(so"2p₁α", 𝐉₀, so"2p₁α") ≈ 3/2
            @test dot(so"2p₋₁α", 𝐉₀, so"2p₋₁α") ≈ -1/2
            @test dot(so"2p₁α", 𝐉₁, so"2p₀α") ≈ -1
            @test dot(so"2p₁α", 𝐉², so"2p₁α") ≈ 3/2*(3/2+1)
        end

        @testset "Coupled basis functions" begin
            @test dot(rso"1s(1/2)", 𝐉², rso"1s(1/2)") ≈ 1/2*(1/2+1)
            @test dot(rso"2p-(1/2)", 𝐉², rso"2p-(1/2)") ≈ 1/2*(1/2+1)
            @test dot(rso"2p-(-1/2)", 𝐉², rso"2p-(-1/2)") ≈ 1/2*(1/2+1)
            @test dot(rso"2p(1/2)", 𝐉², rso"2p(1/2)") ≈ 3/2*(3/2+1)
            @test dot(rso"2p(-1/2)", 𝐉², rso"2p(-1/2)") ≈ 3/2*(3/2+1)
            @test dot(rso"2p(3/2)", 𝐉², rso"2p(3/2)") ≈ 3/2*(3/2+1)
            @test dot(rso"2p(-3/2)", 𝐉², rso"2p(-3/2)") ≈ 3/2*(3/2+1)
        end
    end

    @testset "Tensor acts on subsystem" begin
        @testset "Coupled basis functions" begin
            @test matrix_element((1, half(1), half(3), half(3)),
                                 𝐋₀, (1, half(1), half(3), half(3))) ≈ 1
            @test matrix_element((half(1), 1, half(3), half(3)),
                                 𝐒₀, (half(1), 1, half(3), half(3))) ≈ 1/2

            @test iszero(dot(rso"1s(1/2)", 𝐋₀, rso"1s(1/2)"))
            # See table on page 35 of Lindgren & Morrison (1986)
            @test dot(rso"2p-(1/2)", 𝐋₀, rso"2p-(1/2)") ≈ 2/3 # (1/3*0 + 2/3*1)
            @test dot(rso"2p-(-1/2)", 𝐋₀, rso"2p-(-1/2)") ≈ -2/3 # (1/3*0 + 2/3*(-1))
            @test dot(rso"2p(1/2)", 𝐋₀, rso"2p(1/2)") ≈ 1/3 # (2/3*0 + 1/3*1)
            @test dot(rso"2p(-1/2)", 𝐋₀, rso"2p(-1/2)") ≈ -1/3 # (2/3*0 + 1/3*(-1))
            @test dot(rso"2p(3/2)", 𝐋₀, rso"2p(3/2)") ≈ 1
            @test dot(rso"2p(-3/2)", 𝐋₀, rso"2p(-3/2)") ≈ -1

            @test iszero(dot(rso"1s(1/2)", 𝐋², rso"1s(1/2)"))
            @test dot(rso"2p-(1/2)", 𝐋², rso"2p-(1/2)") ≈ 1*(1+1)
            @test dot(rso"2p-(-1/2)", 𝐋², rso"2p-(-1/2)") ≈ 1*(1+1)
            @test dot(rso"2p(1/2)", 𝐋², rso"2p(1/2)") ≈ 1*(1+1)
            @test dot(rso"2p(-1/2)", 𝐋², rso"2p(-1/2)") ≈ 1*(1+1)
            @test dot(rso"2p(3/2)", 𝐋², rso"2p(3/2)") ≈ 1*(1+1)
            @test dot(rso"2p(-3/2)", 𝐋², rso"2p(-3/2)") ≈ 1*(1+1)

            @test iszero(dot(rso"1s(1/2)", 𝐋⋅𝐒, rso"1s(1/2)"))
            @test dot(rso"2p-(1/2)", 𝐋⋅𝐒, rso"2p-(1/2)") ≈ -1
            @test dot(rso"2p-(-1/2)", 𝐋⋅𝐒, rso"2p-(-1/2)") ≈ -1
            @test dot(rso"2p(1/2)", 𝐋⋅𝐒, rso"2p(1/2)") ≈ 1/2
            @test dot(rso"2p(-1/2)", 𝐋⋅𝐒, rso"2p(-1/2)") ≈ 1/2
            @test dot(rso"2p(3/2)", 𝐋⋅𝐒, rso"2p(3/2)") ≈ 1/2
            @test dot(rso"2p(-3/2)", 𝐋⋅𝐒, rso"2p(-3/2)") ≈ 1/2
        end

        @testset "Immaterial other subsystem" begin
            a = ((0, half(1), half(1)), (), half(1), half(1))
            @test matrix_element(a, 𝐉₀, a) ≈ 1/2
            b = ((0, half(1), half(3)), (), half(1), half(1))
            @test_throws ArgumentError matrix_element(a, 𝐉₀, b)
        end
    end
end
