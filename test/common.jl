import AngularMomentumAlgebra: triangle_range, ∏, powneg1,
    jmⱼ, spin, @δ, couplings

@testset "Common routines" begin
    @testset "Triangle range" begin
        @test triangle_range(0,0) == [0]
        @test triangle_range(0,1) == [1]
        @test triangle_range(1,0) == [1]
        @test triangle_range(1,1) == [0,2]
        @test triangle_range(0,2) == [2]
        @test triangle_range(2,2) == [0,2,4]

        if VERSION ≥ v"1.7"
            # iseven(::Real) was introduced in 1.7
            @test triangle_range(half(1),half(1)) == [1]
            @test triangle_range(half(1),half(3)) == [2]
        end
    end

    @testset "Angular roots" begin
        @test ∏(0) == 1
        @test ∏(1) ≈ √3 atol=1e-14
        @test ∏(0,1) ≈ √3 atol=1e-14
        @test ∏(0,1,2) ≈ √3*√5 atol=1e-14
    end

    @test all(powneg1.(-10:2:10) .== 1)
    @test all(powneg1.(-9:2:10) .== -1)

    @testset "Orbital accessors" begin
        @test jmⱼ(so"1s₀β") == (0,0)
        @test jmⱼ(so"2p₋₁α") == (1,-1)
        @test jmⱼ(rso"1s(-1/2)") == (1/2,-1/2)
        @test jmⱼ(rso"2p(-3/2)") == (3/2,-3/2)

        @test spin(so"1s₀β") == -1/2
        @test spin(so"2p₋₁α") == 1/2
    end

    @testset "Kronecker δ" begin
        a = 1
        b = 2
        c = 3
        d = 4
        # We have to wrap the calls to @δ in function bodies, since they quick-return on inequality
        @test (() -> @δ((a,a)))() == 1
        @test (() -> @δ((a,b)))() == 0
        @test (() -> @δ((a,a), (b,b)))() == 1
        @test (() -> @δ((a,a), (b,c)))() == 0
        @test (() -> @δ((a,a), (b,c), (d,d)))() == 0
        @test (() -> @δ((a,a), (c,c), (d,d)))() == 1
    end

    @testset "Couplings" begin
        @test_throws DomainError couplings(1, 2, 1, 1)
        @test_throws DomainError couplings(1, 0, 1, 2)

        @test couplings(1, 1, 0, 0) == (1:1, 1)
        # min(|1-3|,1+1) = 2
        @test couplings(1, 1, 3, 1) == (2:4, 2)

        @test couplings(1, 1, 3, -1) == (2:4, 0)
        @test couplings(3, 0, 4, 0) == (2:2:6, 0)
        @test couplings(3, 0, 4, 1) == (1:7, 1)
    end
end
