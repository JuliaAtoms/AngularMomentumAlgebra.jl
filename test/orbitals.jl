@testset "Orbitals" begin
    function compare(a, b::Missing)
        @warn "Comparing with Missing" a
        false
    end
    function compare(a::Missing, b)
        @warn "Comparing with Missing" b
        false
    end
    function compare(a, b)
        @test a == b
        a == b
    end
    function compare(a::Union{Tuple,AbstractVector}, b::Union{Tuple,AbstractVector})
        @test all(compare(i,j) for (i,j) in zip(a,b))
        all(compare(i,j) for (i,j) in zip(a,b))
    end
    compare(::Missing, ::Missing) = true

    @testset "SpatialSubSystem" begin
        S = SpatialSubSystem()
        @test quantum_numbers(S, so"1s₀α") == ((1,0),0)
        @test quantum_numbers(S, so"1s₀β") == ((1,0),0)
        @test quantum_numbers(S, so"2p₋₁α") == ((2,1),-1)
        @test quantum_numbers(S, so"3d₂β") == ((3,2),2)
        @test compare(quantum_numbers(S, rso"1s(1/2)"), ((1,0), missing))

        @test other_quantum_numbers(S, so"1s₀α") == (1/2,1/2)
        @test other_quantum_numbers(S, so"1s₀β") == (1/2,-1/2)
        @test other_quantum_numbers(S, so"2p₋₁α") == (1/2,1/2)
        @test other_quantum_numbers(S, so"3d₂β") == (1/2,-1/2)
        @test compare(other_quantum_numbers(S, rso"1s(1/2)"), (1/2,missing))
    end

    @testset "OrbitalAngularMomentumSubSystem" begin
        S = OrbitalAngularMomentumSubSystem()
        @test quantum_numbers(S, so"1s₀α") == (0,0)
        @test quantum_numbers(S, so"1s₀β") == (0,0)
        @test quantum_numbers(S, so"2p₋₁α") == (1,-1)
        @test quantum_numbers(S, so"3d₂β") == (2,2)
        @test compare(quantum_numbers(S, rso"1s(1/2)"), (0, missing))

        @test other_quantum_numbers(S, so"1s₀α") == ((1,1/2),1/2)
        @test other_quantum_numbers(S, so"1s₀β") == ((1,1/2),-1/2)
        @test other_quantum_numbers(S, so"2p₋₁α") == ((2,1/2),1/2)
        @test other_quantum_numbers(S, so"3d₂β") == ((3,1/2),-1/2)
        @test compare(other_quantum_numbers(S, rso"1s(1/2)"), ((1,1/2),missing))
    end

    @testset "SpinSubSystem" begin
        S = SpinSubSystem()
        @test quantum_numbers(S, so"1s₀α") == (1/2, 1/2)
        @test quantum_numbers(S, so"1s₀β") == (1/2, -1/2)
        @test compare(quantum_numbers(S, rso"1s(1/2)"), (1/2,missing))

        @test other_quantum_numbers(S, so"1s₀α") == ((1,0),0)
        @test other_quantum_numbers(S, so"1s₀β") == ((1,0),0)
        @test other_quantum_numbers(S, so"2p₋₁α") == ((2,1),-1)
        @test other_quantum_numbers(S, so"3d₂β") == ((3,2),2)
        @test compare(other_quantum_numbers(S, rso"1s(1/2)"), ((1,0),missing))
    end

    @testset "TotalAngularMomentumSubSystem" begin
        S = TotalAngularMomentumSubSystem()
        @test quantum_numbers(S, so"1s₀α") == ((0,0), (1/2,1/2))
        @test quantum_numbers(S, rso"1s(1/2)") == ((0,1/2,1/2), 1/2)
        @test quantum_numbers(S, rso"2p(-3/2)") == ((1,1/2,3/2), -3/2)

        @test compare(other_quantum_numbers(S, so"1s₀α"), ((), missing))
        @test compare(other_quantum_numbers(S, rso"1s(1/2)"), ((), missing))
        @test compare(other_quantum_numbers(S, rso"2p(-3/2)"), ((), missing))
    end

    @testset "FullSystem" begin
        S = FullSystem()
        @test quantum_numbers(S, so"1s₀α") == (((1,0),0), (1/2,1/2))
        @test quantum_numbers(S, rso"1s(1/2)") == ((1,0,1/2,1/2), 1/2)
        @test quantum_numbers(S, rso"2p(-3/2)") == ((2,1,1/2,3/2), -3/2)

        @test compare(other_quantum_numbers(S, so"1s₀α"), ((), missing))
        @test compare(other_quantum_numbers(S, rso"1s(1/2)"), ((), missing))
        @test compare(other_quantum_numbers(S, rso"2p(-3/2)"), ((), missing))
    end

    @testset "Two orbitals" begin
        S = SpatialSubSystem()
        @test quantum_numbers(S, so"1s₀α", so"1s₀β") == (((1,0),0), ((1,0),0))
        @test other_quantum_numbers(S, so"1s₀α", so"1s₀β") == ((1/2,1/2), (1/2,-1/2))
    end
end
