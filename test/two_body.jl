@testset "Two-body integrals" begin
    a = o"1s"
    b = o"2s"
    c = o"2p"
    d = o"3s"

    @testset "General Slater integrals" begin
        R = GeneralRepulsionIntegral(a,b,c,d)

        @test diff(R,a) isa RepulsionPotential{:general}
        @test iszero(diff(R,a))
        @test diff(R,b) isa RepulsionPotential{:general}
        @test iszero(diff(R,b))
        @test diff(R,conj(c)) isa RepulsionPotential{:general}
        @test iszero(diff(R,conj(c)))
        @test diff(R,conj(d)) isa RepulsionPotential{:general}
        @test iszero(diff(R,conj(d)))

        @test diff(R,conj(a)) == RepulsionPotential{:general}(b,d,c)
        @test diff(R,conj(b)) == RepulsionPotential{:general}(a,c,d)

        @test diff(R,c) == RepulsionPotential{:general}(b,d,conj(a))
        @test diff(R,d) == RepulsionPotential{:general}(a,c,conj(b))

        @test string(R) == "[1s 2s|2p 3s]"
    end

    @testset "Direct Slater integrals" begin
        F = DirectIntegral(a,b)

        @test !isdiagonal(F)
        @test isdiagonal(DirectIntegral(a,a))

        @test diff(F,c) isa DirectPotential
        @test iszero(diff(F,c))
        @test diff(F,conj(c)) isa DirectPotential
        @test iszero(diff(F,conj(c)))

        @test diff(F,a) == DirectPotential(b,b,conj(a))
        @test diff(F,b) == DirectPotential(a,a,conj(b))
        @test diff(F,conj(a)) == DirectPotential(b,b,a)
        @test diff(F,conj(b)) == DirectPotential(a,a,b)

        @test string(F) == "F(1s, 2s)"
    end

    @testset "Exchange Slater integrals" begin
        G = ExchangeIntegral(a,b)

        @test diff(G,c) isa ExchangePotential
        @test iszero(diff(G,c))
        @test diff(G,conj(c)) isa ExchangePotential
        @test iszero(diff(G,conj(c)))

        @test diff(G,a) == ExchangePotential(b,b,conj(a))
        @test diff(G,b) == ExchangePotential(a,a,conj(b))
        @test diff(G,conj(a)) == ExchangePotential(b,b,a)
        @test diff(G,conj(b)) == ExchangePotential(a,a,b)

        @test string(G) == "G(1s, 2s)"
    end

    @testset "Composite Direct–Exchange integrals" begin
        FG = TwoBodyIntegral(a,b,c,d)

        @test !isdiagonal(FG)
        @test isdiagonal(TwoBodyIntegral(a,a,a,a))

        @test diff(FG,a) isa DirectExchangePotentials
        @test iszero(diff(FG,a))
        @test diff(FG,b) isa DirectExchangePotentials
        @test iszero(diff(FG,b))
        @test diff(FG,conj(c)) isa DirectExchangePotentials
        @test iszero(diff(FG,conj(c)))
        @test diff(FG,conj(d)) isa DirectExchangePotentials
        @test iszero(diff(FG,conj(d)))

        @test diff(FG,conj(a)) == DirectExchangePotentials(b,d,c)
        @test diff(FG,conj(b)) == DirectExchangePotentials(a,c,d)
        @test diff(FG,c) == DirectExchangePotentials(b,d,conj(a))
        @test diff(FG,d) == DirectExchangePotentials(a,c,conj(b))

        @test string(FG) == "[1s 2s||2p 3s]"
    end
end
