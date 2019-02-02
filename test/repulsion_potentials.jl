@testset "Repulsion potentials" begin
    a = o"1s"
    b = o"2s"
    o = o"3s"

    @test zero(DirectPotential) == DirectPotential(0,0,0)
    @test zero(ExchangePotential) == ExchangePotential(0,0,0)
    @test iszero(zero(DirectPotential))
    @test zero(DirectExchangePotentials) == DirectExchangePotentials(0,0,0)

    @test DirectPotential(a,a,a) == ExchangePotential(a,a,a)

    RP = RepulsionPotential{:general}(a,b,o)
    J = DirectPotential(a,b,o)
    K = ExchangePotential(a,b,o)
    JK = DirectExchangePotentials(a,b,o)

    @test !isdiagonal(JK)

    @test string(RP) == "[1s|2s]3s"
    @test string(J) == "Ĵ{1s;2s}3s"
    @test string(K) == "K̂{1s;2s}3s"
    @test string(JK) == "[1s||2s]3s"
end
