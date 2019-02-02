@testset "One-body integrals" begin
    a = o"1s"
    b = o"2s"

    Iaa = OneBodyIntegral(a,a)
    Iab = OneBodyIntegral(a,b)

    @test !isdiagonal(Iab)
    @test isdiagonal(Iaa)

    @test diff(Iaa,b) isa OneBodyHamiltonian{Int}
    @test iszero(diff(Iaa,b))

    @test diff(Iaa, a) == OneBodyHamiltonian(conj(a))
    @test diff(Iaa, conj(a)) == OneBodyHamiltonian(a)

    @test iszero(diff(Iab, a))
    @test diff(Iab, b) == OneBodyHamiltonian(conj(a))
    @test iszero(diff(Iab, conj(b)))
    @test diff(Iab, conj(a)) == OneBodyHamiltonian(b)

    @test string(Iaa) == "I(1s)"
    @test string(Iab) == "I(1s, 2s)"
end
