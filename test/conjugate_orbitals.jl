@testset "Conjugated orbitals" begin
    o = o"1s"
    co = conj(o)
    @test co == Conjugate(o)
    @test string(co) == "1s†"
end
