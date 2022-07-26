import AngularMomentumAlgebra: @linearly_combinable, LinearCombination

struct MySym
    s::Symbol
end
@linearly_combinable MySym

@testset "Linear combinations" begin
    x = MySym(:x)
    y = MySym(:y)
    g = MySym(:g)

    w = x + y
    @test w isa LinearCombination
    @test length(w) == 2
    @test eltype(w) == (MySym,Int)
    @test w.Ts == [x,y]
    @test w.coeffs == [1,1]

    @test w == x + y
    @test hash(w) == hash(x + y)

    @test w + g == LinearCombination([x,y,g], [1,1,1])
    @test g + w == LinearCombination([g,x,y], [1,1,1])
    @test w - g == LinearCombination([x,y,g], [1,1,-1])
    @test g - w == LinearCombination([g,x,y], [1,-1,-1])
    @test w + g ≈ (1+1e-9)*LinearCombination([x,y,g], [1,1,1])

    @test iszero(LinearCombination(MySym[], Int[]))
    @test iszero(zero(w))

    m = 2x
    @test m isa LinearCombination
    @test m.Ts == [x]
    @test m.coeffs == [2]
    n = x*2
    @test n isa LinearCombination
    @test n.Ts == [x]
    @test n.coeffs == [2]

    @test -x isa LinearCombination
    @test (-x).coeffs == [-1]
    @test length((-x)) == 1

    z = 4x - 5y
    z′ = -z
    z′′ = 2z

    @test z.Ts == [x,y]
    @test z.coeffs == [4, -5]

    @test z′.Ts == [x,y]
    @test z′.coeffs == [-4, 5]

    @test z′′.Ts == [x,y]
    @test z′′.coeffs == [8, -10]
    
    @test all((z′′/2).coeffs .== z.coeffs)

    @test [(t,c) for (t,c) in z] == [(x,4), (y,-5)]

    @test string(w) == "MySym(:x) + MySym(:y)"
    @test string(-w) == "- MySym(:x) - MySym(:y)"
    @test string(-x) == "- MySym(:x)"
end
