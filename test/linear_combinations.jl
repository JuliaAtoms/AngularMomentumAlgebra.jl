import AngularMomentumAlgebra: @linearly_combinable

struct MySym
    s::Symbol
end
@linearly_combinable MySym

@testset "Linear combinations" begin
    x = MySym(:x)
    y = MySym(:y)

    w = x + y
    @test length(w) == 2
    @test eltype(w) == (MySym,Int)
    @test w.Ts == [x,y]
    @test w.coeffs == [1,1]

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
end
