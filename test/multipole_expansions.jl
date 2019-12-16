using EnergyExpressions
import EnergyExpressions: OrbitalMatrixElement

@testset "Multipole expansions" begin
    @testset "Pure state" begin
        for (a,b,c,d) in [(SpinOrbital(o"1s", 0, half(1)),
                           SpinOrbital(o"1s", 0, -half(1)),
                           SpinOrbital(o"1s", 0, half(1)),
                           SpinOrbital(o"1s", 0, -half(1))),
                          (SpinOrbital(ro"1s", half(1)),
                           SpinOrbital(ro"1s", -half(1)),
                           SpinOrbital(ro"1s", half(1)),
                           SpinOrbital(ro"1s", -half(1))),
                          (SpinOrbital(o"1s", 0, half(1)),
                           SpinOrbital(o"2p", 1, half(1)),
                           SpinOrbital(o"1s", 0, half(1)),
                           SpinOrbital(o"2p", 1, half(1))),
                          (SpinOrbital(ro"1s", half(1)),
                           SpinOrbital(ro"2p", half(3)),
                           SpinOrbital(ro"1s", half(1)),
                           SpinOrbital(ro"2p", half(3)))]
            ome = OrbitalMatrixElement([a,b], CoulombInteraction(), [c,d])
            mome = multipole_expand(ome)
            @test length(mome.terms) == 1
            t = first(mome.terms)
            @test length(t.factors) == 1
            @test first(t.factors).o.k == 0
            @test t.coeff ≈ 1
        end
    end
end
