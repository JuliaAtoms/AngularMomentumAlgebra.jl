using EnergyExpressions
import EnergyExpressions: OrbitalMatrixElement

include("multipole_refs.jl")

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

    @testset "Direct interaction" begin
        orbitals = sos"k[s-d]"
        g = CoulombInteraction()
        for a in orbitals
            for b in orbitals
                ome = OrbitalMatrixElement([a,b], g, [a,b])
                mome = multipole_expand(ome)
                mref = ref_direct_multipoles(a, b)
                @test mome ≈ mref atol=1e-14
            end
        end
    end

    @testset "Exchange interaction" begin
        orbitals = sos"k[s-d]"
        g = CoulombInteraction()
        for a in orbitals
            for b in orbitals
                ome = OrbitalMatrixElement([a,b], g, [b,a])
                mome = multipole_expand(ome)
                if b.m[2] == a.m[2]
                    mref = ref_exchange_multipoles(a, b)
                    @test mome ≈ mref atol=1e-14
                else
                    @test iszero(mome)
                end
            end
        end
    end
end
