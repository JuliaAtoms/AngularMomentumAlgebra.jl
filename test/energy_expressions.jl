@testset "Energy expressions" begin
    @testset "Configurational averages" begin
        @testset "Expansion coefficients" begin
            # Reference data taken from Table 2-1 of: Froese Fischer,
            # Charlotte (1977). The Hartree-Fock Method for Atoms: A
            # Numerical Approach. New York: Wiley.
            @testset "Equivalent electrons" begin
                ref_data = [
                    0 => [],
                    1 => [(2,-2/25)],
                    2 => [(2,-2/63), (4,-2/63)],
                    3 => [(2,-4/195), (4,-2/143), (6,-100/5577)],
                    4 => [(2,-20/1309), (4,-162/17017), (6,-20/2431), (8,-490/41327)]
                ]
                for (ℓ, d) in ref_data
                    data = AngularMomentumAlgebra.f_av(ℓ)
                    @test length(data) == length(d)
                    ℓ == 0 && continue
                    @test first.(data) == first.(d)
                    @test last.(data) ≈ last.(d)
                end
            end
            @testset "Non-equivalent electrons" begin
                ref_data = [
                    (0,0) => [(0, -1/2)],
                    (0,1) => [(1, -1/6)],
                    (0,2) => [(2, -1/10)],
                    (0,3) => [(3, -1/14)],
                    (0,4) => [(4, -1/18)],
                    (1,1) => [(0, -1/6), (2, -1/15)],
                    (1,2) => [(1, -1/15), (3, -3/70)],
                    (1,3) => [(2, -3/70), (4, -2/63)],
                    (1,4) => [(3, -2/63), (5, -5/198)],
                    (2,2) => [(0, -1/10), (2, -1/35), (4, -1/35)],
                    (2,3) => [(1, -3/70), (3, -2/105), (5, -5/231)],
                    (2,4) => [(2, -1/35), (4, -10/693), (6, -5/286)],
                    (3,3) => [(0, -1/14), (2, -2/105), (4, -1/77), (6, -50/3003)],
                    (3,4) => [(1, -2/63), (3, -1/77), (5, -10/1001), (7, -35/2574)],
                    (4,4) => [(0, -1/18), (2, -10/693), (4, -9/1001), (6, -10/1287), (8, -245/21879)]
                ]
                for ((ℓ,ℓ′), d) in ref_data
                    data = AngularMomentumAlgebra.g_av(ℓ,ℓ′)
                    @test length(data) == length(d)
                    @test first.(data) == first.(d)
                    @test last.(data) ≈ last.(d)
                end
            end
        end
    end
end
