using AtomicLevels
using Symbolics

@testset "Lagrange multipliers" begin
    @test lagrange_multipliers([o"1s", o"2p"]) ==
        (LagrangeMultiplier(o"1s")*Braket(o"1s",o"1s") +
         LagrangeMultiplier(o"2p")*Braket(o"2p",o"2p"))
    # Sadly, order matters
    @test lagrange_multipliers([o"1s", o"2s", o"2p"]) ==
        (LagrangeMultiplier(o"1s")*Braket(o"1s",o"1s") +
         LagrangeMultiplier(o"2s")*Braket(o"2s",o"2s") +
         LagrangeMultiplier(o"2p")*Braket(o"2p",o"2p") +
         LagrangeMultiplier(o"1s",o"2s")*Braket(o"1s",o"2s") +
         LagrangeMultiplier(o"2s",o"1s")*Braket(o"2s",o"1s"))

    @test diff(lagrange_multipliers([o"1s", o"2p"]), o"1s", 1) ==
        2LagrangeMultiplier(o"1s")*Ket(o"1s")
    @test diff(lagrange_multipliers([o"1s", o"2p"]), o"2p", 1) ==
        2LagrangeMultiplier(o"2p")*Ket(o"2p")
    @test diff(lagrange_multipliers([o"1s", o"2p"]), o"2s", 1) == 0
end

@testset "Hartree–Fock equations" begin
    @testset "1s 2p 3P" begin
        eng = (DiagonalIntegral(o"1s") +
               DiagonalIntegral(o"2p") +
               DirectSlaterIntegral(0,o"1s", o"2p") -
               1//3*ExchangeSlaterIntegral(0,o"1s", o"2p") +
               lagrange_multipliers([o"1s", o"2p"]))
        @test diff(eng, o"1s", 1) ==
            (-Sym(:𝓛)*Ket(o"1s") +
             2/Sym(:r)*SlaterPotential(0,o"2p",o"2p")*Ket(o"1s") +
             -2//3*1/Sym(:r)*SlaterPotential(1,o"1s",o"2p")*Ket(o"2p") +
             2LagrangeMultiplier(o"1s")*Ket(o"1s"))
        
        @test diff(eng, o"2p", 1) ==
            (-Sym(:𝓛)*Ket(o"2p") +
             2/Sym(:r)*SlaterPotential(0,o"1s",o"1s")*Ket(o"2p") +
             -2//3*1/Sym(:r)*SlaterPotential(1,o"2p",o"1s")*Ket(o"1s") +
             2LagrangeMultiplier(o"2p")*Ket(o"2p"))
    end
end
