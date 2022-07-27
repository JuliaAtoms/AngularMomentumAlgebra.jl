import EnergyExpressions: OrbitalMatrixElement, ContractedOperator
import AngularMomentumAlgebra: radial_integral
@testset "Coulomb interaction" begin
    g = CoulombInteraction()

    @testset "Multipoles" begin
        g3 = CoulombInteractionMultipole(3, g)
        @test string(g3) == "ĝ³"
        @test string(OrbitalMatrixElement([1,2],g3,[1,2])) == "F³(1,2)"
        @test string(OrbitalMatrixElement([1,2],g3,[2,1])) == "G³(1,2)"
        @test string(OrbitalMatrixElement([1,2],g3,[3,4])) == "R³(1,2;3,4)"
        @test string(ContractedOperator((1,),g3,(2,))) == "r⁻¹×Y³(1,2)"

        @test radial_integral([1,2], (3,g), [3,4]) == OrbitalMatrixElement([1,2],g3,[3,4])
    end

    @test system(CoulombTensor) == SpatialSubSystem()
end
