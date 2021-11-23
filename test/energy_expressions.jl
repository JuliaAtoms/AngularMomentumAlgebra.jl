@testset "Energy expressions" begin
    ome(a, op, b) = EnergyExpressions.OrbitalMatrixElement(a, op, b)
    ome(a::AbstractOrbital, op, b::AbstractOrbital) = ome([a], op, [b])
    c(a, op, b) = ContractedOperator(Tuple(e for e in a), op, Tuple(e for e in b))
    c(a::AbstractOrbital, op, b::AbstractOrbital) = c([a], op, [b])
    eq(orb, op, args...) = NBodyEquation(orb, op, args...)
    ceq(orb, op, args...) = eq(Conjugate(orb), op, args...)
    ov(a,b) = OrbitalOverlap(a,b)
    nbt(args...) = EnergyExpressions.NBodyTerm([args...], 1)

    @testset "Average-of-configuration" begin
        for gst = (c"1s²", rc"1s²")
            cfgs = [gst]

            os = sort(unique(reduce(vcat, orbitals.(cfgs))))
            o = os[1]
            v = csfs(cfgs)

            h = FieldFreeOneBodyHamiltonian()
            g = CoulombInteraction()
            H = h + g

            E = Matrix(H, cfgs)
            eqs = diff(E, conj.(os))

            @test length(eqs.equations) == 1
            eq1 = eqs.equations[1]
            @test eq1.orbital == o
            @test size(eq1.equation) == (1,1)
            eq1 = eq1.equation[1]

            @test eq1 == 2eq(o, h) + 2eq(o, c(o, CoulombInteractionMultipole(0, g), o))
        end
    end
end
