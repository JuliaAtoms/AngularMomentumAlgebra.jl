@testset "Energy expressions" begin
    import AngularMomentumAlgebra: integrate_spinor
    import EnergyExpressions: OrbitalMatrixElement, NBodyTerm, NBodyMatrixElement

    ome(a, op, b) = OrbitalMatrixElement(a, op, b)
    ome(a::AbstractOrbital, op, b::AbstractOrbital) = ome([a], op, [b])
    c(a, op, b) = ContractedOperator(Tuple(e for e in a), op, Tuple(e for e in b))
    c(a::AbstractOrbital, op, b::AbstractOrbital) = c([a], op, [b])
    eq(orb, op, args...) = NBodyEquation(orb, op, args...)
    ceq(orb, op, args...) = eq(Conjugate(orb), op, args...)
    ov(a,b) = OrbitalOverlap(a,b)
    nbt(args...) = NBodyTerm([args...], 1)

    @testset "Integrate spinor" begin
        @test integrate_spinor(1) == 1
        @test integrate_spinor(3) == 3

        @testset "Direct interaction" begin
            orbitals = sos"k[s-d]"
            g = CoulombInteraction()
            for a in orbitals
                for b in orbitals
                    mome = integrate_spinor(ome([a,b], g, [a,b]))
                    mref = ref_direct_multipoles(a, b)
                    @test mome ≈ mref atol=1e-14
                end
            end
        end

        @testset "Zero-body identity operator" begin
            I₀ = IdentityOperator{0}()

            me = ome([], I₀, [])
            @test isone(integrate_spinor(me))

            cfgs = [sc"1s₀α kp₀α"]
            oo = ov(so"kp₀α", so"kp₀α")
            E = Matrix(I₀, cfgs, [oo])
            @test E == NBodyMatrixElement[nbt(oo)][:,:]
        end
    end

    @testset "Multi-configurational" begin
        cfgs = spin_configurations([c"1s2", c"1s 2s"])
        g = CoulombInteraction()
        g⁰ = CoulombInteractionMultipole(0, g)
        E = Matrix(g, cfgs)
        @test size(E) == (5,5)
        x,y,z,w = sos"1[s] 2[s]"
        R = (a,b,c,d) -> NBodyMatrixElement([NBodyTerm([ome([a,b],g⁰,[c,d])],1)])
        @test all(E .== NBodyMatrixElement[R(x,y,x,y) 0 R(x,y,x,w) -R(x,y,z,y) 0
                                           0 -R(x,z,z,x)+R(x,z,x,z) 0 0 0
                                           R(x,w,x,y) 0 R(x,w,x,w) -R(x,w,z,y) 0
                                           -R(y,z,y,x) 0 -R(y,z,w,x) R(y,z,y,z) 0
                                           0 0 0 0 -R(y,w,w,y)+R(y,w,y,w)])
    end

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
