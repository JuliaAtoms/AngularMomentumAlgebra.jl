@testset "Tensors" begin
    import AngularMomentumAlgebra: components, component,
    RadialOperator, RadialGradientOperator,
    OrbitalRadialOverlap, OrbitalRadialMatrixElement, radial_integral, integrate_spinor
    import EnergyExpressions: NBodyTerm, NBodyMatrixElement, OrbitalMatrixElement
    @test system(Tensor) == FullSystem()

    𝐂⁵ = SphericalTensor(5)
    𝐋 = OrbitalAngularMomentum()
    𝐒 = SpinAngularMomentum()
    𝐉 = TotalAngularMomentum()
    𝐉₀ = TensorComponent(𝐉, 0)

    @test 𝐂⁵ == SphericalTensor{5}()
    @test rank(𝐂⁵) == 5
    @test components(𝐂⁵) == -5:5
    @test string(𝐂⁵) == "𝐂̂⁽⁵⁾"

    @test rank(𝐉) == 1
    @test components(𝐉) == -1:1
    @test string(𝐉) == "𝐉̂⁽¹⁾"

    @test_throws ArgumentError TensorComponent(𝐂⁵, 7)
    𝐂⁵₃ = TensorComponent(𝐂⁵, 3)
    @test string(𝐂⁵₃) == "𝐂̂⁽⁵⁾₃"

    @test parent(𝐂⁵₃) == 𝐂⁵
    @test component(𝐂⁵₃) == 3
    @test system(𝐂⁵₃) == system(𝐂⁵) == OrbitalAngularMomentumSubSystem()

    @test 𝐂⁵₃ == TensorComponent(𝐂⁵, 3)
    @test hash(𝐂⁵₃) == hash(TensorComponent(𝐂⁵, 3))

    X = TensorProduct(4, 𝐂⁵, 𝐉)
    @test system(X) == (system(𝐂⁵), system(𝐉))
    @test rank(X) == 4
    @test string(X) == "{𝐂̂⁽⁵⁾𝐉̂⁽¹⁾}⁽⁴⁾"

    @test_throws ArgumentError dot(𝐂⁵, 𝐉)

    CJ = dot(SphericalTensor(1),𝐉)
    @test CJ isa TensorScalarProduct
    @test rank(CJ) == 0
    @test system(CJ) == (OrbitalAngularMomentumSubSystem(), TotalAngularMomentumSubSystem())
    @test string(CJ) == "(𝐂̂⁽¹⁾⋅𝐉̂⁽¹⁾)"

    @testset "OrbitalRadialOverlap" begin
        oro = OrbitalRadialOverlap(so"1s₀α", so"1s₀β")
        @test radial_integral(OrbitalOverlap(so"1s₀α", so"1s₀β")) == oro
        @test !iszero(oro)
        @test oro == OrbitalRadialOverlap(so"1s₀α", so"1s₀β")
        @test oro' == OrbitalRadialOverlap(so"1s₀β", so"1s₀α")
        @test hash(oro') == hash(OrbitalRadialOverlap(so"1s₀β", so"1s₀α"))
        @test numbodies(oro) == 0

        @test isdependent(oro, conj(so"1s₀α"))
        @test !isdependent(oro, conj(so"1s₀β"))
        @test !isdependent(oro, so"1s₀α")
        @test isdependent(oro, so"1s₀β")

        @test string(oro) == "⟨1s₀α|1s₀β⟩ᵣ"
    end

    @testset "TensorOperator" begin
        H_cfgs = spin_configurations([c"1s"])
        He_cfgs = spin_configurations([rc"1s2"])

        Tz = T -> TensorOperator{1}(TensorComponent(T, 0))
        oro = (args...) -> NBodyMatrixElement([NBodyTerm([OrbitalRadialOverlap(a, b) for (a,b) in args], 1)])
        orm = (a,o,b) -> NBodyMatrixElement([NBodyTerm([radial_integral(a, o, b)], 1)])

        a = so"1s₀α"
        b = so"1s₀β"
        c = so"2p₀β"
        a′ = so"ks₀α"
        b′ = so"ls₀β"
        ra = rso"1s(1/2)"
        rb = rso"1s(-1/2)"
        ra′ = rso"ks(1/2)"
        rb′ = rso"ls(-1/2)"

        A = Tz(𝐉)
        B = many_electron_scalar_product(𝐉)
        C = many_electron_scalar_product(OrbitalAngularMomentum(), SpinAngularMomentum())

        @test string(A) == "[𝐉̂⁽¹⁾₀]"
        @test string(B) == "[𝐉̂⁽¹⁾]² + 2.0[𝐉̂⁽¹⁾(1)⋅𝐉̂⁽¹⁾(2)]"
        @test string(C) == "[𝐋̂⁽¹⁾⋅𝐒̂⁽¹⁾] + 2.0[𝐋̂⁽¹⁾(1)⋅𝐒̂⁽¹⁾(2)]"

        @test all(Matrix(A, H_cfgs) .≈ 1/2*NBodyMatrixElement[oro((a,a)) 0.0; 0.0 -oro((b,b))])
        @test all(Matrix(A, He_cfgs) .≈ 1/2*NBodyMatrixElement[oro((ra,ra))-oro((rb,rb))])

        @test all(Matrix(B, H_cfgs) .≈ 1/2*(1/2+1)*NBodyMatrixElement[oro((a,a)) 0.0; 0.0 oro((b,b))])
        @test all(Matrix(B, He_cfgs) .≈ [1/2*(1/2+1)*(oro((ra,ra))+oro((rb,rb))) - oro((rb,ra),(ra,rb)) - 1/2*oro((rb,rb),(ra,ra))])

        @test iszero(Matrix(C, H_cfgs))
        @test iszero(Matrix(C, He_cfgs))

        @test integrate_spinor(OrbitalMatrixElement([so"ks₀α"], A, [so"ks₀α"])) ≈ 1/2*oro((so"ks₀α", so"ks₀α"))
        @test integrate_spinor(OrbitalMatrixElement([so"kp₁β"], A, [so"kp₁β"])) ≈ 1/2*oro((so"kp₁β", so"kp₁β"))
        @test integrate_spinor(OrbitalMatrixElement([so"kp₁α"], A, [so"kp₁α"])) ≈ 3/2*oro((so"kp₁α", so"kp₁α"))

        function generate_scalar_product_omes(left, 𝐓, 𝐔, right)
            O = TensorOperator{2}(𝐓⋅𝐔)

            omes = [OrbitalMatrixElement(l, O, r)
                    for l in left, r in right]
            integrate_spinor.(omes)
        end

        for (a,b,a′,b′) = [(a,b,a′,b′),(ra,rb,ra′,rb′)]
            left_products = [[a,a′],[a′,b′],[b′,a′],[b,b′]]
            right_products = [[a,a′],[a,b],[b,a],[b,b′]]

            int_omes = generate_scalar_product_omes(left_products, 𝐒, 𝐒, right_products)
            ref = 1/4*NBodyMatrixElement[oro((a,a),(a′,a′)) 0 0 0
                                         0 -oro((a′,a),(b′,b)) 2oro((a′,b),(b′,a)) 0
                                         0 2oro((b′,a),(a′,b)) -oro((b′,b),(a′,a)) 0
                                         0 0 0 oro((b,b),(b′,b′))]

            @test all(int_omes .≈ ref)
        end


        for (c,(a,b,a′,b′)) = [(√3,(so"2p₋₁α",so"kd₀α",so"2p₀α",so"ld₋₁α")),
                               (√3,(so"2p₋₁α",so"kd₁β",so"2p₀α",so"ld₀β")),
                               (2,(so"kd₂β",so"2p₁α",so"ld₂β",so"2p₁α"))]
            left_products = [[a′,b′],[b′,a′]]
            right_products = [[a,b],[b,a]]

            int_omes = generate_scalar_product_omes(left_products, 𝐋, 𝐋, right_products)
            ref = c*NBodyMatrixElement[oro((a′,a),(b′,b)) 0
                                       0 oro((b′,b),(a′,a))]
            @test all(int_omes .≈ ref)
        end

        @testset "OrbitalRadialMatrixElement" begin
            cfgs = [sc"1s₀α 1s₀β", sc"1s₀α 2p₀β"]

            r̂ = RadialOperator()
            ∂̂ᵣ = k -> RadialGradientOperator(k)

            @test radial_integral([b],r̂,[c]) == OrbitalRadialMatrixElement([b],r̂,[c])
            @test radial_integral([b],2.0,[c]) == 2oro((b,c))

            @test string(orm([b],r̂,[c])) == "⟨1s₀β|r|2p₀β⟩ᵣ"

            val = Matrix(Tz(Dipole()), cfgs)
            ref = 1/√3*NBodyMatrixElement[0 orm([b],r̂,[c])
                                          orm([c],r̂,[b]) 0]
            @test all(val .≈ ref)

            val = Matrix(Tz(Gradient()), cfgs)
            ref = 1/√3*NBodyMatrixElement[0 orm([b],∂̂ᵣ(2),[c])
                                          orm([c],∂̂ᵣ(0),[b]) 0]
            @test all(val .≈ ref)

            val = Matrix(Tz(SphericalTensor(1)), cfgs)
            ref = 1/√3*NBodyMatrixElement[0 oro((b,c))
                                          oro((c,b)) 0]
            @test all(val .≈ ref)
        end
    end
end
