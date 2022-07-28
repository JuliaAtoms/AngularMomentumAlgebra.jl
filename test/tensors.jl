@testset "Tensors" begin
    import AngularMomentumAlgebra: components, component,
    OrbitalRadialOverlap, radial_integral
    import EnergyExpressions: NBodyTerm, NBodyMatrixElement
    @test system(Tensor) == FullSystem()

    𝐂⁵ = SphericalTensor(5)
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
        oro = (args...) -> NBodyMatrixElement([NBodyTerm([OrbitalRadialOverlap(a, b) for (a,b) in args], 1)])

        a = so"1s₀α"
        b = so"1s₀β"
        ra = rso"1s(1/2)"
        rb = rso"1s(-1/2)"

        A = TensorOperator{1}(𝐉₀)
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
    end
end
