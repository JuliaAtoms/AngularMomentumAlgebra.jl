@testset "Tensors" begin
    import AngularMomentumAlgebra: components, component,
    OrbitalRadialOverlap
    @test system(Tensor) == FullSystem()

    𝐂⁵ = SphericalTensor(5)
    𝐉 = TotalAngularMomentum()

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
end
