using AngularMomentumAlgebra.Dipoles

@testset "Wigner–Eckart" begin
    @testset "Dipoles" begin
        x,y,z = 𝐫̂

        os = os"k[s-h]"

        @testset "z" begin
            for o in os
                ℓ = o.ℓ
                ms = -ℓ:ℓ
                for m in ms
                    @test wigner_eckart(ℓ,m,z,ℓ+1,m) ≈ √((ℓ+m+1)*(ℓ-m+1)/((2ℓ+3)*(2ℓ+1)))
                end
            end
        end

        @testset "x" begin
            for o in os
                ℓ = o.ℓ
                ms = -ℓ:ℓ
                for m in ms
                    # TODO Figure out where the minus sign comes from
                    @test wigner_eckart(ℓ,m,x,ℓ+1,m+1) ≈ -0.5*√((ℓ+m+2)*(ℓ+m+1)/((2ℓ+3)*(2ℓ+1)))
                end
            end
        end
    end
end
