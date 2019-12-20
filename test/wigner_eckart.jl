using AngularMomentumAlgebra.Dipoles

@testset "Wigner–Eckart" begin
    @testset "Dipoles" begin
        x,y,z = 𝐫̂

        os = os"k[s-h]"

        @testset "z" begin
            for o in os
                ℓ = o.ℓ
                mℓs = -ℓ:ℓ
                for mℓ in mℓs
                    @test wigner_eckart((ℓ,mℓ),z,(ℓ+1,mℓ)) ≈ √((ℓ+mℓ+1)*(ℓ-mℓ+1)/((2ℓ+3)*(2ℓ+1)))
                end
            end
        end

        @testset "x" begin
            for o in os
                ℓ = o.ℓ
                mℓs = -ℓ:ℓ
                for mℓ in mℓs
                    # TODO Figure out where the minus sign comes from
                    @test wigner_eckart((ℓ,mℓ),x,(ℓ+1,mℓ+1)) ≈ -0.5*√((ℓ+mℓ+2)*(ℓ+mℓ+1)/((2ℓ+3)*(2ℓ+1)))
                end
            end
        end
    end
end
