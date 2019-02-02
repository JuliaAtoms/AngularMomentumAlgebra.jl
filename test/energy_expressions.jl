@testset "Energy expressions" begin
    c1 = spin_configurations(c"1s2")
    c2 = spin_configurations(c"1s 2s")
    O = eltype(c1[1])[1]

    @testset "One-body energy" begin
        h1 = one_body_hamiltonian_matrix(O,c1)
        @test size(h1) == (1,1)
        # TODO: Proper tests of matrix elements
        h2 = one_body_hamiltonian_matrix(O,c2)
        @test size(h2) == (4,4)
        # TODO: Proper tests of matrix elements
    end
    
    @testset "Two-body energy" begin
        HC1 = two_body_hamiltonian_matrix(O,c1)
        @test size(HC1) == (1,1)
        # TODO: Proper tests of matrix elements
        HC2 = two_body_hamiltonian_matrix(O,c2)
        @test size(HC2) == (4,4)
        # TODO: Proper tests of matrix elements
    end
end
