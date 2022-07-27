@testset "Tensor DSL" begin
    import AngularMomentumAlgebra: wrap_show_debug, remove_line_numbers,
    generate_signature,
    strip′, identify_quantum_numbers, recursepm,
    parse_selection_rule, push_selection_rules!,
    generate_iszero, generate_rme, generate_couplings,
    @tensor

    @testset "Macro utilities" begin
        @test wrap_show_debug(:(sin(a))) == :(sin(a))
        ENV["DEBUG_TENSOR_DSL"] = "1"
        @test remove_line_numbers(wrap_show_debug(:(sin(a)))) ==
            remove_line_numbers(:(@show sin(a)))
        delete!(ENV, "DEBUG_TENSOR_DSL")
        @test wrap_show_debug(:(sin(a))) == :(sin(a))

        @test generate_signature(t -> (:ℓ′, t, :ℓ), :(Base.iszero), :TensorType) ==
            :(Base.iszero(ℓ′, tensor::TensorType, ℓ))
        @test generate_signature(t -> (:ℓ′, t, :ℓ), :(Base.iszero), :(SphericalTensor{k} where k)) ==
            :(Base.iszero(ℓ′, tensor::SphericalTensor{k}, ℓ) where k)
        @test generate_signature(t -> (t, :ℓ), :couplings, :TensorType) ==
            :(couplings(tensor::TensorType, ℓ))

        @test_throws ArgumentError strip′(:(ℓ))
        @test strip′(:(ℓ′)) == :(ℓ)


        sr = :(ℓ′ == ℓ)
        @test identify_quantum_numbers(sr) == (:ℓ′, :ℓ)
        bsr = quote
            n′ ~ n # The dipole couples orbitals of different n, but
            # there is no selection rule.
            ℓ′ == ℓ ± 1
        end
        @test identify_quantum_numbers(bsr) == (:((n′,ℓ′)),:((n,ℓ)))

        @test recursepm(:ℓ) == [:ℓ]
        @test recursepm(:(ℓ ± 1)) == [:(ℓ + 1), :(ℓ - 1)]
        @test recursepm(:(ℓ ∓ 1)) == [:(ℓ - 1), :(ℓ + 1)]
        @test recursepm(:(ℓ ± 1 ∓ n)) == [:((ℓ + 1) - n),
                                          :((ℓ - 1) - n),
                                          :((ℓ + 1) + n),
                                          :((ℓ - 1) + n)]
        @test recursepm(:(ℓ + 1)) == [:(ℓ + 1)]
        @test recursepm(:(ℓ + (1 ∓ n))) == [:(ℓ + (1 - n)), :(ℓ + (1 + n))]

        @test recursepm(:(ℓ′), :(ℓ)) == :(ℓ′ == ℓ)
        @test recursepm(:(ℓ′), :(ℓ ± 1)) == :(ℓ′ == ℓ - 1 || ℓ′ == ℓ + 1)
        @test recursepm(:(ℓ′), :(ℓ ± 1 ∓ n)) == :(ℓ′ == (ℓ + 1) - n ||
            ℓ′ == (ℓ - 1) - n ||
            ℓ′ == (ℓ - 1) + n ||
            ℓ′ == (ℓ + 1) + n)

        @test parse_selection_rule(:(ℓ′ == ℓ)) ==
            :(ℓ′ == ℓ || return true)
        @test parse_selection_rule(:(ℓ′ == ℓ ± 1)) ==
            :((ℓ′ == ℓ - 1 || ℓ′ == ℓ + 1) || return true)
        @test parse_selection_rule(:(ℓ′ ∈ abs(ℓ - k):2:(ℓ+k))) ==
            :(ℓ′ ∈ abs(ℓ - k):2:(ℓ+k) || return true)

        srbody = Expr(:block)
        push_selection_rules!(parse_selection_rule, srbody.args, sr)
        @test srbody.args == [:(ℓ′ == ℓ || return true)]

        bsrbody = Expr(:block)
        push_selection_rules!(parse_selection_rule, bsrbody.args, bsr)
        @test remove_line_numbers(bsrbody.args) == [:((ℓ′ == ℓ - 1 || ℓ′ == ℓ + 1) || return true)]
    end

    function function_name(fn)
        if fn.head == :call
            first(fn.args)
        elseif fn.head == :where
            first(first(fn.args).args)
        else
            error("Do not understand this function: $(fn)")
        end
    end

    function function_definition_name(f)
        @test f.head == :function
        function_name(first(f.args))
    end

    @testset "Function generators" begin
        sr = :(ℓ′ == ℓ)
        bsr = quote
            n′ ~ n # The dipole couples orbitals of different n, but
            # there is no selection rule.
            ℓ′ == ℓ ± 1
        end

        lnn = LineNumberNode(123, :here)

        @testset "generate_iszero" begin
            f = generate_iszero(:OrbitalAngularMomentum, sr, lnn)
            @test function_definition_name(f) == :(Base.iszero)
            @test f.args[1].args == [:(Base.iszero), :ℓ′, :(tensor::OrbitalAngularMomentum), :ℓ]
            @test f.args[2].args == [lnn, :(ℓ′ == ℓ || return true), false]

            g = generate_iszero(:Dipole, bsr, lnn)
            @test function_definition_name(g) == :(Base.iszero)
            @test g.args[1].args == [:(Base.iszero), :((n′, ℓ′)), :(tensor::Dipole), :((n,ℓ))]
            @test g.args[2].args[1] == lnn
            @test remove_line_numbers(g.args[2]).args[1] == [:((ℓ′ == ℓ - 1 || ℓ′ == ℓ + 1) || return true), false]
        end

        @testset "generate_rme" begin
            f = generate_rme(:OrbitalAngularMomentum, sr, "The Docs", :(sin(a)), lnn)
            args = remove_line_numbers(f.args)
            @test first(args) == Symbol("@doc")
            @test args[2] == "The Docs"
            ff = args[3]
            @test function_definition_name(ff) == :rme
            @test ff.args[1].args == [:rme, :ℓ′, :(tensor::OrbitalAngularMomentum), :ℓ]
            @test ff.args[2].args == [:(iszero(ℓ′, tensor, ℓ) && return 0), lnn, :(sin(a))]

            g = generate_rme(:Dipole, bsr, "The Docs", :(sin(a)), lnn)
            args = remove_line_numbers(g.args)
            @test first(args) == Symbol("@doc")
            @test args[2] == "The Docs"
            gg = args[3]
            @test function_definition_name(gg) == :rme
            @test gg.args[1].args == [:rme, :(n′,ℓ′), :(tensor::Dipole), :(n,ℓ)]
            @test gg.args[2].args == [:(iszero((n′, ℓ′), tensor, (n, ℓ)) && return 0), lnn, :(sin(a))]
        end

        @testset "generate_couplings" begin
            f = generate_couplings(:OrbitalAngularMomentum, sr, lnn)
            args = remove_line_numbers(f.args)
            @test first(args) == Symbol("@doc")
            @test occursin("couplings", args[2])
            @test occursin("Generate all quantum numbers", args[2])
            ff = args[3]
            @test function_definition_name(ff) == :couplings
            @test ff.args[1].args == [:couplings, :(tensor::OrbitalAngularMomentum), :ℓ]
            @test ff.args[2].args == [lnn, :(ℓ′ = [ℓ]), :ℓ′]

            g = generate_couplings(:Dipole, bsr, lnn)
            args = remove_line_numbers(g.args)
            @test first(args) == Symbol("@doc")
            @test occursin("couplings", args[2])
            @test occursin("Generate all quantum numbers", args[2])
            gg = args[3]
            @test function_definition_name(gg) == :couplings
            @test gg.args[1].args == [:couplings, :(tensor::Dipole), :(n,ℓ)]
            @test gg.args[2].args[1] == lnn
            @test function_name(gg.args[2].args[2]) == :error

            # Does not end in prime
            @test_throws ArgumentError generate_couplings(:MyTensor, :(ℓ == 1), lnn)
            # Does not know how to handle a constant
            @test_throws ArgumentError generate_couplings(:MyTensor, :(ℓ′ == 1), lnn)
            # Unequalities not supported
            @test_throws ArgumentError generate_couplings(:MyTensor, :(ℓ′ < 4), lnn)
        end
    end

    @testset "@tensor macro" begin
        @testset "SphericalTensor" begin
            met = @macroexpand begin
                @tensor(SphericalTensor{k} where k) do
                    ℓ′ ∈ abs(ℓ - k):2:(ℓ+k)

                    raw"""
                    rme(ℓ′,𝐂̂ᵏ,ℓ)

                Calculate the reduced matrix element of the spherical tensor of rank
                `k`:

                ```math
                \begin{aligned}
                \redmatrixel{\ell'}{\tensor{C}^{(k)}}{\ell}
                &=
                \angroot{\ell}
                C_{\ell 0;k,0}^{\ell'0} =
                (-)^{\ell-k}
                \angroot{\ell\ell'}
                \wignerthreej{\ell&k&\ell'\\0&0&0}.
                \end{aligned}
                \tag{V13.2.107}
                ```
                """
                    ∏(ℓ)*clebschgordan(ℓ,0,k,0,ℓ′,0)
                end
            end
            generated_functions = remove_line_numbers(met.args)[1].args
            @test length(generated_functions) == 3
            fa,fb,fc = generated_functions

            @test function_definition_name(fa) == :(Base.iszero)

            @test fb.head == :block
            @test function_definition_name(first(fb.args)) == :rme
            @test first(fb.args[2].args) == Base.Docs.doc!

            @test fc.head == :block
            @test function_definition_name(first(fc.args)) == :couplings
            @test first(fc.args[2].args) == Base.Docs.doc!
        end

        @testset "Gradient" begin
            met = @macroexpand begin
                @tensor(Gradient) do
                    begin
                        n′ ~ n # The gradient couples orbitals of different n, but
                        # there is no selection rule.
                        ℓ′ == ℓ ± 1
                    end

                    raw"""
                   rme((n′,ℓ′), ::Gradient, (n,ℓ))

               Computes the reduced matrix element of `∇` in terms of
               [`RadialGradientMatrixElement`](@ref).
               """
                    if ℓ′ == ℓ+1
                        √(ℓ+1)*RadialGradientMatrixElement(-ℓ)
                    elseif ℓ′==ℓ-1
                        - √(ℓ)*RadialGradientMatrixElement(ℓ+1)
                    end
                end
            end

            generated_functions = remove_line_numbers(met.args)[1].args
            @test length(generated_functions) == 3
            fa,fb,fc = generated_functions

            @test function_definition_name(fa) == :(Base.iszero)

            @test fb.head == :block
            @test function_definition_name(first(fb.args)) == :rme
            @test first(fb.args[2].args) == Base.Docs.doc!

            @test fc.head == :block
            @test function_definition_name(first(fc.args)) == :couplings
            @test function_name(fc.args[1].args[2].args[2]) == :error
            @test first(fc.args[2].args) == Base.Docs.doc!
        end
    end
end
