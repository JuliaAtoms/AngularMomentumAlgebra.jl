using Unicode
using LightGraphs
using TikzGraphs
using LaTeXStrings

preamble = read(joinpath(dirname(@__FILE__), "..", "src", "preamble.tex"), String)

# * Angular roots

struct AngularRoot{T<:Number} <: Symbolic
    v::T
end

@new_number AngularRoot
Base.:(==)(x::AngularRoot, y::AngularRoot) = x.v == y.v

Base.convert(::Type{T}, root::AngularRoot) where {T<:Union{Real,Complex}} =
    √(2convert(T, root.v) + 1)

Base.convert(::Type{Sym}, root::AngularRoot{Sym}) =
    √(2root.v + 1)

Base.show(io::IO, r::AngularRoot) =
    write(io, Unicode.normalize("$(r.v)̂"))
latex(r::AngularRoot) = "\\hat{$(latex(r.v))}"

# * Sum variables
struct SumVariable <: Symbolic
    v::Symbol
end

@new_number SumVariable
Base.:(==)(a::SumVariable, b::SumVariable) = a.v == b.v

Base.show(io::IO, s::SumVariable) = write(io, "Σ$(s.v)")
latex(s::SumVariable) = "$(s.v)_\\Sigma"

const sum_variables = "αβγδεζηθικλμνξοπρστυφχψωΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ"

function get_sum_variables(n::Int, i=0)
    j = firstindex(sum_variables)
    for j′ = 1:i
        j = nextind(sum_variables, j)
    end
    map(1:n) do j′
        s = sum_variables[j]
        j = nextind(sum_variables, j)
        SumVariable(Symbol(s))
    end
end

# * Comparators
@gen_compare_false AngularRoot SumVariable IIIJ VIJ IXJ
# @new_numbers AngularRoot IIIJ VIJ IXJ

# * Couplings

mutable struct Coupling{T<:Union{HalfInteger,Symbolic}, A<:Number}
    J::T
    amplitude::A
end

Coupling{S,A}(J::S) where {S,A} = Coupling{S,A}(J, one(A))
Coupling(J::HalfInteger) = Coupling{HalfInteger,Float64}(J)
Coupling(J::Symbolic) = Coupling{Symbolic,Number}(J)
Coupling(J::Symbol) = Coupling{Symbolic,Number}(Sym(J))
Coupling(J::S,a::A) where {S<:Symbolic,A<:Number} = Coupling{Symbolic,Number}(J,a)

Base.show(io::IO, coupling::Coupling) = show(io, coupling.J)

latex(coupling::Coupling{HalfInteger,Int}) = latexstring("$(coupling.J)")

function latex(coupling::Coupling{Symbolic,Number})
    S = "\$$(latex(coupling.J))\$"
    if coupling.amplitude != 1
        S *= "\\\\\\\\\$$(latex(coupling.amplitude))\$"
    end
    latexstring(S)
end

# function Base.show(io::IO, coupling::Coupling)
#     write(io, "[")
#     show(io, coupling.j₁)
#     write(io, " ")
#     show(io, coupling.j₂)
#     write(io, "]^{")
#     show(io, coupling.J)
#     write(io, "}")
# end

# * Coupling trees

mutable struct CouplingTree{C<:Coupling}
    couplings::Vector{C}
    tree::SimpleDiGraph{Int}
    sum_variables::Vector{SumVariable}
end

CouplingTree{C}(couplings, tree) where C =
    CouplingTree{C}(couplings, tree, Vector{SumVariable}())

function Base.replace!(tree::CouplingTree, p::Pair{Int,<:Coupling})
    tree.couplings[p[1]] = p[2]
    tree
end

layers(tree::CouplingTree) = Int(log2(length(tree.couplings)+1))

momenta(tree::CouplingTree, nodes::Int...) =
    map(nodes) do node
        tree.couplings[node].J
    end

amplitudes(tree::CouplingTree, nodes::Int...) =
    map(nodes) do node
        tree.couplings[node].amplitude
    end

function directed_binary_tree(k::I) where I
    g = BinaryTree(k)
    g′ = SimpleDiGraph{I}(length(vertices(g)))
    for v in vertices(g)
        for n in neighbors(g,v)
            n > v && add_edge!(g′, v, n)
        end
    end
    g′
end

function CouplingTree(couplings::Vector{C}) where C
    # S = 2ⁿ - 1 => n = log2(S+1)
    n = log2(length(couplings)+1)
    isinteger(n) || throw(ArgumentError("Not a balanced binary tree"))
    n = Int(n)
    tree = directed_binary_tree(n)
    CouplingTree{C}(couplings, tree)
end

function swap_couplings!(tree::CouplingTree, a::Int, b::Int)
    tmp = tree.couplings[a]
    tree.couplings[a] = tree.couplings[b]
    tree.couplings[b] = tmp
end

function swap_subtrees!(tree::CouplingTree, a::Int, b::Int)
    A = dfs_tree(tree.tree, a)
    B = dfs_tree(tree.tree, b)
    @assert length(edges(A)) == length(edges(B))
    Astack = Int[a]
    Bstack = Int[b]
    while !isempty(Astack)
        a,b = pop!(Astack),pop!(Bstack)
        swap_couplings!(tree, a, b)
        append!(Astack, neighbors(tree.tree,a))
        append!(Bstack, neighbors(tree.tree,b))
    end
end

# Eq. (21) Williams 1992
function switch!(tree::CouplingTree, c::Int)
    nn = neighbors(tree.tree, c)
    isempty(nn) && return tree
    a,b = nn
    swap_subtrees!(tree, a, b)
    tree.couplings[c].amplitude *= (-1)^(tree.couplings[a].J
                                         + tree.couplings[b].J
                                         - tree.couplings[c].J)
    tree
end

# Eq. (23) Williams 1992
function switch_9j!(tree::CouplingTree, z::Int)
    z ≤ 2^(layers(tree)-2) - 1 ||
        throw(ArgumentError("Cannot perform a 9j switch starting from $(z) (too close to root)"))

    x,y = neighbors(tree.tree, z)
    a,b = neighbors(tree.tree, x)
    c,d = neighbors(tree.tree, y)

    ã,b̃,c̃,d̃,x̃,ỹ,z̃ = momenta(tree, a, b, c, d, x, y, z)

    p,q = get_sum_variables(2, length(tree.sum_variables))
    append!(tree.sum_variables, [p,q])

    x̃̂,ỹ̂,p̂,q̂ = AngularRoot.((x̃,ỹ,p,q))
    ixj = IXJ(ã,b̃,x̃,
              c̃,d̃,ỹ,
              p,q,z̃)
    tree.couplings[z].amplitude *= x̃̂*ỹ̂*p̂*q̂*ixj

    swap_subtrees!(tree, b, c)

    X,Y = amplitudes(tree, x, y)
    replace!(tree, x => Coupling(p, X))
    replace!(tree, y => Coupling(q, Y))

    tree
end

TikzGraphs.plot(tree::CouplingTree; kwargs...) =
    plot(tree.tree, latex.(tree.couplings);
         node_style="draw, rounded corners, align=center",
         prepend_preamble=preamble,
         kwargs...)


export AngularRoot, SumVariable, Coupling, CouplingTree, switch!, switch_9j!
