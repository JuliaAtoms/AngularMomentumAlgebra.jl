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

Base.convert(::Type{T}, root::AngularRoot) where {T<:Union{Real,Complex}} =
    √(2convert(T, root.v) + 1)

Base.convert(::Type{Sym}, root::AngularRoot{Sym}) =
    √(2root.v + 1)

Base.show(io::IO, r::AngularRoot) =
    write(io, Unicode.normalize("$(r.v)̂"))
function Symbolics.latex(r::AngularRoot)
    lv,d = Symbolics.latex(r.v)
    "\\hat{$(lv)}",d
end

# * Sum variables
struct SumVariable <: Symbolic
    v::Symbol
end

@new_number SumVariable

Base.show(io::IO, s::SumVariable) = write(io, "Σ$(s.v)")
Symbolics.latex(s::SumVariable) = "$(s.v)_\\Sigma",0

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

# * Couplings

mutable struct Coupling{T<:Union{HalfInteger,Symbolic}, A<:Number}
    J::T
    amplitude::A
end

Coupling{S,A}(J::S) where {S,A} = Coupling{S,A}(J, one(A))
Coupling(J::HalfInteger) = Coupling{HalfInteger,Float64}(J)
Coupling(J::Integer) = Coupling(HalfInteger(J))
Coupling(J::Symbolic) = Coupling{Symbolic,Number}(J)
Coupling(J::Symbol) = Coupling{Symbolic,Number}(Sym(J))
Coupling(J::S,a::A) where {S<:Symbolic,A<:Number} = Coupling{Symbolic,Number}(J,a)

Base.zero(::Type{Coupling{S,A}}) where {S,A} = Coupling{S,A}(zero(S))

Base.:(==)(a::Coupling{S,A},b::Coupling{S,A}) where {S,A} =
    a.J == b.J && a.amplitude == b.amplitude

function Base.show(io::IO, coupling::Coupling)
    if coupling.amplitude != 1
        show(io, coupling.amplitude)
    end
    show(io, coupling.J)
end

Symbolics.latex(coupling::Coupling{HalfInteger,Int}) = latexstring("$(coupling.J)"),0

function Symbolics.latex(coupling::Coupling{Symbolic,Number})
    S = "\$$(latex(coupling.J))\$"
    if coupling.amplitude != 1
        S *= "\\\\\\\\\$$(latex(coupling.amplitude))\$"
    end
    latexstring(S),0
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

function directed_binary_tree(k::I) where {I<:Integer}
    g = BinaryTree(k)
    g′ = SimpleDiGraph{I}(length(vertices(g)))
    for v in vertices(g)
        for n in neighbors(g,v)
            n > v && add_edge!(g′, v, n)
        end
    end
    g′
end

function CouplingTree{C}(k::I) where {C<:Coupling,I<:Integer}
    tree = directed_binary_tree(k)
    couplings = zeros(C, length(vertices(tree)))
    CouplingTree{C}(couplings, tree)
end

function CouplingTree(couplings::Vector{C}) where {C<:Coupling}
    # S = 2ⁿ - 1 => n = log2(S+1)
    n = log2(length(couplings)+1)
    isinteger(n) || throw(ArgumentError("Not a balanced binary tree"))
    n = Int(n)
    tree = directed_binary_tree(n)
    CouplingTree{C}(couplings, tree)
end

CouplingTree(couplings::Vector{N}) where {N<:Number} =
    CouplingTree(Coupling.(couplings))

function Base.replace!(tree::CouplingTree, p::Pair{Int,<:Coupling})
    tree.couplings[p[1]] = p[2]
    tree
end

layers(tree::CouplingTree) = Int(log2(length(tree.couplings)+1))
layer(tree::CouplingTree, n::Int) = ceil(Int, log2(n+1))

final_momentum(tree::CouplingTree) = tree.couplings[1].J

function view_layer(tree::CouplingTree, i::Int)
    (i < 1 || i > layers(tree)) &&
        throw(BoundsError("Trying to access layer $(i) outside valid range"))
    a = 2^(i-1)
    b = 2^i - 1
    view(tree.couplings, a:b)
end

momenta(tree::CouplingTree, nodes::Int...) =
    map(nodes) do node
        tree.couplings[node].J
    end

amplitudes(tree::CouplingTree, nodes::Int...) =
    map(nodes) do node
        tree.couplings[node].amplitude
    end

function propagate_layer!(tree::CouplingTree{C}, i::Int) where C
    v = view_layer(tree, i)
    v′ = view_layer(tree, i+1)
    for i in eachindex(v)
        v′[2i-1] = v[i]
    end
end

function expand_layers!(tree::CouplingTree{C}, req_layers::Int) where C
    new_layers = layers(tree) - req_layers
    new_layers ≥ 0 && return tree
    tree.tree = directed_binary_tree(req_layers)
    append!(tree.couplings, zeros(C, length(vertices(tree.tree)) - length(tree.couplings)))
    num_layers = layers(tree)
    for i = (num_layers+new_layers):(num_layers-1)
        propagate_layer!(tree, i)
    end
    tree
end

function expand_tree_upwards(tree::CouplingTree{C}, req_layers::Int) where C
    new_layers = req_layers - layers(tree)
    new_layers ≤ 0 && return tree
    new_tree = CouplingTree{C}(req_layers)
    insert!(new_tree, 2^new_layers, tree)
    new_tree
end

function traverse(op::Function, tree::CouplingTree{C}, n::Int = 1) where C
    stack = Int[n]
    while !isempty(stack)
        n = pop!(stack)
        op(n) || break
        append!(stack, neighbors(tree.tree,n))
    end
end

valid_coupling(j₁::HalfInteger, j₂::HalfInteger, J::HalfInteger) =
    abs(j₁-j₂) ≤ J && J ≤ (j₁ + j₂)

function Base.isvalid(tree::CouplingTree{C}) where {C<:Coupling{HalfInteger,<:Number}}
    valid = true
    traverse(tree) do n
        nn = neighbors(tree.tree, n)
        if !isempty(nn)
            a,b = nn
            valid = valid_coupling(momenta(tree, a, b, n)...)
        else
            true
        end
    end
    valid
end

function modify_subtrees!(op::Function,
                          A::CouplingTree{C}, a::Integer,
                          B::CouplingTree{C}, b::Integer) where C
    Astack = Int[a]
    Bstack = Int[b]
    while !isempty(Astack) && !isempty(Bstack)
        a,b = pop!(Astack),pop!(Bstack)
        op((a,b))
        append!(Astack, neighbors(A.tree,a))
        append!(Bstack, neighbors(B.tree,b))
    end
end

function Base.insert!(tree::CouplingTree{C}, n::Integer, sub_tree::CouplingTree{C}) where C
    apex = sub_tree.couplings[1]
    insertion_point = tree.couplings[n]

    apex == insertion_point ||
        throw(ArgumentError("Apex node [$(apex)] does not match with insertion point [$(insertion_point)]"))

    # If trying to insert a too small subtree, i.e. at a node that is
    # high enough such that the sub-tree does not reach the root row
    # or beyond, it is an ill-defined operation, since we do not know
    # what the rows below should be filled with.
    layers(sub_tree) + layer(tree,n) ≤ layers(tree) &&
        throw(ArgumentError("Cannot insert subtree with $(layers(sub_tree)) layer(s) at node $(n) at layer $(layer(tree,n)) of $(layers(tree)); ill-defined operation"))

    expand_layers!(tree, layer(tree,n) + layers(sub_tree) - 1)

    modify_subtrees!(tree, n, sub_tree, 1) do (a,b)
        tree.couplings[a] = sub_tree.couplings[b]
    end
    tree
end

function swap_couplings!(tree::CouplingTree, a::Int, b::Int)
    tmp = tree.couplings[a]
    tree.couplings[a] = tree.couplings[b]
    tree.couplings[b] = tmp
end

function swap_subtrees!(tree::CouplingTree, a::Int, b::Int)
    A = dfs_tree(A.tree, a)
    B = dfs_tree(A.tree, b)
    @assert length(edges(A)) == length(edges(B))
    modify_subtrees!(tree, a, tree, b) do ab
        a,b = ab
        swap_couplings!(tree, a, b)
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
    plot(tree.tree, Symbolics.latex.(tree.couplings);
         node_style="draw, rounded corners, align=center",
         prepend_preamble=preamble,
         kwargs...)

function UnicodeFun.to_superscript(h::HalfInteger)
    if isinteger(h)
        to_superscript(convert(Int, h))
    else
        to_superscript(convert(Int, 2h)) * "⁽²⁾"
    end
end

function Base.show(io::IO, tree::CouplingTree{C}, n::Int=1) where {C<:Coupling{HalfInteger,<:Number}}
    nn = neighbors(tree.tree, n)
    if isempty(nn)
        show(io, tree.couplings[n])
    else
        amplitude = tree.couplings[n].amplitude
        amplitude != 1 && show(io, tree.couplings[n].amplitude)
        write(io, "[")
        show(io, tree, nn[1])
        write(io, " ")
        show(io, tree, nn[2])
        write(io, "]")
        write(io, to_superscript(tree.couplings[n].J))
    end
end


export AngularRoot, SumVariable, Coupling,
    CouplingTree, expand_tree_upwards, final_momentum,
    switch!, switch_9j!
