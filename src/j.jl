using WignerSymbols
using Formatting

# * Pretty-printing

function show_singleline(io::IO, s::String,e::String,
                         lines::Tuple...)
    n = 0

    for line in lines
        for v in last(line)
            n = max(n, length(string(v)))
        end
    end
    f = FormatExpr(join(["{$i:$(n)s}" for i in 1:length(first(lines))], " "))
    write(io, s, " ")
    for (i,line) in enumerate(lines)
        write(io, format(f, string.(line)...))
        i != lastindex(lines) && write(io, " ; ")
    end
    write(io, " ", e)
end

function show_multiline(io::IO, lines::Pair{<:NTuple{2,String},<:Tuple}...)
    n = 0

    for line in lines
        for v in last(line)
            n = max(n, length(string(v)))
        end
    end
    f = FormatExpr(join(["{$i:$(n)s}" for i in 1:length(last(first(lines)))], " "))
    for (i,line) in enumerate(lines)
        brackets = first(line)
        write(io, brackets[1], format(f, string.(last(line))...), brackets[2])
        i != lastindex(lines) && write(io, "\n")
    end
end

# * 3j

struct IIIJ{A,B,C,D,E,F}
    j₁::A
    j₂::B
    j₃::C
    m₁::D
    m₂::E
    m₃::F
end

Base.convert(::Type{T}, iiij::IIIJ{R,R,R,R,R,R}) where {T<:Number,R<:Real} =
    wigner3j(T,
             iiij.j₁, iiij.j₂, iiij.j₃,
             iiij.m₁, iiij.m₂, iiij.m₃)

Base.show(io::IO, iiij::IIIJ) =
    show_singleline(io, "(", ")",
                    (iiij.j₁, iiij.j₂, iiij.j₃),
                    (iiij.m₁, iiij.m₂, iiij.m₃))

Base.show(io::IO, ::MIME"text/plain", iiij::IIIJ) =
    show_multiline(io,
                   ("⎛","⎞")=>(iiij.j₁, iiij.j₂, iiij.j₃),
                   ("⎝","⎠")=>(iiij.m₁, iiij.m₂, iiij.m₃))

function triangle_inequality(iiij::IIIJ)
    println("|$(iiij.j₁)-$(iiij.j₂)| ≤ $(iiij.j₃) ≤ $(iiij.j₁)+$(iiij.j₂)")
end

# * 6j

struct VIJ{A,B,C,D,E,F}
    j₁::A
    j₂::B
    j₃::C
    j₄::D
    j₅::E
    j₆::F
end

Base.convert(::Type{T}, vij::VIJ{R,R,R,R,R,R}) where {T<:Number,R<:Real} =
    wigner6j(T,
             vij.j₁, vij.j₂, vij.j₃,
             vij.j₄, vij.j₅, vij.j₆)

Base.show(io::IO, vij::VIJ) =
    show_singleline(io, "{", "}",
                    (vij.j₁, vij.j₂, vij.j₃),
                    (vij.j₄, vij.j₅, vij.j₆))

Base.show(io::IO, ::MIME"text/plain", vij::VIJ) =
    show_multiline(io,
                   ("⎰","⎱")=>(vij.j₁, vij.j₂, vij.j₃),
                   ("⎱","⎰")=>(vij.j₄, vij.j₅, vij.j₆))

# * 9j

struct IXJ{A,B,C,D,E,F,G,H,I}
    j₁::A
    j₂::B
    j₃::C
    j₄::D
    j₅::E
    j₆::F
    j₇::G
    j₈::H
    j₉::I
end

# # wigner9j not yet implemented
# Base.convert(::Type{T}, ixj::IXJ{R,R,R,R,R,R}) where {T<:Number,R<:Real} =
#     wigner9j(T,
#              ixj.j₁, ixj.j₂, ixj.j₃,
#              ixj.j₄, ixj.j₅, ixj.j₆,
#              ixj.j₇, ixj.j₈, ixj.j₉)

Base.show(io::IO, ixj::IXJ) =
    show_singleline(io, "{", "}",
                    (ixj.j₁, ixj.j₂, ixj.j₃),
                    (ixj.j₄, ixj.j₅, ixj.j₆),
                    (ixj.j₇, ixj.j₈, ixj.j₉))

Base.show(io::IO, ::MIME"text/plain", ixj::IXJ) =
    show_multiline(io,
                   ("⎧","⎫")=>(ixj.j₁, ixj.j₂, ixj.j₃),
                   ("⎨","⎬")=>(ixj.j₄, ixj.j₅, ixj.j₆),
                   ("⎩","⎭")=>(ixj.j₇, ixj.j₈, ixj.j₉))

export IIIJ, VIJ, IXJ
