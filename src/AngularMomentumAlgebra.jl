module AngularMomentumAlgebra

using AtomicLevels
using UnicodeFun
using LinearAlgebra
import AtomicLevels: HalfInteger
using Parameters

include("symbolics.jl")
include("j.jl")
include("common.jl")
include("clebsch_gordan.jl")
include("kronecker.jl")
include("tensors.jl")
include("slater_integrals.jl")
include("couplings.jl")

end # module
