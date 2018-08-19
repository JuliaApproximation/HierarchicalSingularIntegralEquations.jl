
__precompile__()
module HierarchSingularIntegralEquations
    using Base, BandedMatrices, ApproxFun, SingularIntegralEquations, DualNumbers, RecipesBase,
            LinearAlgebra, Random, SpecialFunctions, LowRankApprox, InteractiveUtils

import Base: values, getindex, setindex!, *, +, -, ==, <, <=, >,
                >=, /, ^, \, ∪, transpose, convert, Array, Vector, Matrix

import Base.Broadcast: broadcasted, DefaultArrayStyle

import LinearAlgebra: ldiv!, mul!, rank, cond, qr



include("LinearAlgebra/LinearAlgebra.jl")
include("Operators/Operators.jl")
include("greensfun.jl")
include("fractals.jl")
include("clustertree.jl")

# if isdir(Pkg.dir("TikzGraphs"))
#     include("introspect.jl")
# end

include("plot.jl")
include("show.jl")



end #module
