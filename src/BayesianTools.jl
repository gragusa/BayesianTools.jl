module BayesianTools

using Reexport

include("Links/Links.jl")
include("ProductDistributions/ProductDistributions.jl")
include("crossmethods.jl")
#export ProductDistribution, link, invlink


end # module
