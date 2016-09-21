using ProductDistributions
using Base.Test

using ProductDistributions

const T = typeof(1.0)


d = ProductDistribution(Normal(0,1), Uniform(0,1))

n = 10000

x = Array(T, 2, n)

rand!(d, x)

μ = mean(x, 2)
σ²= var(x, 2)

@test isa(μ, Array{T, 2})
@test isa(σ², Array{T, 2})

@test !insupport(d, [0., 2.])
@test insupport(d, [0., 1.])

@test logpdf(d, [0., 0.5]) .== logpdf(d.marginals[1], 0.) + logpdf(d.marginals[2], 0.5)
