using Base.Test
using BayesianTools
using StatsFuns

const T = typeof(1.0)

using BayesianTools.ProductDistributions

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

@test logpdf(d, [0., 0.5]) .== logpdf(d.marginals[1].m, 0.) + logpdf(d.marginals[2].m, 0.5)

@test pdf(d, [0., 0.5]) ≈ pdf(d.marginals[1].m, 0.)*pdf(d.marginals[2].m, 0.5)


using BayesianTools.Links

for d in [Uniform(0,1), Beta(2,2)]
    @test link(d, 0.1) == StatsFuns.logit(0.1)
    @test invlink(d, link(d, 0.1)) ≈ 0.1
end

for d in [Normal(0,1), TDist(2)]
    @test link(d, 0.1) == 0.1
    @test invlink(d, link(d, 0.1)) ≈ 0.1
end
