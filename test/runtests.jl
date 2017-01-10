using Base.Test
using BayesianTools
using StatsFuns

const T = typeof(1.0)

using BayesianTools.ProductDistributions

d = ProductDistribution(Normal(0,1), Uniform(0,1), InverseGamma(0.01, 0.01))

@test length(d) == 3
@test insupport(d, [0., 1., 1.])
@test insupport(d, [0., .1, 21.])
@test !insupport(d, [0., 2., .1])
@test !insupport(d, [0., -2., .1])
@test !insupport(d, [0., .1, -2])

@test logpdf(d, [0., 0.5, 0.5]) .== logpdf(d.marginals[1], 0.) + logpdf(d.marginals[2], 0.5) + logpdf(d.marginals[3], 0.5)
@test pdf(d, [0., 0.5, 0.5]) ≈ pdf(d.marginals[1], 0.)*pdf(d.marginals[2], 0.5)*pdf(d.marginals[3], 0.5)


using BayesianTools.Links

for d in [Uniform(0,1), Beta(2,2)]
    @test link(d, 0.1) == StatsFuns.logit(0.1)
    @test invlink(d, link(d, 0.1)) ≈ 0.1
end

for d in [Normal(0,1), TDist(2)]
    @test link(d, 0.1) == 0.1
    @test invlink(d, link(d, 0.1)) ≈ 0.1
end
d = ProductDistribution(Normal(0,1), Uniform(0,1), InverseGamma(0.01, 0.01))
@test link(d, invlink(d, [0.0, -1., 2.])) == [0.0, -1., 2.]
