using Test
using BayesianTools
using StatsFuns
using Random
const T = typeof(1.0)

using BayesianTools.ProductDistributions

d1 = ProductDistribution(Normal(0,1), Uniform(0,1), InverseGamma(0.01, 0.01))

@test length(d1) == 3
@test insupport(d1, [0., 1., 1.])
@test insupport(d1, [0., .1, 21.])
@test !insupport(d1, [0., 2., .1])
@test !insupport(d1, [0., -2., .1])
@test !insupport(d1, [0., .1, -2])

@test logpdf(d1, [0., 0.5, 0.5]) .== logpdf(d1.marginals[1], 0.) + logpdf(d1.marginals[2], 0.5) + logpdf(d1.marginals[3], 0.5)
@test pdf(d1, [0., 0.5, 0.5]) ≈ pdf(d1.marginals[1], 0.)*pdf(d1.marginals[2], 0.5)*pdf(d1.marginals[3], 0.5)

d1 = Improper(0, +Inf)
d2 = Improper(-Inf, +Inf)
d3 = Improper(-Inf, 5)
d4 = Improper(5, +Inf)
d5 = Improper()

@test d2===d5
@test_throws AssertionError Improper(0,1)
@test_throws AssertionError Improper(+Inf,-Inf)

@test !insupport(d1, -1)
@test insupport(d1, 1)

@test minimum(d1) == 0
@test minimum(d2) == -Inf
@test minimum(d3) == -Inf
@test minimum(d4) == 5

@test maximum(d1) == +Inf
@test maximum(d2) == +Inf
@test maximum(d3) == 5.0
@test maximum(d4) == +Inf

@test pdf(d1, 100) == 1.0
@test logpdf(d1, 100) == 0.0
@test pdf(d1, -100) == 0.0
@test logpdf(d1, -100) == Inf

@test pdf(d2, 100) == 1.0
@test logpdf(d2, 100) == 0.0
@test pdf(d2, -100) == 1.0
@test logpdf(d2, -100) == 0.0

@test pdf(d4, 7) == 1.0
@test logpdf(d4, 7) == 0.0
@test pdf(d4, 4) == 0.0
@test logpdf(d4, 4) == Inf

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
@test link(d, invlink(d, [0.0, -1., 2.])) ≈ [0.0, -1., 2.] atol=1e-15

d = ProductDistribution(Normal(0,1), Uniform(0,1), InverseGamma(0.01, 0.01), Improper(0,+Inf))
@test link(d, invlink(d, [0.0, -1., 2., -2.])) ≈ [0.0, -1., 2., -2.] atol=1e-15

d = ProductDistribution(Normal(), Uniform())
Random.seed!(1)
x2 = hcat([rand( d ) for i in 1:1000]...)
@test maximum(x2[2,:])<=1.0
@test maximum(x2[1,:])>=0.0

@test(params(d) == (0.0, 1.0, 0.0, 1.0))
