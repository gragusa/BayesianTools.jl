module ProductDistributions

import Distributions: ContinuousDistribution, ContinuousMultivariateDistribution,
## Univariate Distributions
logpdf, pdf, rand!
using Reexport
@reexport using Distributions
@reexport using BayesianTools

import Base: first, start, next, done

struct ProductDistribution{T} <: ContinuousMultivariateDistribution
    marginals::T
end

Base.length(d::ProductDistribution) = length(d.marginals)

function ProductDistribution(args...)
    ProductDistribution((args...,))
end

ProductDistribution(x::Distribution) = ProductDistribution([x])

function Distributions.insupport(d::ProductDistribution, x::AbstractVector{T}) where T<:Real
    insup = true
    for i in eachindex(d.marginals)
        insupport(d.marginals[i], x[i]) || (insup = false; break)
    end
    insup
end

function Distributions.logpdf(d::ProductDistribution, x::AbstractVector{T}) where T<:Real
    if Distributions.insupport(d, x)
        l = zero(T)
        @inbounds for i in 1:length(d.marginals)
            l += logpdf(d.marginals[i], x[i])
        end
    else
        l = convert(T, -Inf)
    end
    l
end

function Distributions.pdf(d::ProductDistribution, x::AbstractVector{T}) where T<:Real
    exp(logpdf(d, x))
end

function Distributions.rand!(d::ProductDistribution, x::DenseMatrix{T}) where T<:Real
    P, N = size(x)
    @assert P == length(d)
    @inbounds for i = 1:N
        for p = 1:P
            x[p, i] = rand(d.marginals[p])
        end
    end
    x
end

function Distributions.rand!{T<:Real}(d::ProductDistribution, x::AbstractVector{T})
    P = length(x)
    @assert P == length(d)
    for p = 1:P
        x[p] = rand(d.marginals[p])
    end
    x
end

Base.start(d::ProductDistribution) = 1
Base.first(d::ProductDistribution) = d.marginals[1]
Base.done(d::ProductDistribution, state) = state == length(d) + 1
Base.next(d::ProductDistribution, state) = (d.marginals[state], state + 1)



export ProductDistribution

end
