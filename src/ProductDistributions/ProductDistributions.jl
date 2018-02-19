module ProductDistributions

import Distributions: ContinuousDistribution, ContinuousMultivariateDistribution,
## Univariate Distributions
logpdf, pdf, rand!
using Reexport
@reexport using Distributions
@reexport using BayesianTools

import Base: first, start, next, done

immutable ProductDistribution{T} <: ContinuousMultivariateDistribution
    marginals::T
end

Base.length(d::ProductDistribution) = length(d.marginals)

function ProductDistribution(args...)
    ProductDistribution((args...,))
end

ProductDistribution(x::Distribution) = ProductDistribution([x])

function Distributions.insupport{T<:Real}(d::ProductDistribution, x::AbstractVector{T})
    insup = true
    for i in eachindex(d.marginals)
        insupport(d.marginals[i], x[i]) || (insup = false; break)
    end
    insup
end

function Distributions.logpdf{T<:Real}(d::ProductDistribution, x::AbstractVector{T})
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

function Distributions.pdf{T<:Real}(d::ProductDistribution, x::AbstractVector{T})
    exp(logpdf(d, x))
end

function Distributions.rand!{T<:Real}(d::ProductDistribution, x::DenseMatrix{T})
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
