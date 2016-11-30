module ProductDistributions

import Distributions: ContinuousDistribution, ContinuousMultivariateDistribution,
## Univariate Distributions
logpdf, pdf, rand!
using Reexport
@reexport using Distributions

import Base: first, start, next, done

immutable ProductDistribution <: ContinuousMultivariateDistribution
    marginals
end

ProductDistribution(z::AbstractArray) = ProductDistribution(z...)

function ProductDistribution(args...)
    @assert all(map(x->isa(x, Distribution{Univariate,Continuous}), args))
    ProductDistribution(map(x -> Marginal(x), args))
end

immutable Marginal{T}
    m::T
end

Base.length(d::ProductDistribution) = length(d.marginals)






# function ProductDistribution(z::AbstractArray)
#     ProductDistribution(map(x -> Marginal(x), z))
# end


function Distributions.insupport{T<:Real}(d::ProductDistribution, x::AbstractVector{T})
    all(map((i,x) -> insupport(d.marginals[i].m,x), 1:length(x), x))
end

function Distributions.logpdf{T<:Real}(d::ProductDistribution, x::AbstractVector{T})
    if Distributions.insupport(d, x)
        l = zero(T)
        @inbounds for i = 1:length(d.marginals)
            l += logpdf(d.marginals[i].m, x[i])
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
    @inbounds for i = 1:N
        for p = 1:P
            x[p, i] = rand(d.marginals[p].m)
        end
    end
    x
end

Base.start(d::ProductDistribution) = 1
Base.first(d::ProductDistribution) = d.marginals[1].m
Base.done(d::ProductDistribution, state) = state == length(d) + 1
Base.next(d::ProductDistribution, state) = (d.marginals[state].m, state + 1)


export ProductDistribution

end
