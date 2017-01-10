using Distributions
import Distributions: ContinuousUnivariateDistribution, pdf, logpdf

immutable Improper{T} <: ContinuousUnivariateDistribution
    l::T
    u::T
    function Improper(l::Number, u::Number)
        @assert any((!isfinite(l), !isfinite(u))) "At least one of the boud must be Inf"
        @assert l < u "Lower bound must be smaller than upper bound"
        new{Float64}(float(l), float(u))
    end
end

Improper() = Improper{Float64}(-Inf, +Inf)
Improper(a, b) = Improper{Float64}(a, b)

function insupport(d::Improper, x::Real)
    (x >= d.l && x <= d.u) ? true : false
end

Distributions.pdf(d::Improper, x::Real) = insupport(d, x) ? 0.0 : 1.0
Distributions.logpdf(d::Improper, x::Real) = insupport(d, x) ? +Inf : 0.0
Base.minimum(d::Improper) = d.l
Base.maximum(d::Improper) = d.u
export Improper