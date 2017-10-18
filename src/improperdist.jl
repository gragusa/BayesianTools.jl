using Distributions
import Distributions: ContinuousUnivariateDistribution, pdf, logpdf

struct Improper{T} <: ContinuousUnivariateDistribution
    l::T
    u::T
    function Improper{T}(l::Number, u::Number) where T<:Number
        @assert any((!isfinite(l), !isfinite(u))) "At least one of the bounds must be Inf"
        @assert l < u "Lower bound must be smaller than upper bound"
        new{Float64}(float(l), float(u))
    end
end

Improper() = Improper{Float64}(-Inf, +Inf)
Improper(a, b) = Improper{Float64}(a, b)

function insupport(d::Improper, x::Real)
    (x >= d.l && x <= d.u) ? true : false
end

Distributions.pdf(d::Improper, x::Real) = insupport(d, x) ? 1.0 : 0.0
Distributions.logpdf(d::Improper, x::Real) = insupport(d, x) ? 0.0 : Inf
Base.minimum(d::Improper) = d.l
Base.maximum(d::Improper) = d.u
export Improper
