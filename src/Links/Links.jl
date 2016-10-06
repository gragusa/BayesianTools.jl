module Links
#################### TransformDistribution ####################

import Distributions: Continuous, ContinuousDistribution,
ContinuousUnivariateDistribution, ContinuousMultivariateDistribution,
## Univariate Distributions
Truncated, Cauchy, Gumbel, Laplace, Logistic, NoncentralT, Normal,
NormalCanon, TDist, Beta, BetaPrime, Chi, Chisq, Erlang, Exponential,
FDist, Frechet, Gamma, InverseGamma, InverseGaussian, Kolmogorov,
LogNormal, NoncentralChisq, NoncentralF, Rayleigh, Weibull, KSOneSided,
NoncentralBeta

import StatsFuns: logit

using Reexport
@reexport using Distributions

typealias TransformDistribution{T<:ContinuousUnivariateDistribution}
  Union{T, Truncated{T}}

function link(d::TransformDistribution, x::Real)
  a, b = minimum(d), maximum(d)
  lowerbounded, upperbounded = isfinite(a), isfinite(b)
  if lowerbounded && upperbounded
    logit((x - a) / (b - a))
  elseif lowerbounded
    log(x - a)
  elseif upperbounded
    log(b - x)
  else
    x
  end
end

function invlink(d::TransformDistribution, x::Real)
  a, b = minimum(d), maximum(d)
  lowerbounded, upperbounded = isfinite(a), isfinite(b)
  if lowerbounded && upperbounded
    (b - a) * invlogit(x) + a
  elseif lowerbounded
    exp(x) + a
  elseif upperbounded
    b - exp(x)
  else
    x
  end
end

function logjacobian(d::TransformDistribution, x::Real)
    # x - a >= 0
    # b - x >= 0
    a, b = minimum(d), maximum(d)
    lowerbounded, upperbounded = isfinite(a), isfinite(b)
    if lowerbounded && upperbounded
      - log(b - x)
    elseif lowerbounded
      - log(x - a)
    elseif upperbounded
      log(abs((b - a)/((a - x)*(x - b))))
    else
      0.0
    end
end


#################### RealDistribution ####################

typealias RealDistribution
          Union{Cauchy, Gumbel, Laplace, Logistic, NoncentralT, Normal,
                NormalCanon, TDist}

link(d::RealDistribution, x::Real) = x
invlink(d::RealDistribution, x::Real) = x
logjacobian(d::RealDistribution, x::Real) = zero(eltype(x))
# logpdf(d::RealDistribution, x::Real, transform::Bool) = logpdf(d, x)


#################### PositiveDistribution ####################

typealias PositiveDistribution
          Union{BetaPrime, Chi, Chisq, Erlang, Exponential, FDist, Frechet,
                Gamma, InverseGamma, InverseGaussian, Kolmogorov, LogNormal,
                NoncentralChisq, NoncentralF, Rayleigh, Weibull}

link(d::PositiveDistribution, x::Real) = log(x)

invlink(d::PositiveDistribution, x::Real) = exp(x)

logjacobian(d::PositiveDistribution, x::Real) = log(x)

# function  logpdf(d::PositiveDistribution, x::Real, transform::Bool)
#   lp = logpdf(d, x)
#   transform ? lp + log(x) : lp
# end


#################### UnitDistribution ####################

typealias UnitDistribution
          Union{Beta, KSOneSided, NoncentralBeta}

link(d::UnitDistribution, x::Real) = logit(x)

invlink(d::UnitDistribution, x::Real) = invlogit(x)

## logjacobian is d link / d x = log (1/(1-x))
logjacobian(d::UnitDistribution, x::Real) = -log(x)-log(1.0-x)



invlogit(x::Real) = 1.0 / (exp(-x) + 1.0)
invlogit(x::AbstractArray) = map(invlogit, x)




export link, invlink, logjacobian

end
