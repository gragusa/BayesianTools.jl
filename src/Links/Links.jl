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

# function logpdf(d::TransformDistribution, x::Real, transform::Bool)
#   lp = logpdf(d, x)
#   if transform
#     a, b = minimum(d), maximum(d)
#     lowerbounded, upperbounded = isfinite(a), isfinite(b)
#     if lowerbounded && upperbounded
#       lp += log((x - a) * (b - x) / (b - a))
#     elseif lowerbounded
#       lp += log(x - a)
#     elseif upperbounded
#       lp += log(b - x)
#     end
#   end
#   lp
# end


#################### RealDistribution ####################

typealias RealDistribution
          Union{Cauchy, Gumbel, Laplace, Logistic, NoncentralT, Normal,
                NormalCanon, TDist}

link(d::RealDistribution, x::Real) = x

invlink(d::RealDistribution, x::Real) = x

# logpdf(d::RealDistribution, x::Real, transform::Bool) = logpdf(d, x)


#################### PositiveDistribution ####################

typealias PositiveDistribution
          Union{BetaPrime, Chi, Chisq, Erlang, Exponential, FDist, Frechet,
                Gamma, InverseGamma, InverseGaussian, Kolmogorov, LogNormal,
                NoncentralChisq, NoncentralF, Rayleigh, Weibull}

link(d::PositiveDistribution, x::Real) = log(x)

invlink(d::PositiveDistribution, x::Real) = exp(x)

# function  logpdf(d::PositiveDistribution, x::Real, transform::Bool)
#   lp = logpdf(d, x)
#   transform ? lp + log(x) : lp
# end


#################### UnitDistribution ####################

typealias UnitDistribution
          Union{Beta, KSOneSided, NoncentralBeta}

link(d::UnitDistribution, x::Real) = logit(x)

invlink(d::UnitDistribution, x::Real) = invlogit(x)

# function logpdf(d::UnitDistribution, x::Real, transform::Bool)
#   lp = logpdf(d, x)
#   transform ? lp + log(x * (1.0 - x)) : lp
# end

invlogit(x::Real) = 1.0 / (exp(-x) + 1.0)
invlogit(x::AbstractArray) = map(invlogit, x)

export link, invlink

end
