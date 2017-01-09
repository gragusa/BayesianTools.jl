# BayesianTools.jl
[![Build Status](https://travis-ci.org/gragusa/ProductDistributions.jl.svg?branch=master)](https://travis-ci.org/gragusa/ProductDistributions.jl)
[![Coverage Status](https://coveralls.io/repos/gragusa/ProductDistributions.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/gragusa/ProductDistributions.jl?branch=master)
[![codecov.io](http://codecov.io/github/gragusa/ProductDistributions.jl/coverage.svg?branch=master)](http://codecov.io/github/gragusa/ProductDistributions.jl?branch=master)

`BayesianTools.jl` is a Julia package with methods useful for Monte Carlo Markov Chain simulations. The package has two submodules: 

- `ProductDistributions`: defines a `ProductDistribution` type and related methods useful for defining and evaluating independent priors
- `Link`: usuful rescale MC proposals to the parameter space of the underlying prior

## Installation

The package is not registered, so it must be cloned:
```julia
Pkg.clone("https://github.com/gragusa/BayesianTools.jl.git")
```

## Usage

### ProductDistributions

The following code defines the distribution (up to a scaling constant) that results from multiplying a normal and a Beta
```julia
using BayesianTools.ProductDistributions
p = ProductDistribution(Normal(0,1), Beta(1.,1.))
```
To check whether an `Array{Float64}` is in the support of `p`
```julia
insupport(p, [.1,2.]) ## false
insupport(p, [.1,1.]) ## true
```
The likelihood and loglikelihood are
```julia
loglikelihood(p, [.1,.5]) # = loglikelihood(Normal(0,1), .1) + loglikelihood(Beta(1.,1.), .5)
likelihood(p, [.1,.5]) # = likelihood(Normal(0,1), .1) * likelihood(Beta(1.,1.), .5)
```




