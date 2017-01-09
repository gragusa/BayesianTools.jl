# BayesianTools.jl
[![Build Status](https://travis-ci.org/gragusa/BayesianTools.jl.svg?branch=master)](https://travis-ci.org/gragusa/BayesianTools.jl)
[![Coverage Status](https://coveralls.io/repos/gragusa/BayesianTools.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/gragusa/BayesianTools.jl?branch=master)
[![codecov.io](http://codecov.io/github/gragusa/BayesianTools.jl/coverage.svg?branch=master)](http://codecov.io/github/gragusa/BayesianTools.jl?branch=master)ssss

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

It is also possible to draw a sample from `p`
```julia
rand!(p, Array{Float64}(2,100))
```

### Example
This is a common use case for the `ProductDistributions` part of `BayesianTools`, MCMC simulations.

![](https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20y%20%26%20%5Csim%20N%28%5Cbeta_0%20&plus;%20%5Cbeta_1%20x%2C%20%5Csigma%5E2%29%20%5C%5C%20%5Cbeta_0%2C%20%5Cbeta_1%20%26%5Csim%20N%280%2C%201000%29%20%5C%5C%20%5Csigma%5E2%20%26%20%5Csim%20%5Cmathrm%7BinvGamma%7D%280.001%2C%200.001%29%20%5Cend%7Balign*%7D)

```julia
using BayesianTools.ProductDistributions

p = ProductDistribution(Normal(0,10), Normal(0,10), Inverse

f(theta) = loglikelihood(y, x, theta) + loglikelihood(p, theta)

## Load RDataset
using RDataset

## Estimation ssss


```


