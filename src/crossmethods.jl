using .Links
using .ProductDistributions
BayesianTools.Links.link(p::ProductDistribution, x::Vector) = map(BayesianTools.Links.link, p, x)
BayesianTools.Links.invlink(p::ProductDistribution, x::Vector) = map(BayesianTools.Links.invlink, p, x)
BayesianTools.Links.logjacobian(p::ProductDistribution, x::Vector) = map(BayesianTools.Links.logjacobian, p, x)
