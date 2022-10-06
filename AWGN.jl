using Random, Distributions, SpecialFunctions
Random.seed!(1048)
qfunc(x::Float64)::Float64 = 0.5*erfc(x/sqrt(2)) #0と1に分布のピークがある．
qfuncinv(y::Float64)::Float64 = erfcinv(2*y)*sqrt(2)
bertosigma(ber::Float64)::Float64 = 0.5/qfuncinv(ber)
#%%
d = Normal(0.0, bertosigma(0.01))
noise = rand(d, (10,20))
