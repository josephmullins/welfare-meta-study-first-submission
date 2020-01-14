include("BaselineModel.jl")
include("EstimationBaseline.jl")
include("SetupBaseline.jl")
using NLopt


Γ = exp.(0 .+ -0.1*(1:18)) #<- we need to figure this out and decide what we're doing!!!
#Γ = 0.1*ones(18)
#pars = (αc=1.,gN = 0. *ones(2), gF = ones(2)*0.2,wq = 2.,σC = 1.,σH = 1.,αWR = 0.5,αF = 2., αA = 1.,αH=1.,Γ=Γ,β=0.)
pars = (αc=1.,gN = 0. *ones(2), gF = ones(2)*0.2,wq = 2.,σC = 1.,σH = 1.,αWR = 0.5,αF = 2., αA = 1. *ones(num_sites),αH=1. *ones(num_sites),Γ=Γ,β=0.)

#GetStaticProbs(pars,budget.CTJF[1,:,:,:,:],30.,1,5,3)
#moms = GetMoments(pars,budget.CTJF[1,:,:,:,:],site_features.prices[1,1],0,4,site_features.π0[1,:,:],3)
CriterionP(pars,site_list,budget,moments,wghts,site_features)
#αc,gN,gF,wq,σC,σH,αWR,αF = x[1:8]
#αH = x[9:(8+num_sites)]
#αA = x[(9+num_sites):(8+2*num_sites)]
x0 = [[1.,0.,0.2,2.,1.,1.,0.5,2.]; ones(8); ones(8)]
np = length(x0)
lower_bounds = [[0.,-Inf,-Inf,0.,0.,0.,0.,-Inf];-Inf*ones(8);-Inf*ones(8)]
upper_bounds = Inf*ones(np)
g = zeros(np)

df = ForwardDiff.gradient(x->test3(x),[1.,2.])


Criterion(x0,g,site_list,budget,moments,wghts,site_features,Γ)


opt = Opt(:LD_LBFGS,np)
lower_bounds!(opt,lower_bounds)
upper_bounds!(opt,upper_bounds)
min_objective!(opt,(x,g)->Criterion(x,g,site_list,budget,moments,wghts,site_features,Γ))
res = optimize(opt,x0)

pars=GetParameters(res[2],num_sites,Γ)
moms=GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)

for s in site_list
    m1 = getfield(moments,s)
    m2 = getfield(moms,s)
    scatter(m1,m2)
end
