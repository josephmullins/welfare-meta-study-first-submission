include("BaselineModel.jl")
include("EstimationBaseline.jl")
include("SetupBaseline.jl")
using NLopt

# don't forget that this needs to be estimated
Γ = exp.(0 .+ -0.1*(1:18)) #<- we need to figure this out and decide what we're doing!!!
#Γ = 0.1*ones(18)
#pars = (αc=1.,gN = 0. *ones(2), gF = ones(2)*0.2,wq = 2.,σC = 1.,σH = 1.,αWR = 0.5,αF = 2., αA = 1.,αH=1.,Γ=Γ,β=0.)
#pars = (αc=1. ,gN = 0. *ones(2), gF = ones(2)*0.2,wq = 2.,σC = 1.,σH = 1.,αWR = 0.5,αF = 2., αA = 1. *ones(num_sites),αH=1. *ones(num_sites),Γ=Γ,β=0.)

pars = parameters()
vlist = [:αc, :αH, :αA, :σH, :σC, :αWR,:αWR2, :αF, :β,:wq]


# # first an experiment to see what it takes to fit one site well

# i=2
# vsite1 = [:αc,:αH,:αA,:αF,:β] #,:αWR,:αWR2] #what if we went with just one alphaH
# vsite2 = [:αc,:αH,:αA,:σH,:σC,:αF,:β,:wq] # :αWR,:αWR2]#,:β]
# vsite2 = [:αc,:αH,:αA,:σH,:σC,:αF,:wq,:αWR,:αWR2] # :αWR,:αWR2]#,:β]
#
# pars,moms1,moms0 = FitSite(vsite1,site_list,budget,moments,wghts,site_features,i)
# pars,moms1,moms0 = FitSite(pars,vsite2,site_list,budget,moments,wghts,site_features,i)
#
# colors = ["blue","green","red"]
# T = site_features.T[i]
# for a=1:site_features.n_arms[i]
#     figure("AFDC")
#     #subplot(2,4,i)
#     title(String(site_list[i]))
#     plot(moms0[1:T,a],color=colors[a])
#     plot(moms1[1:T,a],color=colors[a],linestyle="--")
#     figure("LFP")
#     #subplot(2,4,i)
#     title(String(site_list[i]))
#     plot(moms0[T+1:2*T,a],color=colors[a])
#     plot(moms1[T+1:2*T,a],color=colors[a],linestyle="--")
# end


vlist0 = [:αc, :αH, :αA, :αF]
opt,x0 = GetOptimization(pars,vlist0,site_list,budget,moments,wghts,site_features)
#np = length(x0)

#g = zeros(np)

#Criterion(x0,g,pars,vlist0,site_list,budget,moments,wghts,site_features)

res = optimize(opt,x0)
pars = UpdatePars(res[2],pars,vlist0)

opt,x0 = GetOptimization(pars,vlist,site_list,budget,moments,wghts,site_features)
res2 = optimize(opt,x0)
pars2 = UpdatePars(res2[2],pars,vlist)

# vlist2 = [vlist;:wq]
# opt,x0=GetOptimization(pars2,vlist2,site_list,budget,moments,wghts,site_features)
# res2 = optimize(opt,x0)
# pars2 = UpdatePars(res2[2],pars,vlist2)


moms=GetMomentsAll(pars2,site_list,budget,moments,wghts,site_features)
#InspectTreatFit(moms,moments,site_features,site_list)
InspectModelFit(moms,data_moments,site_features,site_list)

break
opt,x0 = GetOptimization(pars,[:αH],site_list,budget,moments,wghts,site_features)
res = optimize(opt,x0)
pars = UpdatePars(res[2],pars,[:αH])

vlist1 = [:αc, :gN, :gF, :αA, :σH, :σC, :αWR,:αWR2, :αF,:β]
opt,x0 = GetOptimization(pars,vlist1,site_list,budget,moments,wghts,site_features)
res = optimize(opt,x0)
pars = UpdatePars(res[2],pars,vlist1)

vlist2 = [:αWR,:αWR2]
opt,x0 = GetOptimization(pars,vlist2,site_list,budget,moments,wghts,site_features)
res = optimize(opt,x0)
pars = UpdatePars(res[2],pars,vlist2)

break
#     αc,gN,gF,wq,σC,σH,αWR,αF = x[1:8]
x0 = [[3.,0.,0.2,3.,1.,1.,1.,-2.]; 2*ones(8); 2*ones(8); zeros(8)]
np = length(x0)
lower_bounds = [[0.,-Inf,-Inf,0.,0.,0.,0.,-Inf];-Inf*ones(8);-Inf*ones(8);-Inf*ones(8)]
upper_bounds = Inf*ones(np)
g = zeros(np)



opt = Opt(:LD_LBFGS,np)
opt2 = Opt(:LN_NELDERMEAD,np)
lower_bounds!(opt,lower_bounds)
upper_bounds!(opt,upper_bounds)
min_objective!(opt,(x,g)->Criterion(x,g,site_list,budget,moments,wghts,site_features,Γ))
lower_bounds!(opt2,lower_bounds)
upper_bounds!(opt2,upper_bounds)
min_objective!(opt2,(x,g)->Criterion(x,site_list,budget,moments,wghts,site_features,Γ))
res = optimize(opt,x0)

pars=GetParameters(res[2],num_sites,Γ)
moms=GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)

for s in site_list
    m1 = getfield(moments,s)
    m2 = getfield(moms,s)
    scatter(m1,m2)
end

InspectModelFit(moms,moments,site_features,site_list)

# model prediction: shouldn't the difference in treatment ewffect of work requirements decrease when welfare use decreases??
# seems a little at odds here, so let's revist
