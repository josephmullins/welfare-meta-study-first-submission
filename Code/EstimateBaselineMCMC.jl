include("BaselineModel.jl")
include("MCMCRoutines.jl")
include("SetupBaseline.jl")
using NLopt
Γ = readdlm("FirstStageSIPP/Gamma_est")[:]
gF = readdlm("FirstStageSIPP/gF")[:] #<- we need to figure this out and decide what we're doing!!!

xm0,mpars = ModPars(Γ,gF)
xh0,hpars = HyperPars()
mlist = [:αc, :αH, :αA, :σH, :σC, :αWR,:αWR2, :αF, :wq,:gN,:β]
hlist = keys(hpars.pos)

ll = LogLikeMod(xm0,mpars,budget,moments,site_features)

#Pm,Ph,Amhist,Ahhist,Lhist = GetChainMod(500,xm0,xh0,mlist,hlist,mpars,hpars,budget,moments,site_features)

# first get maximum likelihood estimates as a way to initiate the markov chain
opt = Opt(:LD_LBFGS,length(xm0))
xm0[mpars.pos.wq] = 3.
lower_bounds!(opt,[0.;-Inf*ones(16);0.;0.;-Inf*ones(24);0.1;-Inf;0])
upper_bounds!(opt,[Inf;Inf*ones(16);Inf;Inf;Inf*ones(24);Inf;Inf;1])
max_objective!(opt,(x,g)->LogLikeMod(x,g,mpars,budget,moments,site_features))
mres = optimize(opt,xm0) #<- this gives us a starting point to run the chain

# next get estimates of the hyperparameters
#(:αH, :σαH, :αA, :σαA, :αF, :σαF, :αWR, :σαWR, :αWR2, :σαWR2)
opt2 = Opt(:LD_LBFGS,length(xh0))
lower_bounds!(opt2,[-Inf,0.,-Inf,0.,-Inf,0.,-Inf,0.,-Inf,0])
upper_bounds!(opt2,Inf*ones(10))
max_objective!(opt2,(x,g)->LogLikeHyper(x,g,mres[2],mpars,hpars))
hres = optimize(opt2,xh0) #<- this gives us a starting point to run the chain

#mlist0 = [:αc,:αA,:αH,:αWR,:αWR2,:wq]
Pm,Ph,Amhist,Ahhist,Lhist = GetChainMod(5000,mres[2],hres[2],mlist,hlist,mpars,hpars,budget,moments,site_features)

Pm2,Ph2,Amhist2,Ahhist2,Lhist2 = GetChainMod(5000,mres[2],hres[2],mlist,hlist,mpars,hpars,budget,moments,site_features)

Pm3,Ph3,Amhist3,Ahhist3,Lhist3 = GetChainMod(50000,mres[2],hres[2],mlist,hlist,mpars,hpars,budget,moments,site_features)

#Pm2,Ph2,Amhist2,Ahhist2,Lhist2 = GetChainMod(10000,Pm[:,end],Ph[:,end],mlist,hlist,mpars,hpars,budget,moments,site_features)
