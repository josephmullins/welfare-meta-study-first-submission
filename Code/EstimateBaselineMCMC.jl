include("BaselineModel.jl")
include("MCMCRoutines.jl")
include("SetupBaseline.jl")

Γ = readdlm("FirstStageSIPP/Gamma_est")[:]
gF = readdlm("FirstStageSIPP/gF")[:] #<- we need to figure this out and decide what we're doing!!!

xm0,mpars = ModPars(Γ,gF)
xh0,hpars = HyperPars()
mlist = [:αc, :αH, :αA, :σH, :σC, :αWR,:αWR2, :αF, :β,:wq,:gN]
hlist = keys(hpars.pos)

ll = LogLikeMod(xm0,mpars,budget,moments,site_features)

Pm,Ph,Amhist,Ahhist,Lhist = GetChainMod(500,xm0,xh0,mlist,hlist,mpars,hpars,budget,moments,site_features)
