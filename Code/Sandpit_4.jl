using DataFrames
using PyPlot
include("BaselineModel.jl")
include("EstimationRoutines.jl")
using DelimitedFiles


Earnings = reshape(readdlm("earnings"),4,73)
budget1 = reshape(readdlm("budget"),4,3,3,18,18*4,2,2,2)
Budget_Ageout = reshape(readdlm("budget_ageout"),4,73,2,2)
Foodstamps_receipt=reshape(readdlm("foodstamps_receipt"),4,3,3,18,18*4,2,2,2)
#findmax(budget1[:,3,:,:,:,:,:,:])
TE_moms = CSV.read("../Data/ChildTreatmentEffectsScore.csv")[1:16,:]
TE_moms[:Site] = [ones(2);2*ones(2);3*ones(7);4*ones(5)]
TE_index = convert(Array{Int64,2},TE_moms[:,[:Site,:Treatment,:AgeMin,:AgeMax]])
TE_moms0 = convert(Array{Float64,1},TE_moms.FacScore)

Q_moms = CSV.read("../Data/QuarterlyMoms.csv")
E_mom = convert(Array{Float64,1},Q_moms.LFP)/100
A_mom = convert(Array{Float64,1},Q_moms.Participation)/100
A2_mom = convert(Array{Float64,1},Q_moms.Receipt)

τ = 0.1*ones(4,3)
XG0 = [1082.75 1351.5; 94.33 200.25; 857.67 1089.67] #<- average subsidy expense for each program
XGm = XG0/12 #<- monthly
Xcm = [57.3 71;21. 20] #<- monthly payment from CTJF and FTP members
τ[1:2,1:2] = XGm[1:2,:]./(XGm[1:2,:] .+ Xcm)
XGmom = [(XG0[:,2] .- XG0[:,1])./XG0[:,1]; XG0[3,1]/XG0[1,1]]

TimeLimit_Ind=[false true true; false true true; false false false; false false false]

TimeLimits=[0 7 7; 0 8 8; 0 0 0; 0 0 0]
lengths = [16,18,12,12]
Work_Reqs_Ind=[false true true; false true true; false true false;false true false];

# set up the model
Mod1=initialize_model();

# set up moments and weights
moms0 = [E_mom;A_mom;A2_mom;XGmom;TE_moms0]
N1 = length(E_mom)
N2 = length(TE_moms0)
# forget about targeting A2 for now.
wghts = [ones(N1)/sqrt(abs(mean(E_mom))); ones(N1)/sqrt(abs(mean(A_mom))); zeros(N1)/sqrt(abs(mean(A2_mom))); ones(4)/sqrt(abs(mean(XGmom))); ones(N2)/sqrt(abs(mean(TE_moms0)))]

# set up the parameter object
# :αc, :αθ, :αH, :αA, :β, :δI, :δθ, :ϵ, :τ, :pc, :wq, :αWR)
αc = 0.1
αθ = 0.1
αH = 0.1*ones(4)
αA = 0.1*ones(4)
β = 0.98
δI = [-10.,0.]
δθ = 0.9
ϵ = 1.5
τ = [0.3,0.3]
pc = [0.,0.]
wq = 3. *12
αWR = 0.1
np = (αc = 1, αθ = 1, αH = 4, αA = 4, β = 1, δI = 2, δθ = 1, ϵ = 1, τ = 2, pc = 2, wq = 1, αWR = 1)
lb = (αc = 0, αθ = 0, αH = -Inf*ones(4),αA = -Inf*ones(4),β = 0.4, δI = [-5,-0.2],δθ = 0, ϵ = 0,τ = zeros(2),pc = [-2,-0.2],wq = 0.1, αWR = 0)
ub = (αc = Inf, αθ = Inf, αH = Inf*ones(4),αA = Inf*ones(4),β = 1, δI = [1,0.2],δθ = 1.5, ϵ = Inf,τ = 0.99*ones(2),pc = [1,0.2],wq = Inf,αWR = Inf)

pars = Parameters(np,lb,ub,αc,αθ,αH,αA,β,δI,δθ,ϵ,τ,pc,wq,αWR)

labor_block = [:αc,:αH,:αA,:β,:wq,:αWR]
#labor_block2 = [:αH,:αA,:β,:αWR] #<- ok, I wonder if there's a better way to do this,
wghts_alt = copy(wghts)
wghts_alt[2*N1+1:end] .= 0

# ----- Find initial guess of preference parameters
opt1,x0 = GetOptimization(pars,Mod1,labor_block,moms0,wghts_alt,5,lengths,TE_index,maxevals=500)
println(" ------ Round 1 ------")
res1 = NLopt.optimize(opt1,x0)
println(" ------ Round 2 ------")
res2 = NLopt.optimize(opt1,res1[2])
println(" ------ Round 3 ------")
res3 = NLopt.optimize(opt1,res2[2])

# ---- Find initial guess of production parameters
pars.δI[1] = -2.
wghts_child = copy(wghts);
wghts_child[1:3*N1] .= 0
child_block = [:δI,:δθ,:ϵ,:τ,:pc,:wq] #<- use wq?
opt2,x0 = GetOptimization(pars,Mod1,child_block,moms0,wghts_child,5,lengths,TE_index,maxevals=500,no_solve=true) #<- just pick production pars to hit treatment effects, nothing else
res4 = NLopt.optimize(opt2,x0)
res5 = NLopt.optimize(opt2,res4[2])


# see if we can get FTP only
# wghts_ftp = wghts_alt
# wghts_ftp[1:32] .= 0
# wghts_ftp[69:end] .= 0
# wghts_ftp[(N1+33):(N1+69)] = wghts[(N1+33):(N1+69)]
# opt2,x0 = GetOptimization(pars,Mod1,labor_block,moms0,wghts_ftp,5,lengths,TE_index)
# res2 = optimize(opt2,x0)

vlist = [:αc,:αθ,:αH,:αA,:β,:δI,:δθ,:ϵ,:τ,:pc,:wq,:αWR]
opt,x0 = GetOptimization(pars,Mod1,vlist,moms0,wghts,5,lengths,TE_index)
res = optimize(opt,x0)
writedlm("new_est_JM.csv",res[2])


break

wghts_alt2 = copy(wghts)
wghts_alt1 = copy()


#E,A,A2,XG,skill_moms = MomentsBaseline(Mod1,5,lengths,TE_index)
