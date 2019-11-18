cd("/Users/FilipB/github/welfare-meta-study/Code")
using DataFrames
using PyPlot
using Revise
using Optim
includet("BaselineModel.jl")
includet("EstimationRoutines.jl")
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
wghts = [ones(N1)/sqrt(abs(mean(E_mom))); ones(N1)/sqrt(abs(mean(A_mom))); zeros(N1)/sqrt(abs(mean(A2_mom))); ones(4)/sqrt(abs(mean(XGmom))); ones(N2)/sqrt(abs(mean(TE_moms0)))]

# set up the parameter object
# :αc, :αθ, :αH, :αA, :β, :δI, :δθ, :ϵ, :τ, :pc, :wq, :αWR)

#Param1=[0.1 0.1 0.1 0.1  0.98 -10. 0.0 0.9 0.9 0.1 0.1 0.0 0.0 36.0 0.1]
Param1=CSV.read("new_est.csv", header=false)
Param1=Param1[:,1]

αc = Param1[1]# 0.1
αθ = Param1[2]# 0.1
αH = [Param1[3],Param1[4],Param1[5],Param1[6]]#Param1[3].*ones(4)# 0.1*ones(4)
αA = [Param1[7],Param1[8],Param1[9],Param1[10]]#Param1[4].*ones(4)# 0.1*ones(4)
β = Param1[11]# 0.98
δI = [Param1[12],Param1[13]]#[Param1[6],Param1[7]]# [-10.,0.]
δθ = Param1[14]# 0.9
ϵ = Param1[15]# 0.9
τ1 = [Param1[16],Param1[17]]# [0.1,0.1]
pc = [Param1[18],Param1[19]]# [0.,0.]
wq = Param1[20]# 3. *12
αWR = Param1[21]# 0.1

Bnds=10000

np = (αc = 1, αθ = 1, αH = 4, αA = 4, β = 1, δI = 2, δθ = 1, ϵ = 1, τ = 2, pc = 2, wq = 1, αWR = 1)
lb = (αc = 0, αθ = 0, αH = -Bnds*ones(4),αA = -Bnds*ones(4),β = 0, δI = [-15,-5],δθ = 0, ϵ = 0,τ = zeros(2),pc = -5*ones(2),wq = 0.1, αWR = 0)
ub = (αc = Bnds, αθ = Bnds, αH = Bnds*ones(4),αA = Bnds*ones(4),β = 1, δI = 5*ones(2),δθ = 1.5, ϵ = Bnds,τ = 0.99*ones(2),pc = 5*ones(2),wq = Bnds,αWR = Bnds)

vlist = [:αc,:αθ,:αH,:αA,:β,:δI,:δθ,:ϵ,:τ,:pc,:wq,:αWR]

x1=[αc, αθ,  αH[1], αH[2],  αH[3], αH[4],
    αA[1],αA[2], αA[3],αA[4], β, δI[1], δI[2],  δθ, ϵ, τ[1],τ[2], pc[1],pc[2], wq, αWR]

pars = Parameters(np,lb,ub,αc,αθ,αH,αA,β,δI,δθ,ϵ,τ1,pc,wq,αWR)

labor_block = [:αc,:αH,:αA,:β]
labor_block2 = [:αH,:αA,:β] #<- ok, I wonder if there's a better way to do this,
wghts_alt = copy(wghts)
wghts_alt[2*N1+1:end] .= 0
xj=2