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
N_control = convert(Array{Float64,1},TE_moms.N_control)

Q_moms = CSV.read("../Data/QuarterlyMoms.csv")
E_mom = convert(Array{Float64,1},Q_moms.LFP)/100
A_mom = convert(Array{Float64,1},Q_moms.Participation)/100
A2_mom = convert(Array{Float64,1},Q_moms.Receipt)
Y_mom = convert(Array{Float64,1},Q_moms.TotInc)

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
wghts = [ones(N1)/abs(mean(E_mom)); ones(N1)/abs(mean(A_mom)); zeros(N1)/abs(mean(A2_mom)); 5*ones(4)/abs(mean(XGmom)); N_control/(mean(N_control)*abs(mean(TE_moms0)))]
wghts_labor = copy(wghts);
wghts_labor[(3*N1+1):end] .= 0
wghts_TE = copy(wghts);
wghts_TE[1:3*N1] .= 0

# set up the parameter object
# :αc, :αθ, :αH, :αA, :β, :δI, :δθ, :ϵ, :τ, :pc, :wq, :αWR)
αc = 3. *ones(4)
αθ = 0.
αH = 1. *ones(4)
αA = 1. *ones(4)
β = 0.98
σA = 1.
δI = [-2.,0.]
δθ = 0.9
ϵ = 1.5
τ = [0.5,0.7]
#pc = [1.5,0.]
pc = [3.]
wq = 3. *12
αWR = 0.5
np = (αc = 4, αθ = 1, αH = 4, αA = 4, β = 1, σA = 1, δI = 2, δθ = 1, ϵ = 1, τ = 2, pc = 1, wq = 1, αWR = 1)
lb = (αc = zeros(4), αθ = 0, αH = -Inf*ones(4),αA = -Inf*ones(4),β = 0.4, σA = 0.01, δI = [-5,-0.2],δθ = 0, ϵ = 0,τ = zeros(2),pc = [0.],wq = 0.1, αWR = 0)
ub = (αc = Inf*ones(4), αθ = Inf, αH = Inf*ones(4),αA = Inf*ones(4),β = 1, σA = 10, δI = [3,0.1],δθ = 1.5, ϵ = Inf,τ = 0.99*ones(2),pc = [10.],wq = Inf,αWR = Inf)

pars = Parameters(np,lb,ub,αc,αθ,αH,αA,β,σA,δI,δθ,ϵ,τ,pc,wq,αWR)

labor_block = [:αc,:αH,:αA,:β,:σA,:wq,:αWR]
labor_block2 = [:αc,:αH,:αA,:β,:αWR] #<- ok, I wonder if there's a better way to do this,
#TE_block = [:δI,:δθ,:ϵ,:τ,:pc]
TE_block = [:δI,:δθ,:ϵ,:pc]
vlist = [:αc,:αθ,:αH,:αA,:β,:σA,:δI,:δθ,:ϵ,:τ,:pc,:wq,:αWR]

println("Just getting Labor Parameters")
opt,x0 = GetOptimization(pars,Mod1,labor_block,moms0,wghts_labor,5,lengths,TE_index,maxevals = 200,LBFGS = 1)
res = optimize(opt,x0)

#opt,x0 = GetOptimization(pars,Mod1,labor_block,moms0,wghts_labor,5,lengths,TE_index,maxevals = 500,LBFGS = 0)
#res = optimize(opt,x0)

break

function CriterionSpecial(x)
    Mod1.αc[1] = x[1]
    Mod1.αA[1] = x[2]
    Mod1.αH[1] = x[3]
    Mod1.β = x[4]
    Mod1.σA = x[5]
    Random.seed!(1212)
    SolveModel!(Mod1)
    E,A,A2,XG,skill = MomentsBaseline(Mod1,5,lengths,TE_index)
	mom_sim = [E;A;A2;XG;skill]
	Q = sum(wghts_labor.* (mom_sim .- moms0).^2)
    println(Q)
    return Q
end

opt_special = Opt(:LN_NELDERMEAD,5)
lower_bounds!(opt_special,[0,0,0,0.4,0.01])
upper_bounds!(opt_special,[10.,10.,10.,1.,10.])
min_objective!(opt_special,(x,g)->CriterionSpecial(x))
x0 = [Mod1.αc[1],Mod1.αA[1],Mod1.αH[1],Mod1.β,Mod1.σA]
res = optimize(opt_special,x0)

E,A,A2,XG,th,Y = MomentsBaseline(Mod1,5,lengths,TE_index);
D0 = Q_moms[:,[:Year,:Quarter,:Site,:Treatment,:LFP,:Participation]]
D0.Case = "Data"
D1 = Q_moms[:,[:Year,:Quarter,:Site,:Treatment]]
D1.LFP = E*100
D1.Participation = A*100
D1.Case = "Model"
CSV.write("ModelFit.csv",[D0;D1])

break
println("Just getting Production Parameters")
opt,x0 = GetOptimization(pars,Mod1,TE_block,moms0,wghts_TE,5,lengths,TE_index,maxevals = 200,LBFGS = 1,solve=false)
res = optimize(opt,x0)

opt,x0 = GetOptimization(pars,Mod1,TE_block,moms0,wghts_TE,5,lengths,TE_index,maxevals = 500,LBFGS = 0,solve=false)
res = optimize(opt,x0)

println("Starting Full")
println("FYI - subsidy parameters are:")
display(Mod1.τ)
pars.αθ = 0.1
opt,x0 = GetOptimization(pars,Mod1,labor_block,moms0,wghts_labor,5,lengths,TE_index,maxevals=200,LBFGS=1)
res = optimize(opt,x0)

opt,x0 = GetOptimization(pars,Mod1,vlist,moms0,wghts,5,lengths,TE_index,maxevals=200,LBFGS=1)
res = optimize(opt,x0)

opt,x0 = GetOptimization(pars,Mod1,vlist,moms0,wghts,5,lengths,TE_index,tightness = 1e-5,maxevals=10000)
res = optimize(opt,x0)
