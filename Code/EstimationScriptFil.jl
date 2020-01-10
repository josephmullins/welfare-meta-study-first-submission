using DataFrames
using PyPlot
using Revise
using ForwardDiff
using Optim
include("BaselineModel.jl")
includet("EstimationRoutines.jl")
includet("Budget_Function_Code.jl")
using DelimitedFiles
cd("/Users/FilipB/github/welfare-meta-study/Code")

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
moms0 = [E_mom;A_mom;A2_mom;XGmom;TE_moms0;Y_mom]
N1 = length(E_mom)
N2 = length(TE_moms0)
# forget about targeting A2 for now.
wghts = [ones(N1)/abs(mean(E_mom));
        ones(N1)/abs(mean(A_mom));
        zeros(N1)/abs(mean(A2_mom));
        5*ones(4)/abs(mean(XGmom));
        N_control/(mean(N_control)*abs(mean(TE_moms0)));
        zeros(N1)/abs(mean(Y_mom))]

wghts_Full = [ones(N1)/abs(mean(E_mom));
        ones(N1)/abs(mean(A_mom));
        ones(N1)/abs(mean(A2_mom));
        5*ones(4)/abs(mean(XGmom));
        N_control/(mean(N_control)*abs(mean(TE_moms0)));
        ones(N1)/abs(mean(Y_mom))]


ENames= fill("E",N1)
ANames= fill("A",N1)
A2Names= fill("A2",N1)
XGNames= fill("XG",4)
TENames= fill("TE",16)
YNames= fill("Y",N1)

Names_all=[ENames; ANames; A2Names; XGNames; TENames; YNames]

wghts_labor = copy(wghts);
wghts_labor[(3*N1+1):end] .= 0
wghts_TE = copy(wghts);
wghts_TE[1:3*N1] .= 0


Par1=CSV.read("Latest_Params.csv")
Par1=Par1[:,2]

# set up the parameter object
# :αc, :αθ, :αH, :αA, :β, :δI, :δθ, :ϵ, :τ, :pc, :wq, :αWR)
αc =  3. *ones(4)
αθ =  0.
αH =  1. *ones(4)
αA =  1. *ones(4)
β =  0.98
σA = 1.
δI =  [-2.,0.]
δθ =  0.9
ϵ =  1.5
τ =  [0.5,0.7]
#pc = [1.5,0.]
pc =  [3.]
wq =  3. *12
αWR =  0.5
np = (αc = 4, αθ = 1, αH = 4, αA = 4, β = 1, σA = 1, δI = 2, δθ = 1, ϵ = 1, τ = 2, pc = 1, wq = 1, αWR = 1)
lb = (αc = zeros(4), αθ = 0, αH = -Inf*ones(4),αA = -Inf*ones(4),β = 0.4, σA = 0.01, δI = [-5,-0.2],δθ = 0, ϵ = 0,τ = zeros(2),pc = [0.],wq = 0.1, αWR = 0)
ub = (αc = Inf*ones(4), αθ = Inf, αH = Inf*ones(4),αA = Inf*ones(4),β = 1, σA = 10, δI = [3,0.1],δθ = 1.5, ϵ = Inf,τ = 0.99*ones(2),pc = [10.],wq = Inf,αWR = Inf)

pars = Parameters(np,lb,ub,αc,αθ,αH,αA,β,σA,δI,δθ,ϵ,τ,pc,wq,αWR)



αc[4]

function pars_to_vec(pars)

    x1=[pars.αc[1], pars.αc[2],  pars.αc[3], pars.αc[4],
    pars.αθ,
    pars.αH[1], pars.αH[2],  pars.αH[3], pars.αH[4],
    pars.αA[1],pars.αA[2], pars.αA[3],pars.αA[4],
    pars.β,
    pars.σA,
    pars.δI[1], pars.δI[2],
    pars.δθ, pars.ϵ,
    pars.τ[1],pars.τ[2], pars.pc[1], pars.wq, pars.αWR]

    return x1
end

function pars_to_vec_LBlock(pars,Par1) # feed in best guesses only in labor block

    x1=[Par1[1], Par1[2],  Par1[3], Par1[4],
    pars.αθ,
    Par1[6], Par1[7], Par1[8], Par1[9],
    Par1[10],Par1[11], Par1[12],Par1[13],
    Par1[14],
    Par1[15],
    pars.δI[1], pars.δI[2],
    pars.δθ, pars.ϵ,
    pars.τ[1],pars.τ[2], pars.pc[1], Par1[23], Par1[24]]

    return x1
end

x0=pars_to_vec(pars)

x1=pars_to_vec_LBlock(pars,Par1)



Names_Pars=["alpha_c","alpha_c","alpha_c","alpha_c",
            "alpha_theta",
            "alpha_H","alpha_H","alpha_H","alpha_H",
            "alpha_A","alpha_A","alpha_A","alpha_A",
            "beta",
            "sigma_A",
            "delta_I","delta_I",
            "delta_theta","epsilon",
            "tau","tau","pc","wq","alpha_WR"]



labor_block = [:αc,:αH,:αA,:β,:σA,:wq,:αWR]
labor_block2 = [:αc,:αH,:αA,:β,:αWR] #<- ok, I wonder if there's a better way to do this,
#TE_block = [:δI,:δθ,:ϵ,:τ,:pc]
TE_block = [:δI,:δθ,:ϵ,:pc]
vlist = [:αc,:αθ,:αH,:αA,:β,:σA,:δI,:δθ,:ϵ,:τ,:pc,:wq,:αWR]

Criterion(x0,pars,Mod1,vlist,moms0,wghts_labor,5,lengths,TE_index, true, true)
@time Criterion(x1,pars,Mod1,vlist,moms0,wghts_labor,5,lengths,TE_index, true, true)

function pars_to_criterion(x2)
    return Criterion(x2,pars,Mod1,vlist,moms0,wghts_labor,5,lengths,TE_index, true, true)

end

pars_to_criterion(x1)



x1=pars_to_vec(pars)
Criterion(x0,pars,Mod1,vlist,moms0,wghts_Full,5,lengths,TE_index, true, true)

Criterion(x1,pars,Mod1,vlist,moms0,wghts_Full,5,lengths,TE_index, true, true)



#ForwardDiff.gradient(pars_to_criterion,x1)

# automatic differentiation with dual numbers does not work due to type issue--is this fixable?

println("Just getting Labor Parameters")

opt,x0 = GetOptimization(pars,Mod1,labor_block,moms0,wghts_labor,5,lengths,TE_index,maxevals = 200,LBFGS = 1)
res = NLopt.optimize(opt,x0)




# finally let's polish off with a local optimizer or two

opt,x0 = GetOptimization(pars,Mod1,labor_block,moms0,wghts_labor,5,lengths,TE_index,maxevals = 1000,SBPLX = 1)
res2 = NLopt.optimize(opt,x0)

opt,x0 = GetOptimization(pars,Mod1,labor_block,moms0,wghts_labor,5,lengths,TE_index,maxevals = 1000)
res3 = NLopt.optimize(opt,x0)

# OK, Nelder-Mead seemed to do better than Subplex here, and it didn't seem like it was done optimizing.

opt,x0 = GetOptimization(pars,Mod1,labor_block,moms0,wghts_labor,5,lengths,TE_index,maxevals = 1000)
res4 = NLopt.optimize(opt,x0)

# the other nelder mead code is helpful as I can control the simplex
# but it doesn't support bounds and should only be used for polishing

A=Optim.optimize(pars_to_criterion,x1,NelderMead(;initial_simplex=Optim.AffineSimplexer(0.0, 0.01)), Optim.Options(iterations=200))


# Next let's check and store the pars

pars
x0


x1=pars_to_vec(pars)

Criterion(x1,pars,Mod1,vlist,moms0,wghts_Full,5,lengths,TE_index, true, true)
Criterion(x1,pars,Mod1,vlist,moms0,wghts,5,lengths,TE_index, true, true)
Criterion(x1,pars,Mod1,vlist,moms0,wghts_labor,5,lengths,TE_index, true, true)

DF_Pars=DataFrame([Names_Pars,x1])

CSV.write("Latest_Params.csv",DF_Pars)


# Next let's store all the model and simulation data in one place

E,A,A2,XG,skill,Y =MomentsBaseline(Mod1,5,lengths,TE_index)

Mom_sims=[E;A;A2;XG;skill;Y]
Mom_delta=(Mom_sims.-moms0).^2
Mom_delta_wght=Mom_delta.*wghts_labor

DFNames=[:Moment_Category,:Data_Moment,:Model_Moment,:weights,:Moment_Deviation,:Moment_Deviation_Weighted]

DF2=DataFrame([Names_all, moms0,Mom_sims,wghts_labor,Mom_delta,Mom_delta_wght])
names!(DF2,[:Moment_Category,:Data_Moment,:Model_Moment,:weights,:Moment_Deviation,:Moment_Deviation_Weighted])
DF2

CSV.write("Moment_Fit.csv",DF2)

clf()
X2=Fit(Mod1)
gcf()
savefig("LFP_Part_Fit.pdf")

clf()
X2=Fit_A2_Y(Mod1)
gcf()
savefig("TotInc_Receipt_Fit.pdf")

#opt,x0 = GetOptimization(pars,Mod1,labor_block,moms0,wghts_labor,5,lengths,TE_index,maxevals = 500,LBFGS = 0)
#res = optimize(opt,x0)


#=

Reminders:

E is LFP (aka L)
A is program participation
A2 is program receipt
Xg is govt expenditure
Y is earnings
Te is treatment effects on children


=#
