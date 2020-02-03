# this assumes we have run the script "SetupBaseline.jl"
include("ProductionEstimation.jl")
include("MCMCRoutines.jl")
#measures = [:PB,:Repeat]#,:Read,:BelowRead] #,:BPI]
measures = [:Achievement,:Repeat,:Math,:Read] #,:Math]#,:Math,:Read] #<- could we possibly get more data here?
#measures = [:BPI,:PB] #,:Suspend]
#measures = [:Achievement,:Math]

D = CSV.read("../Data/ChildTreatmentEffects.csv")
SE = CSV.read("../Data/ChildTreatmentEffectsSEs.csv")
#D.Achievement = D.Achievement
#D.Math = D.Math*10/17 #<- Table 12.3 bounds standard error for WJ math between 15 and 20
D.AchieveBelowAverage = -D.AchieveBelowAverage
D.Suspend = -D.Suspend
D.BelowMath = -D.BelowMath
D.BelowRead = -D.BelowRead
D.Repeat = -D.Repeat
D.BPI = -D.BPI

D.N_treat = coalesce.(D.N_treat,0)
for v in measures; SE[v] = coalesce.(SE[v],Inf); end

# create a TEmoms object that stores all the info we want
moms_collect = []
for s in site_list
    m = D[D.Site.==String(s),:]
    se = SE[SE.Site.==String(s),:]
    TE_ = m[:,measures]
    se_ = se[:,measures]
    tevar = convert(Array{Float64,2},se[:,measures]).^2
    N = convert(Array{Float64,1},m.N_control+m.N_treat)
    tevar = 4 ./N
    wght = ones(1,length(measures)) ./tevar
    Idrop = [ismissing.(TE_[v])[:] for v in measures] #convert(Array{Bool,2},ismissing.(TE))
    TE = zeros(size(TE_)) #convert(Array{Float64,2},coalesce.(TE,0.))
    SE_ = zeros(size(TE_))
    for i=1:length(measures)
        TE[:,i] = convert(Array{Float64,1},coalesce.(TE_[measures[i]],0.) ./se[measures[i]])
        wght[Idrop[i],i] .= 0.
        SE_[:,i] = sqrt.(4 ./N)
    end
    if s.==:MFIPRA
        wght .= 0.
    end
    arm = convert(Array{Int64,1},m.Treatment) .+ 1
    a0 = convert(Array{Int64,1},m.AgeMin)
    a1 = convert(Array{Int64,1},m.AgeMax)
    append!(moms_collect,[(TE=TE,SE=SE_,wght=wght,arm=arm,a0=a0,a1=a1)])
end


TEmoms = (;zip(site_list,moms_collect)...)

# Initialize the parameters

x0,pars_bayes = ParsBayes(length(measures))
#pars_prod = ProdPars(length(measures))

CP = GetChoiceProbsAll(pars2,site_list,budget,site_features);

#ProductionCriterion(pars_prod,pars2,CP,site_list,budget,TEmoms,site_features)
#vlist = [:gN,:gF,:δI,:δθ]
vlist = [:δI,:δθ,:gN,:gF,:λ] #,:gN,:gF] #,:λ]
#opt,x0 = GetOptimization(vlist,pars_prod,pars,CP,site_list,budget,TEmoms,site_features)
P,Ahist,Lhist = GetChain(5000,x0,vlist,pars_bayes,CP,site_list,budget,TEmoms,site_features)

for i=1:8
    subplot(2,4,i)
    hist(P[i,2000:end],density=true)
end
