# this assumes we have run the script "SetupBaseline.jl"
include("ProductionEstimation.jl")
#measures = [:Achievement,:AchieveBelowAverage,:Repeat,:Math,:BelowMath,:Read,:AboveRead]
#measures = [:Achievement,:AchieveBelowAverage,:Repeat,:BelowRead]
#measures = [:Achievement,:AchieveBelowAverage,:Repeat,:PB,:BPI,:Math,:BelowMath,:Read,:AboveRead]
#measures = [:Achievement,:AchieveBelowAverage,:Repeat]
#measures = [:AchieveBelowAverage,:PB,:BPI,:Math,:BelowMath,:Read,:BelowRead,:Repeat]
measures = [:Achievement,:AchieveBelowAverage,:Math,:BelowMath,:Read,:BelowRead]
D = CSV.read("../Data/ChildTreatmentEffects.csv")
D.Achievement = D.Achievement*10
D.AchieveBelowAverage = -D.AchieveBelowAverage
D.BelowMath = -D.BelowMath
D.BelowRead = -D.BelowRead

D.N_treat = coalesce.(D.N_treat,0)


# create a TEmoms object that stores all the info we want
moms_collect = []
for s in site_list
    m = D[D.Site.==String(s),:]
    TE_ = m[:,measures]
    wght = convert(Array{Float64,1},m.N_control+m.N_treat)/1000
    wght = wght.*ones(1,size(TE_)[2])
    Idrop = [ismissing.(TE_[v]) for v in measures] #convert(Array{Bool,2},ismissing.(TE))
    TE = zeros(size(TE_)) #convert(Array{Float64,2},coalesce.(TE,0.))
    for i=1:length(measures)
        TE[:,i] = convert(Array{Float64,1},coalesce.(TE_[measures[i]],0.))
        wght[Idrop[i],i] .= 0.
    end
    arm = convert(Array{Int64,1},m.Treatment) .+ 1
    a0 = convert(Array{Int64,1},m.AgeMin)
    a1 = convert(Array{Int64,1},m.AgeMax)
    append!(moms_collect,[(TE=TE,wght=wght,arm=arm,a0=a0,a1=a1)])
end
TEmoms = (;zip(site_list,moms_collect)...)

# Initialize the parameters


pars_prod = ProdPars(length(measures))

CP = GetChoiceProbsAll(pars2,site_list,budget,site_features);

#ProductionCriterion(pars_prod,pars2,CP,site_list,budget,TEmoms,site_features)
#vlist = [:gN,:gF,:δI,:δθ]
vlist = [:δI,:δθ,:λ]
opt,x0 = GetOptimization(vlist,pars_prod,pars,CP,site_list,budget,TEmoms,site_features)

res = optimize(opt,x0)
pars_prod = UpdatePars(res[2],pars_prod,vlist)




# opt,x0 = GetOptimization([:λ],pars_prod,pars2,CP,site_list,budget,TEmoms,site_features)
# res2 = optimize(opt,x0)
# pars_prod = UpdatePars(res2[2],pars_prod,[:λ])
InspectTreatFitProduction!(pars_prod,pars2,CP,site_list,budget,TEmoms,site_features)

break

for i=1:8
    sname = site_list[i]
    m = data_moments[sname]
    for a=2:site_features.n_arms[i]
        figure("Income")
        inc = mean(m[a].Inc) - mean(m[1].Inc)
        scatter(inc*ones(size(TEmoms[sname].TE,1)),TEmoms[sname].TE[:,a-1])
        figure("LFP")
        inc = mean(m[a].LFP) - mean(m[1].LFP)
        scatter(inc*ones(size(TEmoms[sname].TE,1)),TEmoms[sname].TE[:,a-1])
    end
end




break
i=1
sname = site_list[i]
#println(sname)
moms = getfield(TEmoms,sname) #<- structure: age0,age1,1
nmom = size(moms.TE)[1]
π0 = site_features.π0[i,:,:]
#wght = getfield(wghts,sname)
year_meas = site_features.year_meas[i]
Y = getfield(budget,sname)
P = getfield(CP,sname)
TH = zeros(Real,site_features.n_arms[i],NK,17)
price = site_features.prices[i,1]
TH[1,:,:] =  GetChildOutcomesStatic(year_meas,P.pA[1],P.pWork[1],P.pF[1],Y[1,:,:,:,:],pars_prod,price,0.1)

Abar = size(Y)[1]
#moms_control = GetMoments(pars_site,Y[1,:,:,:,:],price,0,T,π0,year_meas)
#Qn += sum(wght[:,1].*(moms[:,1] .- moms_control).^2)
for a = 2:site_features.n_arms[i]
    WR = site_features.work_reqs[i,a]s
    price = site_features.prices[i,a]
    abar = min(Abar,a)
    if site_features.time_limits[i,a]==1
        Y_I = budget[Symbol(sname,"_I")]
        TLlength = site_features.TLlength[i,a]
        TH[a,:,:] = GetChildOutcomesDynamic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],pars_prod,price,0.1)
    else
        TH[a,:,:] = GetChildOutcomesStatic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[abar,:,:,:,:],pars_prod,price,0.1)
    end
end
TE = TreatmenEffects(TH,π0,moms.a0,moms.a1,moms.arm)
