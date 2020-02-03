# this assumes we have run the script "SetupBaseline.jl"
include("ProductionEstimation.jl")
#measures = [:Achievement,:AchieveBelowAverage,:Repeat,:Math,:BelowMath,:Read,:AboveRead,:PB,:BPI]
#measures = [:Achievement,:AchieveBelowAverage,:Repeat,:BelowRead]
#measures = [:Achievement,:AchieveBelowAverage,:Repeat,:PB,:BPI,:Math,:BelowMath,:Read,:AboveRead]
#measures = [:Achievement,:AchieveBelowAverage,:Repeat]
#measures = [:AchieveBelowAverage,:PB,:BPI,:Math,:BelowMath,:Read,:BelowRead,:Repeat]
#measures = [:Achievement,:AchieveBelowAverage,:Math,:BelowMath,:Read,:BelowRead]
#measures = [:PB,:Repeat]#,:Read,:BelowRead] #,:BPI]
measures = [:Achievement,:Repeat] #,:Math]#,:Math,:Read] #<- could we possibly get more data here?
measures = [:BPI,:PB] #,:Suspend]
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
    for i=1:length(measures)
        TE[:,i] = convert(Array{Float64,1},coalesce.(TE_[measures[i]],0.) ./se[measures[i]])
        wght[Idrop[i],i] .= 0.
    end
    if s.==:MFIPRA
        wght .= 0.
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
vlist = [:δI] #,:gN,:gF] #,:gN,:gF] #,:λ]
opt,x0 = GetOptimization(vlist,pars_prod,pars,CP,site_list,budget,TEmoms,site_features)

res = optimize(opt,x0)
pars_prod = UpdatePars(res[2],pars_prod,vlist)

V = GetVariance(res[2],vlist,pars_prod,pars,CP,site_list,budget,TEmoms,site_features)

se = sqrt.(diag(V))

println(vlist)
display([res[2] .- 1.96*se res[2] res[2] .+ 1.96*se])

T = SE[:,measures]
for v in measures
    T[v] = SE[v] .* sqrt.(4 ./(D.N_control+D.N_treat))
end

break
V = Symmetric(V)
wage = 7.50
B = zeros(100)
for b=1:100
    p = copy(res[2])
    p[3:end] .+= rand(MultivariateNormal(V[3:end,3:end]))
    pp = UpdatePars(p,pars_prod,vlist)
    th = 0
    for t=1:3
        th = pp.δI[1]*log((30*wage + (112-30)*pars2.wq)/(112*pars2.wq)) - pp.gN[1] + pp.δθ*th
    end
    B[b] = th
end



break

B = zeros(2,100)
for b=1:100
    moms_collect = []
    for s in site_list
        m = D[D.Site.==String(s),:]
        se = SE[SE.Site.==String(s),:]
        TE_ = m[:,measures]
        se_ = se[:,measures]
        tevar = convert(Array{Float64,2},se[:,measures]).^2
        N = convert(Array{Float64,1},m.N_control+m.N_treat)
        tevar = 2*tevar./N
        wght = 1 ./tevar
        Idrop = [ismissing.(TE_[v]) for v in measures] #convert(Array{Bool,2},ismissing.(TE))
        TE = zeros(size(TE_)) #convert(Array{Float64,2},coalesce.(TE,0.))
        for i=1:length(measures)
            for k=1:size(TE_)[1]
                if ~ismissing(TE_[k,i])
                    TE[k,i] = TE_[k,i] + rand(Normal(0,sqrt(tevar[k,i])))
                end
            end
            #TE[:,i] = convert(Array{Float64,1},coalesce.(TE_[measures[i]],0.))
            #wght[Idrop[i],i] .= 0.
        end
        arm = convert(Array{Int64,1},m.Treatment) .+ 1
        a0 = convert(Array{Int64,1},m.AgeMin)
        a1 = convert(Array{Int64,1},m.AgeMax)
        append!(moms_collect,[(TE=TE,wght=wght,arm=arm,a0=a0,a1=a1)])
    end
    TEmoms = (;zip(site_list,moms_collect)...)
    opt,x0 = GetOptimization(vlist,pars_prod,pars,CP,site_list,budget,TEmoms,site_features)
    res2 = optimize(opt,x0)
    B[:,b] = res2[2]
end



break


TEmod = GetTreatmentEffects(pars_prod,pars2,CP,site_list,budget,TEmoms,site_features)

# opt,x0 = GetOptimization([:λ],pars_prod,pars,CP,site_list,budget,TEmoms,site_features)
#
# res = optimize(opt,x0)
# pars_prod = UpdatePars(res[2],pars_prod,[:λ])



# opt,x0 = GetOptimization([:λ],pars_prod,pars2,CP,site_list,budget,TEmoms,site_features)
# res2 = optimize(opt,x0)
# pars_prod = UpdatePars(res2[2],pars_prod,[:λ])
InspectTreatFitProduction!(pars_prod,pars,CP,site_list,budget,TEmoms,site_features)

break

LFP = []
LFP2 = []
INC = []
INC2 = []
TE = []
TEm = []
for i=1:8
    sname = site_list[i]
    println(sname)
    m = data_moments[sname]
    m2 = moms[sname]
    for a=2:site_features.n_arms[i]
        ii = TEmoms[sname].arm.==a
        figure("LFP")
        lfp = (mean(m[a].LFP) - mean(m[1].LFP))*ones(sum(ii))
        lfp2 = (mean(m2[a].LFP) - mean(m2[1].LFP))*ones(sum(ii))
        te = TEmoms[sname].TE[ii,1]
        tem = TEmod[sname][ii]
        scatter(lfp,te,color="blue")
        scatter(lfp2,te,color="red")
        scatter(lfp2,tem,color="pink")
        figure("Income")
        inc = (mean(m[a].Inc) - mean(m[1].Inc))*ones(sum(ii))
        inc2 =(mean(m2[a].Inc) - mean(m2[1].Inc) - lfp2[1]*pars.wq*30)*ones(sum(ii))
        scatter(inc,te,color="blue")
        scatter(inc2,te,color="red")
        scatter(inc2,tem,color="orange")
        append!(LFP,lfp)
        append!(LFP2,lfp)
        append!(INC,inc)
        append!(INC2,inc2)
        append!(TE,te)
    end
end

B = zeros(2,100)
for b=1:100
    LFP = []
    LFP2 = []
    INC = []
    INC2 = []
    TE = []
    TEm = []
    WGHT = []
    for i=1:8
        sname = site_list[i]
        #println(sname)
        m = data_moments[sname]
        m2 = moms[sname]
        for a=2:site_features.n_arms[i]
            ii = TEmoms[sname].arm.==a
            wght = TEmoms[sname].wght[ii,1]
            lfp = (mean(m[a].LFP) - mean(m[1].LFP))*ones(sum(ii))
            lfp2 = (mean(m2[a].LFP) - mean(m2[1].LFP))*ones(sum(ii))
            te = TEmoms[sname].TE[ii,1]
            for t=1:length(te)
                if wght[t]>0
                    se = sqrt(1/wght[t])
                    te[t] += rand(Normal(0.,se))
                end
            end
            tem = TEmod[sname][ii]
            inc = (mean(m[a].Inc) - mean(m[1].Inc))*ones(sum(ii))
            inc2 =(mean(m2[a].Inc) - mean(m2[1].Inc) - lfp2[1]*pars.wq*30)*ones(sum(ii))
            append!(LFP,lfp)
            append!(LFP2,lfp)
            append!(INC,inc)
            append!(INC2,inc2)
            append!(TE,te)
            append!(WGHT,wght)
        end
        #scatter(INC[WGHT.>0]./LFP[WGHT.>0],TE[WGHT.>0])
        B[1,b] = sum(TE.*INC2.*WGHT)/sum(WGHT)
        B[2,b] = 100*sum(TE.*LFP.*WGHT)/sum(WGHT)
    end
end

meas = [:Achievement,:AchieveBelowAverage,:Repeat,:Suspend,:Math,:Read,:BPI,:PB]


for i=1:8
    sname = site_list[i]
    println(sname)
    m = data_moments[sname]
    m2 = moms[sname]
    for a=2:site_features.n_arms[i]
        ii = TEmoms[sname].arm.==a
        figure("LFP")
        lfp = mean(m[a].LFP) - mean(m[1].LFP)
        lfp2 = mean(m2[a].LFP) - mean(m2[1].LFP)
        scatter(lfp,mean(TEmoms[sname].TE[ii,1]),color="blue")
        scatter(lfp2,mean(TEmoms[sname].TE[ii,1]),color="red")
        figure("Income")
        inc = mean(m[a].Inc) - mean(m[1].Inc)
        inc2 = mean(m2[a].Inc) - mean(m2[1].Inc) - lfp2*pars.wq*30
        scatter(inc,mean(TEmoms[sname].TE[ii,1]),color="blue")
        scatter(inc2,mean(TEmoms[sname].TE[ii,1]),color="red")
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
