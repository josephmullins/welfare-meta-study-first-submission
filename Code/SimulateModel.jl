using Revise
using Distributions
using Random
using Distributed
using SharedArrays
using Statistics
using StatsBase


cd("/Users/FilipB/github/welfare-meta-study/Code")
includet("BaselineModel.jl")
includet("EstimationBaseline.jl")
includet("SetupBaseline.jl")



Γ = exp.(0 .+ -0.1*(1:18))


best_pars=readdlm("current_ests")
pars = parameters()
vlist = [:αc, :αH, :αA, :σH, :σC, :αWR,:αWR2, :αF, :β,:wq]

pars = UpdatePars(best_pars,pars,vlist)
site_list
sample_size=[4803, 1405+1410,15683,3208,6009,4433,4554,8322]

site_features.yb
budget.CTJF[1,:,:,:,:]
1

# This block of code sets up some things I use to test the functions below
αHT = [0;pars.αHT]
yb = site_features.yb[3]
T = site_features.T[3]
years = (yb+1-1991):(yb-1991+T)
pos = sum(site_features.T[1:1-1])
#αH = pars.αH[(pos+1):(pos+site_features.T[i])]
pars_site = GetSitePars(pars,3)
sname = site_list[3]
Y = getfield(budget,sname)
moms = getfield(moments,sname)
#π0_a = site_features.π0[2,2,:]
#sum(π0_a)
#π01a = site_features.π0[:,:,:]
#π01 = site_features.π0[1,:,:]
year_meas = site_features.year_meas[3]
Abar = size(Y)[3]
moms_model = zeros(size(moms))
WR = site_features.work_reqs[3,2]
price = site_features.prices[3,2]
abar = min(Abar,2)
#Y_I = budget[Symbol(sname,"_I")]
TLlength = site_features.TLlength[3,2]
NK =3# size(π01)[3]
pA1,pWork1,pF1 = GetStaticProbs(pars_site,Y[abar,:,:,:,:],price,WR,T,NK)
#pA2,pWork2,pF2 = GetDynamicProbs(pars_site,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,NK,TLlength)
# end test code prep



#=

CTJF: 4969
    - CTJF samples ALL kids at last year

FTP: 3698/1108
    - FTP gives choice of last year for all kids or all sample for focal kid

MFIPLR: 879
    - Seems to be for all kids, not at last year
MFIPRA: 587
    - Does not seem to list sample sizes?

LAGAIN: 372+374
    - Asks if EVER used child care


NEWWSA: 580/650 LFA/HCD
NEWWSG: 394/380 LFA/HCD
NEWWSR: 492/405 LFA/HCD, 312 LFA no degree

    - All NEWWS were in past week for focal child

None of these specified less than 9!
=#

#kid_sample_size=[4969, ]

(π0[1,:,:]*ones(17))'*[1,2,3]

mean_nk=zeros(8)

for i in 1:8
    mean_nk[i]=(π0[i,:,:]*ones(17))'*[1,2,3]
end

C = CSV.read("../Data/ChildCareMoms_estimated.csv")

formal_sample_size=zeros(8,3)

site_str = ["CTJF","FTP","LAGAIN","MFIP-LR","MFIP-RA","NEWWS-A","NEWWS-G","NEWWS-R"]

for i in 1:8
    for j in 1:(site_features.n_arms[i])
        c = C[(C.Site.==site_str[i]) .& (C.Arm.==j-1),:]
        formal_sample_size[i,j]=c.N[1]
    end
end

formal_sample_size[1,3]
1

 # old was [4803,1405+1410,3208,6009]

 # I got the numbers for NEWWS from line 139, and for LAGAIN from 374




function draw_kids(site,site_features)
    # get NK
    π0_nk = site_features.π0[site,:,:]
    prob_nk=zeros(3)
    prob_nk[1]=sum(π0_nk[1,:])
    for i in 2:3
        prob_nk[i]=sum(π0_nk[i,:])+prob_nk[i-1]
    end
    nk_draw=rand(Uniform())
    nk=0

    A1=zeros(3)
    for i in 1:2
        A1[i]=ifelse(nk_draw>=prob_nk[i],1.0,0)
    end
    nk=convert(Int,sum(A1))+1



    # get age
    π0_a = site_features.π0[site,nk,:]
    if sum(π0_a)>1.0001
        println("problem: pi")
    end
    π0_a=π0_a./sum(π0_a)
    prob_age=zeros(17)
    prob_age[1]=π0_a[1]
    for i in 2:17
        prob_age[i]=π0_a[i]+prob_age[i-1]
    end

    age_draw=rand(Uniform())
    A2=zeros(17)
    for i in 1:16 # will never be above 16 by construction
        A2[i]=ifelse(age_draw>=prob_age[i],1.0,0) # A2[i]=1 if weakly older than age i
    end

    age=convert(Int,sum(A2))+1 # sum A2 vector to get age INDEX

    return nk, age
end
1


function decisions_static(t,nk, age,pA,pWork,pF,budget_mat)

    A=0
    W=0
    F=0

    part=rand(Uniform())
    if part<=pA[t, nk, age]
        A=1
    end

    working=rand(Uniform())
    if working<=pWork[t,nk, age, A+1]
        W=1
    end

    formal=rand(Uniform())
    if formal<=pF[t,nk, age, A+1] && W==1
        F=1
    end

    Y=budget_mat[t, nk,A+1,W+1] # pretty sure it's time, numkids, part, work

    return A,W,F, Y
end
1

function decisions_dynamic(t,nk,age,tl,pA,pWork, pF,budget_mat)
    A=0
    W=0
    F=0

    part=rand(Uniform())
    if part<=pA[t, nk, age,tl]
        A=1
    end

    working=rand(Uniform())
    if working<=pWork[t,nk, age, tl,A+1]
        W=1
    end

    formal=rand(Uniform())
    if formal<=pF[t,nk, age, tl,A+1] && W==1
        F=1
    end

    Y=budget_mat[t, nk,A+1,W+1] # pretty sure it's time, numkids, part, work

    return A,W,F,Y
end
1


function Moments_Simulated_Static(nsims, pA,pWork,pF, sample_size_site,T, site,site_features,kid_length, budget_mat)

    A_mat=zeros(T,nsims, sample_size_site)
    A_sims=zeros(T,nsims)
    A=zeros(T)

    Y_mat=zeros(T,nsims, sample_size_site)
    Y_sims=zeros(T,nsims)
    Y=zeros(T)


    W_mat=zeros(T,nsims, sample_size_site)
    W_mat_2=zeros(T,nsims, sample_size_site) # for childcare calculation
    W_sims=zeros(T,nsims)
    W=zeros(T)


    F_mat=zeros(T,nsims, sample_size_site)



    # note: only give benefits if kids aged 0-9, aka indices 1-10
    @inbounds @simd for n in 1:nsims
        @inbounds @simd  for i in 1:sample_size_site
                                nk, age=draw_kids(site,site_features)
                                a=age
                                    @inbounds  for t in 1:T
                                                    A_mat[t,n,i],W_mat[t,n,i],F0,Y_mat[t,n,i]=decisions_static(t,nk, age,pA,pWork,pF,budget_mat)

                                                        if (a>=1) & (a<=10) &  (n<=kid_length)
                                                            F_mat[t,n,i]=F0
                                                            W_mat_2[t,n,i]=W_mat[t,n,i]
                                                        end
                                                    #a+=1
                                                end
        end
    end
    # aggregate the moments
    @inbounds @simd for t in 1:T
        @inbounds @simd for n in 1:nsims
            A_sims[t,n]=mean(A_mat[t,n,:])
            W_sims[t,n]=mean(W_mat[t,n,:])
            Y_sims[t,n]=mean(Y_mat[t,n,:])
        end
            A[t]=mean(A_sims[t,:])
            W[t]=mean(W_sims[t,:])
            Y[t]=mean(Y_sims[t,:])
    end
    F=sum(F_mat)/sum(W_mat_2)

    return A,W,F,Y
end
1

function Moments_Simulated_Dynamic(nsims, pA,pWork,pF, sample_size_site,T, site, tlsite,site_features,kid_length, budget_mat, budget_mat_inel)
    A_mat=zeros(T,nsims, sample_size_site)
    A_sims=zeros(T,nsims)
    A=zeros(T)

    Y_mat=zeros(T,nsims, sample_size_site)
    Y_sims=zeros(T,nsims)
    Y=zeros(T)

    W_mat=zeros(T,nsims, sample_size_site)
    W_mat_2=zeros(T,nsims, sample_size_site)
    W_sims=zeros(T,nsims)
    W=zeros(T)

    F_mat=zeros(T,nsims, sample_size_site)
    F0=0

    # child care only for kids in ages 0:9
    @inbounds @simd for n in 1:nsims
         @inbounds @simd for i in 1:sample_size_site
                            tl=1
                            nk, age=draw_kids(site,site_features)
                            a=age
                                @inbounds for t in 1:T
                                    if tl==tlsite+1
                                        A_mat[t,n,i],W_mat[t,n,i],F0,Y_mat[t,n,i]=decisions_dynamic(t,nk,age,tl,pA,pWork, pF, budget_mat_inel)

                                    else
                                        A_mat[t,n,i],W_mat[t,n,i],F0,Y_mat[t,n,i]=decisions_dynamic(t,nk,age,tl,pA,pWork, pF, budget_mat)

                                    end
                                            if (a>=1) & (a<=10) & (n<=kid_length)
                                                    F_mat[t,n,i]=F0
                                                    W_mat_2[t,n,i]=W_mat[t,n,i]
                                                end
                                            if A_mat[t,n,i]==1
                                                tl=min(tl+1,tlsite+1)
                                            end
                                            #a+=1 # age 0-9 at study ENTRY
                                        end
        end
    end
    @inbounds @simd for t in 1:T
        @inbounds @simd for n in 1:nsims
                            A_sims[t,n]=mean(A_mat[t,n,:])
                            W_sims[t,n]=mean(W_mat[t,n,:])
                            Y_sims[t,n]=mean(Y_mat[t,n,:])
                        end
                            A[t]=mean(A_sims[t,:])
                            W[t]=mean(W_sims[t,:])
                            Y[t]=mean(Y_sims[t,:])
                    end

    F=sum(F_mat)/sum(W_mat_2)

    return A,W,F,Y

end
1

function Get_Simulated_MomentsAll(pars1,site_list,budget,moments,wghts,site_features,nsims,sample_size,formal_sample_size; seed=1234, wrapped=false)
    Random.seed!(seed)
    moms_collect = []
    Qn = 0
    αHT = [0;pars1.αHT]
    for i=1:length(site_list)
        NK=3
        yb = site_features.yb[i]
        T = site_features.T[i]
        years = (yb+1-1991):(yb-1991+T)
        pos = sum(site_features.T[1:i-1])
        #αH = pars.αH[(pos+1):(pos+site_features.T[i])]
        pars_site = GetSitePars(pars1,i)
        sname = site_list[i]
        Y = getfield(budget,sname)
        moms = getfield(moments,sname)
        #π1 = site_features.π0[2,:,:]
        #π0 = site_features.π0
        wght = getfield(wghts,sname)
        year_meas = site_features.year_meas[i]
        Abar = size(Y)[1]
        moms_model = zeros(size(moms))
        arms=site_features.n_arms[i]
        ssize2=convert(Int,(round(sample_size[i]/arms))) # will replace later with actual arm sample sizes
        if wrapped==true
            moms_model=[]
        end
        price = site_features.prices[i,1]
        WR = site_features.work_reqs[i,1]
        pA1,pWork1,pF1 = GetStaticProbs(pars_site,Y[1,:,:,:,:],price,WR,T,NK) # last 3 are WR, T, NK


        A,W,F,Inc=Moments_Simulated_Static(nsims, pA1,pWork1,pF1,  ssize2,T, i,site_features, formal_sample_size[i,1],Y[1,:,:,:,:])

        moms_control = [A;W;F;Inc]
        Qn += sum(wght[:,1].*(moms[:,1] .- moms_control).^2)

        if wrapped==true
            append!(moms_model, [(Part=A,LFP=W,Inc=Inc,Care=F)])
        else
            moms_model[:,1]=moms_control
        end


        for a = 2:site_features.n_arms[i]
            WR = site_features.work_reqs[i,a]
            price = site_features.prices[i,a]
            abar = min(Abar,a)
            ssize2=convert(Int,(round(sample_size[i])))
            if site_features.time_limits[i,a]==1
                Y_I = budget[Symbol(sname,"_I")]
                TLlength = site_features.TLlength[i,a]
                pA2,pWork2,pF2 = GetDynamicProbs(pars_site,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,NK,TLlength)
                println(nsims)
                A,W,F,Inc=Moments_Simulated_Dynamic(nsims,pA2,pWork2,pF2,  ssize2,T,i,TLlength,site_features, formal_sample_size[i,abar],Y[abar,:,:,:,:],Y_I[abar,:,:,:,:])#(nsims, pA,pF,pWork, π1, sample_size_site,T, site, tlsite,site_features)
                moms_treatment = [A;W;F;Inc]

                if wrapped==true
                    append!(moms_model, [(Part=A,LFP=W,Inc=Inc,Care=F)])
                else
                    moms_model[:,a]=moms_treatment
                end
            else
                pA3,pWork3,pF3 = GetStaticProbs(pars_site,Y[abar,:,:,:,:],price,WR,T,NK)
                A,W,F,Inc=Moments_Simulated_Static(nsims,pA3,pWork3,pF3, ssize2,T, i,site_features, formal_sample_size[i,a],Y[abar,:,:,:,:])
                moms_treatment = [A;W;F;Inc]
                if wrapped==true
                    append!(moms_model, [(Part=A,LFP=W,Inc=Inc,Care=F)])
                else
                    moms_model[:,a]=moms_treatment
                end
            end
            TE = moms[:,a] .- moms[:,1]
            TE_mod = moms_treatment .- moms_control
            Qn += sum(wght[:,a].*(TE .- TE_mod).^2)
        end
        append!(moms_collect,[moms_model])
    end
    return Qn, (;zip(site_list,moms_collect)...)
end

fit_sim, moms_sim_optimized=
    Get_Simulated_MomentsAll(pars,site_list,budget,moments,wghts,
    site_features,1,sample_size, formal_sample_size; wrapped=false)

CriterionP(pars,site_list,budget,moments,wghts,site_features)

# check with optimized params
moments.CTJF

opt,x0 = GetOptimization(pars,vlist,site_list,budget,moments,wghts,site_features)
res2 = optimize(opt,x0)
pars2 = UpdatePars(res2[2],pars,vlist)
fit_sim, moms_sim_optimized=
    Get_Simulated_MomentsAll(pars2,site_list,budget,moments,wghts,
    site_features,5,sample_size, formal_sample_size; wrapped=false)




CriterionP(pars2,site_list,budget,moments,wghts,site_features)
moms_diff_optimized=GetMomentsAll(pars2,site_list,budget,moments,wghts,site_features)

pars.αc
pars2.αc



#=

bootstrapped SE code

=#


length1=1
success_vec=[]
param_mat=zeros(length(res2[2]), length1)

@time for s in 1:length1
    fit_sim, moms_sim_optimized_bootstrap=Get_Simulated_MomentsAll(pars2,site_list,budget,moments,wghts,site_features,1,sample_size, formal_sample_size; wrapped=false, seed=s)
    opt2,x2 = GetOptimization(pars2,vlist,site_list,budget,moms_sim_optimized_bootstrap,wghts,site_features)
    @time res3 = optimize(opt2,x2)
    println(res3[3])
    param_mat[:,s].=res3[2]
    append!(success_vec,[res3[3]])
end

param_mat


mean_sim_params=zeros(length(param_mat[:,1]))
ses=zeros(length(param_mat[:,1]))
for i in 1:length(param_mat[:,1])
    mean_sim_params[i]=mean(param_mat[i,1:10])
    ses[i]=std(param_mat[i,1:10])
end





names1=[]
names1_index=[]
names1_number=[]
for i in 1:length(vlist)
    l=pars.np[vlist[i]]
    for j in 1:l
        append!(names1,[string(vlist[i])])
        append!(names1_index,j)
        append!(names1_number,i)
    end
end


df = DataFrame()
df.names=names1
df.par_number=names1_number
df.names_index=names1_index
df.estimated_params=res2[2]
df.baseline_params=best_pars[:]
df.mean_sim_params=mean_sim_params
df.standard_deviations=ses


for i in 1:length1
    df.i=param_mat[:,i]
    rename!(df, Dict(:i =>"sim_$i"))
end
*CSV.write("Bootstrap/Bootstrap_params.csv",df)
*df2=CSV.read("Bootstrap/Bootstrap_params.csv")




function InspectModelFit2(model_moments,data_moments,site_features,site_list;Part=true, inc=false)
    colors = ["blue","green","red"]
    for i=1:8
        T = site_features.T[i]
        sname = site_list[i]
        moms0 = getfield(data_moments,sname)
        moms1 = getfield(model_moments,sname)
        for a=1:site_features.n_arms[i]
            v=ifelse(Part==false,:LFP,:Part)
            v=ifelse(inc==true,:Inc,:Part)
                figure(String(v))
                subplot(2,4,i)
                PyPlot.title(String(sname))
                PyPlot.plot(moms0[a][v],color=colors[a])
                PyPlot.plot(moms1[a][v],color=colors[a],linestyle="--")
            #end
        end
    end
end


fit_sim, moms_sim_optimized=Get_Simulated_MomentsAll(pars2,site_list,budget,moments,wghts,site_features,1,sample_size, formal_sample_size; wrapped=true)
CriterionP(pars2,site_list,budget,moments,wghts,site_features)
moms_diff_optimized=GetMomentsAll(pars2,site_list,budget,moments,wghts,site_features)


InspectModelFit2(moms_diff_optimized,data_moments,site_features,site_list)
# compare participation
clf()
Participation_Diff=InspectModelFit2(moms_diff_optimized,data_moments,site_features,site_list)
gcf()
clf()
Participation_Sim=InspectModelFit2(moms_sim_optimized,data_moments,site_features,site_list)
gcf()

# compare lfp
clf()
LFP_Diff=InspectModelFit2(moms_diff_optimized,data_moments,site_features,site_list, Part=false)
gcf()
clf()
LFP_Sim=InspectModelFit2(moms_sim_optimized,data_moments,site_features,site_list, Part=false)
gcf()


# compare inc
clf()
LFP_Diff=InspectModelFit2(moms_diff_optimized,data_moments,site_features,site_list, inc=true)
gcf()
clf()
LFP_Sim=InspectModelFit2(moms_sim_optimized,data_moments,site_features,site_list, inc=true)
gcf()

moms_diff_optimized.MFIPLR[1].Care






moms_sim_optimized.MFIPLR[1].Care
