using Revise
using Distributions
using Random
using Distributed
using SharedArrays

nprocs()

cd("/Users/FilipB/github/welfare-meta-study/Code")
includet("BaselineModel.jl")
includet("EstimationBaseline.jl")
includet("SetupBaseline.jl")

Γ = exp.(0 .+ -0.1*(1:18))
pars = parameters()
vlist = [:αc, :αH, :αA, :σH, :σC, :αWR,:αWR2, :αF, :β]

sample_size=[4803, 1405+1410,15683,3208,6009,4433,4554,8322]

 # old was [4803,1405+1410,3208,6009]

 # I got the numbers for NEWWS from line 139, and for LAGAIN from 374


# This block of code sets up some things I use to test the functions below
αHT = [0;pars.αHT]
yb = site_features.yb[1]
T = site_features.T[1]
years = (yb+1-1991):(yb-1991+T)
pos = sum(site_features.T[1:1-1])
#αH = pars.αH[(pos+1):(pos+site_features.T[i])]
pars_site = GetSitePars(pars,1)
sname = site_list[1]
Y = getfield(budget,sname)
moms = getfield(moments,sname)
π0_a = site_features.π0[2,2,:]
sum(π0_a)
π01a = site_features.π0[:,:,:]
π01 = site_features.π0[1,:,:]
year_meas = site_features.year_meas[1]
Abar = size(Y)[1]
moms_model = zeros(size(moms))
WR = site_features.work_reqs[1,2]
price = site_features.prices[1,2]
abar = min(Abar,2)
Y_I = budget[Symbol(sname,"_I")]
TLlength = site_features.TLlength[1,2]
NK = size(π01)[1]
pA1,pWork1,pF1 = GetStaticProbs(pars_site,Y[abar,:,:,:,:],price,WR,T,NK)
pA2,pWork2,pF2 = GetDynamicProbs(pars_site,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,NK,TLlength)
# end test code prep

function draw_kids(site,site_features)
    # get NK
    π0_nk = site_features.π0[site,:,:]
    prob_nk=zeros(3)
    prob_nk[1]=sum(π0_nk[1,:])
    for i in 2:3
        prob_nk[i]=sum(π0_nk[i,:])+prob_nk[i-1]
    end
    if sum(prob_nk)<0.99
        println("problem: pi")
    end
    nk_draw=rand(Uniform())
    nk=0

    A1=zeros(3)
    for i in 2:3
        A1[i]=ifelse(nk_draw>=prob_nk[i],1.0,0)
    end
    nk=convert(Int,sum(A1))+1
    if nk>3
        println("nk broken")
        println(nk)
        println(A1)
    end

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
            # im
    end

    age=convert(Int,sum(A2))+1 # sum A2 vector to get age INDEX

    return nk, age
end
draw_kids(1,site_features)
π0_nk =π01[1,:,:]
π01a

#
function decisions_static(t,nk, age,pA,pWork,pF)

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

    return A,W,F
end
decisions_static(2,3,2,pA1,pWork1,pF1)

function decisions_dynamic(t,nk,age,tl,pA,pWork, pF)
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

    return A,W,F
end
decisions_dynamic(1,2,3,1,pA2,pWork2,pF2)



A_mat1=SharedArray{Float64}(10,10)
1+1

function Moments_Simulated_Static(nsims, pA,pWork,pF, sample_size_site,T, site,site_features)

    A_mat=zeros(T,nsims, sample_size_site)
    A_sims=zeros(T,nsims)
    A=zeros(T)


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
                                                    A_mat[t,n,i],W_mat[t,n,i],F0=decisions_static(t,nk, age,pA,pWork,pF)
                                                        if (a>=1) & (a<=10)
                                                            F_mat[t,n,i]=F0
                                                            W_mat_2[t,n,i]=W_mat[t,n,i]
                                                        end
                                                    a+=1
                                                end
        end
    end
    # aggregate the moments
    @inbounds @simd for t in 1:T
        @inbounds @simd for n in 1:nsims
            A_sims[t,n]=mean(A_mat[t,n,:])
            W_sims[t,n]=mean(W_mat[t,n,:])
        end
            A[t]=mean(A_sims[t,:])
            W[t]=mean(W_sims[t,:])
    end
    F=sum(F_mat)/sum(W_mat_2)

    return A,W,F
end
@time A,W,F=Moments_Simulated_Static(2,pA1,pWork1,pF1,  sample_size[1],4,1,site_features)
@time Moments_Simulated_Static(2,pA1,pWork1,pF1, sample_size[1],4,1,site_features)
mom_d=[A;W;F]

function Moments_Simulated_Dynamic(nsims, pA,pWork,pF, sample_size_site,T, site, tlsite,site_features)
    A_mat=zeros(T,nsims, sample_size_site)
    A_sims=zeros(T,nsims)
    A=zeros(T)

    W_mat=zeros(T,nsims, sample_size_site)
    W_mat_2=zeros(T,nsims, sample_size_site)
    W_sims=zeros(T,nsims)
    W=zeros(T)

    F_mat=zeros(T,nsims, sample_size_site)

    # child care only for kids in ages 0:9
    @inbounds @simd for n in 1:nsims
         @inbounds @simd for i in 1:sample_size_site
                            tl=1
                            nk, age=draw_kids(site,site_features)
                            a=age
                                @inbounds for t in 1:T
                                        A_mat[t,n,i],W_mat[t,n,i],F0=decisions_dynamic(t,nk,age,tl,pA,pWork, pF)
                                            if (a>=1) & (a<=10)
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
                        end
                            A[t]=mean(A_sims[t,:])
                            W[t]=mean(W_sims[t,:])
                    end

    F=sum(F_mat)/sum(W_mat_2)

    return A,W,F

end
@time A,W,F=Moments_Simulated_Dynamic(1,pA2,pF2,pWork2,  sample_size[1],3,1,2,site_features)
@time Moments_Simulated_Dynamic(2,pA2,pF2,pWork2,  sample_size[1],4,1,2,site_features)
mom_d=[A;W;F]
site_features.T
site_features.TLlength
1
function Get_Simulated_MomentsAll(pars1,site_list,budget,moments,wghts,site_features,nsims,sample_size; seed=1234, wrapped=false)
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
        pA1,pWork1,pF1 = GetStaticProbs(pars_site,Y[1,:,:,:,:],price,0,T,NK) # last 3 are WR, T, NK


        A,W,F=Moments_Simulated_Static(nsims, pA1,pWork1,pF1,  ssize2,T, i,site_features)
        moms_control = [A;W;F]
        Qn += sum(wght[:,1].*(moms[:,1] .- moms_control).^2)

        if wrapped==true
            append!(moms_model, [(Part=A,LFP=W,Care=F)])
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
                A,W,F=Moments_Simulated_Dynamic(nsims,pA2,pWork2,pF2,  ssize2,T,i,TLlength,site_features)#(nsims, pA,pF,pWork, π1, sample_size_site,T, site, tlsite,site_features)
                moms_treatment = [A;W;F]

                if wrapped==true
                    append!(moms_model, [(Part=A,LFP=W,Care=F)])
                else
                    moms_model[:,a]=moms_treatment
                end
            else
                pA3,pWork3,pF3 = GetStaticProbs(pars_site,Y[abar,:,:,:,:],price,WR,T,NK)
                A,W,F=Moments_Simulated_Static(nsims,pA3,pWork3,pF3, ssize2,T, i,site_features)
                moms_treatment = [A;W;F]
                if wrapped==true
                    append!(moms_model, [(Part=A,LFP=W,Care=F)])
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

# check with baseline params
pars = parameters()
q,GS1a=Get_Simulated_MomentsAll(pars,site_list,budget,moments,wghts,site_features,1,sample_size; wrapped=true)
CriterionP(pars,site_list,budget,moments,wghts,site_features)
G1a=GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)


# check with optimized params


opt,x0 = GetOptimization(pars,vlist,site_list,budget,moments,wghts,site_features)
res2 = optimize(opt,x0)
pars2 = UpdatePars(res2[2],pars,vlist)
q1, GS1=Get_Simulated_MomentsAll(pars2,site_list,budget,moments,wghts,site_features,10,sample_size; wrapped=true)
CriterionP(pars2,site_list,budget,moments,wghts,site_features)
G1=GetMomentsAll(pars2,site_list,budget,moments,wghts,site_features)



















G1=GetMomentsAll(pars2,site_list,budget,moments,wghts,site_features)
G1.CTJF

function InspectModelFit2(model_moments,data_moments,site_features,site_list;Part=true)
    colors = ["blue","green","red"]
    for i=1:8
        T = site_features.T[i]
        sname = site_list[i]
        moms0 = getfield(data_moments,sname)
        moms1 = getfield(model_moments,sname)
        for a=1:site_features.n_arms[i]
            #for v in [:Part,:LFP]
            v=:LFP
            if Part==false
                v=:LFP
            end
            v=ifelse(Part==false,:Part,:LFP)
                figure(String(v))
                subplot(2,4,i)
                PyPlot.title(String(sname))
                PyPlot.plot(moms0[a][v],color=colors[a])
                PyPlot.plot(moms1[a][v],color=colors[a],linestyle="--")
            #end
        end
    end
end
G1=GetMomentsAll(pars2,site_list,budget,moments,wghts,site_features)
G1a=GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)


# compare participation--old pars
clf()
A=InspectModelFit2(G1a,data_moments,site_features,site_list)
gcf()
clf()
A=InspectModelFit2(GS1a,data_moments,site_features,site_list)
gcf()

#new pars
clf()
A=InspectModelFit2(G1,data_moments,site_features,site_list)
gcf()
clf()
A=InspectModelFit2(GS1,data_moments,site_features,site_list)
gcf()

# compare lfp
clf()
A=InspectModelFit2(G1a,data_moments,site_features,site_list, Part=false)
gcf()
clf()
A=InspectModelFit2(GS1a,data_moments,site_features,site_list, Part=false)
gcf()

#new pars
clf()
A=InspectModelFit2(G1,data_moments,site_features,site_list, Part=false)
gcf()
clf()
A=InspectModelFit2(GS1,data_moments,site_features,site_list, Part=false)
gcf()


1





GS1.FTP
1















1


for i=1:length(site_list)
        for a = 1:site_features.n_arms[i]
            Arm="Control"
            if a>1
                Arm="Treatment"
            end
        end
end









function GetMomentsAll2(pars,site_list,budget,moments,wghts,site_features)
    moms_collect = []
    αHT = [0;pars.αHT]
    for i=1:length(site_list)
        yb = site_features.yb[i]
        T = site_features.T[i]
        years = (yb+1-1991):(yb-1991+T)
        pos = sum(site_features.T[1:i-1])
        #αH = pars.αH[(pos+1):(pos+site_features.T[i])]
        pars_site = GetSitePars(pars,i)
        sname = site_list[i]
        Y = getfield(budget,sname)
        moms = getfield(moments,sname)
        π0 = site_features.π0[i,:,:]
        year_meas = site_features.year_meas[i]
        Abar = size(Y)[1]
        moms_model = zeros(size(moms))
        for a = 1:site_features.n_arms[i]
            WR = site_features.work_reqs[i,a]
            price = site_features.prices[i,a]
            abar = min(Abar,a)
            if site_features.time_limits[i,a]==1
                Y_I = budget[Symbol(sname,"_I")]
                TLlength = site_features.TLlength[i,a]
                moms_model[:,a] = GetMomentsTimeLims(pars_site,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,π0,TLlength,year_meas)
            else
                moms_model[:,a] = GetMoments(pars_site,Y[abar,:,:,:,:],price,WR,T,π0,year_meas)
            end
        end
        append!(moms_collect,[moms_model])
    end
    return (;zip(site_list,moms_collect)...)
end
moms_model[:,a]=
1

















T = site_features.T[1]
sname = site_list[1]
moms0 = getfield(data_moments,sname)
moms1 = getfield(GS1,sname)
moms1 = getfield(G1,sname)


for a=1:site_features.n_arms[1]
    for v in [:Part,:LFP,:Inc]
        figure(String(v))
        subplot(2,4,1)
        title(String(sname))
        plot(moms0[a][v],color=colors[a])
        plot(moms1[a][v],color=colors[a],linestyle="--")
    end
end





vlist0 = [:αc, :αH, :αA, :αF]
opt,x0 = GetOptimization(pars,vlist0,site_list,budget,moments,wghts,site_features)
res = optimize(opt,x0)
pars = UpdatePars(res[2],pars,vlist0)
Get_Simulated_MomentsAll(pars,site_list,budget,moments,wghts,site_features,10,sample_size)
CriterionP(pars,site_list,budget,moments,wghts,site_features)

# with further optimization the gap in the criterion functions grows
opt,x0 = GetOptimization(pars,vlist,site_list,budget,moments,wghts,site_features)
res2 = optimize(opt,x0)
pars2 = UpdatePars(res2[2],pars,vlist)
Get_Simulated_MomentsAll(pars2,site_list,budget,moments,wghts,site_features,10,sample_size)
CriterionP(pars2,site_list,budget,moments,wghts,site_features)



# same experiment as Jo
opt,x0 = GetOptimization(pars,[:αH],site_list,budget,moments,wghts,site_features)
res = optimize(opt,x0)
pars = UpdatePars(res[2],pars,[:αH])
@time Get_Simulated_MomentsAll(pars,site_list,budget,moments,wghts,site_features,10,sample_size)
CriterionP(pars,site_list,budget,moments,wghts,site_features)
momsnew=GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)


vlist1 = [:αc, :gN, :gF, :αA, :σH, :σC, :αWR,:αWR2, :αF,:β]
opt,x0 = GetOptimization(pars,vlist1,site_list,budget,moments,wghts,site_features)
res = optimize(opt,x0)
pars = UpdatePars(res[2],pars,vlist1)
@time Get_Simulated_MomentsAll(pars,site_list,budget,moments,wghts,site_features,10,sample_size)
CriterionP(pars,site_list,budget,moments,wghts,site_features)
momsnew=GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)


vlist2 = [:αWR,:αWR2]
opt,x0 = GetOptimization(pars,vlist2,site_list,budget,moments,wghts,site_features)
res = optimize(opt,x0)
pars = UpdatePars(res[2],pars,vlist2)
@time Get_Simulated_MomentsAll(pars,site_list,budget,moments,wghts,site_features,100,sample_size)
@time Get_Simulated_MomentsAll(pars,site_list,budget,moments,wghts,site_features,10,sample_size)
CriterionP(pars,site_list,budget,moments,wghts,site_features)
momsnew=GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)



pars = parameters()
q1, GS1=Get_Simulated_MomentsAll(pars,site_list,budget,moments,wghts,site_features,10,sample_size)
CriterionP(pars,site_list,budget,moments,wghts,site_features)
moments
GS2=GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)

GS1.CTJF










moments.CTJF










GS2.CTJF






function distance_new(x, pars)
    pars2 = UpdatePars(x,pars,vlist)
    q1, GS1=Get_Simulated_MomentsAll(pars2,site_list,budget,moments,wghts,site_features,10,sample_size)
    return q1
end


opt,x0 = GetOptimization(pars,vlist,site_list,budget,moments,wghts,site_features)
res2 = optimize(opt,x0)
pars2 = UpdatePars(res2[2],pars,vlist)
CriterionP(pars2,site_list,budget,moments,wghts,site_features)


Get_Simulated_MomentsAll(pars2,site_list,budget,moments,wghts,site_features,10,sample_size)

distance_new(res2[2], pars2)
d3(x)=distance_new(x, pars2)
d3(res2[2])

Optim.optimize(d3, res2[2])
