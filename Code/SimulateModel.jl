using Revise
using Distributions
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

function draw_kids(π0, site,site_features)
    # get NK
    π01 = site_features.π0[site,:,:]
    prob_nk=zeros(3)
    prob_nk[1]=sum(π01[1,:])
    for i in 2:3
        prob_nk[i]=sum(π01[i,:])+prob_nk[i-1]
    end
    prob_nk
    nk_draw=rand(Uniform())
    nk=0

    A1=zeros(3)
    for i in 2:3
        A1[i]=ifelse(nk_draw>=prob_nk[i],1.0,0)
    end
    nk=convert(Int,sum(A1))+1
    if nk>3
        println(nk)
        println(A1)
    end

    # get age
    π02 = site_features.π0[site,nk,:]
    π02=π02./sum(π02)
    prob_age=zeros(17)
    prob_age[1]=π02[1]
    for i in 2:17
        prob_age[i]=π02[i]+prob_age[i-1]
    end
    prob_age
    age_draw=rand(Uniform())
    A2=zeros(17)
    for i in 1:17
        A2[i]=ifelse(age_draw>=prob_age[i],1.0,0) # A2[i]=1 if weakly older than age i
    end
    A2
    age=convert(Int,sum(A2))+1 # sum A2 vector to get age INDEX

    return nk, age
end
draw_kids(π0,1,site_features)



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





function Moments_Simulated_Static(nsims, pA,pF,pWork, π1, sample_size_site,T, site,site_features)

    A_mat=zeros(T,nsims, sample_size_site)
    A_sims=zeros(T,nsims)
    A=zeros(T)

    W_mat=zeros(T,nsims, sample_size_site)
    W_mat_2=zeros(T,nsims, sample_size_site) # for childcare calculation
    W_sims=zeros(T,nsims)
    W=zeros(T)

    F_mat=zeros(T,nsims, sample_size_site)
    F0=0
    F=0
    # note: only give benefits if kids aged 0-9, aka indices 1-10
    for n in 1:nsims
        for i in 1:sample_size_site
            nk, age=draw_kids(π1,site,site_features)
            a=age
            for t in 1:T
                A_mat[t,n,i],W_mat[t,n,i],F0=decisions_static(t,nk, age,pA,pWork,pF)
                a+=1
                if (a>=1) & (a<=10)
                    F_mat[t,n,i]=F0
                    W_mat_2[t,n,i]=W_mat[t,n,i]
                end
            end
        end
    end
    # aggregate the moments
    for t in 1:T
        for n in 1:nsims
        A_sims[t,n]=mean(A_mat[t,n,:])
        W_sims[t,n]=mean(W_mat[t,n,:])
        end
        A[t]=mean(A_sims[t,:])
        W[t]=mean(W_sims[t,:])
    end
    F=sum(F_mat)/sum(W_mat_2)

    return A,W,F
end
@time A,W,F=Moments_Simulated_Static(2,pA1,pF1,pWork1, π01, sample_size[1],4,1,site_features)
@time Moments_Simulated_Static(2,pA1,pF1,pWork1, π01, sample_size[1],4,1,site_features)


function Moments_Simulated_Dynamic(nsims, pA,pF,pWork, π1, sample_size_site,T, site, tlsite,site_features)
    A_mat=zeros(T,nsims, sample_size_site)
    A_sims=zeros(T,nsims)
    A=zeros(T)

    W_mat=zeros(T,nsims, sample_size_site)
    W_mat_2=zeros(T,nsims, sample_size_site)
    W_sims=zeros(T,nsims)
    W=zeros(T)

    F_mat=zeros(T,nsims, sample_size_site)
    F_0=0
    F=0


    # child care only for kids in ages 0:9
    for n in 1:nsims
        for i in 1:sample_size_site
            tl=1
            nk, age=draw_kids(π1,site,site_features)
            a=age
            for t in 1:T
                A_mat[t,n,i],W_mat[t,n,i],F0=decisions_dynamic(t,nk,age,tl,pA,pWork, pF)
                a+=1
                if (a>=1) & (a<=10)
                    F_mat[t,n,i]=F0
                    W_mat_2[t,n,i]=W_mat[t,n,i]
                end
                if A_mat[t,n,i]==1
                    tl=min(tl+1,tlsite)
                end
            end
        end
    end
    for t in 1:T
        for n in 1:nsims
        A_sims[t,n]=mean(A_mat[t,n,:])
        W_sims[t,n]=mean(W_mat[t,n,:])
        end
        A[t]=mean(A_sims[t,:])
        W[t]=mean(W_sims[t,:])
    end

    F=sum(F_mat)/sum(W_mat_2)

    return A,W,F

end
@time A,W,F=Moments_Simulated_Dynamic(2,pA2,pF2,pWork2, π01, sample_size[1],4,1,3,site_features)
mom_d=[A;W;F]
@time Moments_Simulated_Dynamic(2,pA2,pF2,pWork2, π01, sample_size[1],4,1,3,site_features)



function Get_Simulated_MomentsAll(pars,site_list,budget,moments,wghts,site_features,nsims,sample_size)
    moms_collect = []
    Qn = 0
    αHT = [0;pars.αHT]
    for i=1:length(site_list)
        NK=3
        yb = site_features.yb[i]
        T = site_features.T[i]
        years = (yb+1-1991):(yb-1991+T)
        pos = sum(site_features.T[1:i-1])
        #αH = pars.αH[(pos+1):(pos+site_features.T[i])]
        pars_site = GetSitePars(pars,i)
        sname = site_list[i]
        Y = getfield(budget,sname)
        moms = getfield(moments,sname)
        π1 = site_features.π0[i,:,:]
        wght = getfield(wghts,sname)
        year_meas = site_features.year_meas[i]
        Abar = size(Y)[1]
        moms_model = zeros(size(moms))
        price = site_features.prices[i,1]
        pA1,pWork1,pF1 = GetStaticProbs(pars_site,Y[1,:,:,:,:],price,0,T,NK) # last 3 are WR, T, NK
        ssize2=convert(Int,(round(sample_size[i]))) # will replace later with actual arm sample sizes
        A,W,F=Moments_Simulated_Static(nsims, pA1,pF1,pWork1, π1, ssize2,T, i,site_features)
        moms_control = [A;W;F]
        Qn += sum(wght[:,1].*(moms[:,1] .- moms_control).^2)
        moms_model[:,1]=moms_control


        for a = 2:site_features.n_arms[i]
            WR = site_features.work_reqs[i,a]
            price = site_features.prices[i,a]
            abar = min(Abar,a)
            ssize2=convert(Int,(round(sample_size[i])))
            if site_features.time_limits[i,a]==1
                Y_I = budget[Symbol(sname,"_I")]
                TLlength = site_features.TLlength[i,a]
                pA2,pWork2,pF2 = GetDynamicProbs(pars_site,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,NK,TLlength)
                A,W,F=Moments_Simulated_Dynamic(2,pA2,pF2,pWork2, π01, ssize2,T,i,TLlength,site_features)#(nsims, pA,pF,pWork, π1, sample_size_site,T, site, tlsite,site_features)
                moms_treatment = [A;W;F]
                moms_model[:,a]=moms_treatment
            else
                pA1,pWork1,pF1 = GetStaticProbs(pars_site,Y[abar,:,:,:,:],price,WR,T,NK) # last 3 are WR, T, NK
                A,W,F=Moments_Simulated_Static(nsims, pA1,pF1,pWork1, π1, ssize2,T, i,site_features)
                moms_treatment = [A;W;F]
                moms_model[:,a]=moms_treatment
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
Get_Simulated_MomentsAll(pars,site_list,budget,moments,wghts,site_features,10,sample_size)
CriterionP(pars,site_list,budget,moments,wghts,site_features)
moments
GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)


# check with optimized params
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
