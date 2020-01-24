using Revise
using Distributions

includet("BaselineModel.jl")
includet("EstimationBaseline.jl")
includet("SetupBaseline.jl")

Γ = exp.(0 .+ -0.1*(1:18))
pars = parameters()
vlist = [:αc, :αH, :αA, :σH, :σC, :αWR,:αWR2, :αF, :β]



sample_size=[4803, 1405+1410,15683,4433,4554,8322]

 # old was [4803,1405+1410,3208,6009]

 # I got the numbers for NEWWS from line 139, and for LAGAIN from 374

# Below: scratch notes to help me get under the hood of Jo's code

sname=site_list[1]
Y = getfield(budget,sname)

moms=GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)
moms.CTJF


site_features.yb

site_list

moms_collect = []
αHT = [0;pars.αHT]
π0
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
sum(π01)
sum(π01[3,:])
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


arm_selector=rand(Uniform())

NK = size(π0)[1]
pA,pWork,pF = GetDynamicProbs(pars,Y,Y_I,price,WR,T,NK,TLlength)

# this is taken from the criterion function
Qn = 0
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
    π1 = site_features.π0[i,:,:]
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
            moms_model[:,a] = GetMomentsTimeLims(pars_site,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,π1,TLlength,year_meas)
        else
            moms_model[:,a] = GetMoments(pars_site,Y[abar,:,:,:,:],price,WR,T,π1,year_meas)
        end
    end
    append!(moms_collect,[moms_model])
end
return q1=(;zip(site_list,moms_collect)...)
q1

moms_collect
#simulation pseudocode


function draw_kids(π0, site)
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
    for i in 1:3
        A1[i]=ifelse(nk_draw>=prob_nk[i],1.0,0)
    end
    nk=convert(Int,sum(A1))+1

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
draw_kids(π0,1)


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


# this is for each arm





function Moments_Simulated_Static(nsims, pA,pF,pWork, π1, sample_size_site,T, site)

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
            nk, age=draw_kids(π1,site)
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

@time Moments_Simulated_Static(2,pA1,pF1,pWork1, π01, sample_size[1],4,1)
@time Moments_Simulated_Static(2,pA1,pF1,pWork1, π01, sample_size[1],4,1)
1+1

function Moments_Simulated_Dynamic(nsims, pA,pF,pWork, π1, sample_size_site,T, site, tlsite)
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
            nk, age=draw_kids(π1,site)
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


@time Moments_Simulated_Dynamic(2,pA2,pF2,pWork2, π01, sample_size[1],4,1,3)
@time Moments_Simulated_Dynamic(2,pA2,pF2,pWork2, π01, sample_size[1],4,1,3)

convert(Int,(round(1.5)))





for k in 1:length(nsims)
    #arm_selector=rand(Uniform())
    #arm=0 # treatment or control
    #if arm_selector>=0.5
    #    arm=1
    #end

    # draw random to select number of kids
    # for everyone in sample size

            for i in 1:length(sample_size[1]) # for everyone in sample size
            nk,age=draw_kids(π0,1)

        for t in 1:T
            # draw uniform random to decide if work or not
            # draw uniform random to decide if participate or not
                # update time limit if participate
                # draw uniform random  to decide if formal or not
            # It seems the age index is fixed and only the time index changes?
            end
        end

end

test=randn(2,2)
test2=zeros(2)
