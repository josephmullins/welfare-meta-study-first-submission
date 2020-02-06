using Distributions
using ForwardDiff

function ProdPars(nmeas)
    pos = (δI = 1:2, δθ = 3, gN = 4:5,gF=6:7,λ = 8:(8+nmeas-2))
    δI = 0.1 *ones(2)
    δθ = 1.
    Qδ(x) = exp.(log.(x) .+ 0.3*quantile.(Normal(),rand(2)))
    μδ = log.([0.5,0.5])
    Fδ(x) = -0.5*sum(( ((log.(x) .- μδ)/0.5).^2)) - 2*log(0.5) #<- log prior
    Qθ(x) = rand() #exp.(log.(x) .+ 0.1*quantile.(Normal(),rand(2)))
    Fθ(x) = 0. #<- constant density, uniform prior
    gN = zeros(2) #log.([1.5,1.5])
    gF = zeros(2) #log.([0.7,0.7])
    QN(x) = x .+ 0.25*quantile.(Normal(),rand(2))
    QF(x) = x .+ 0.25*quantile.(Normal(),rand(2))
    # prior at zero, lots of uncertainty
    FN(x) = -0.5*sum(( ((x)/1.).^2)) - 2*log(1.)
    FF(x) = -0.5*sum(( ((x)/1.).^2)) - 2*log(1.)
    #λ = zeros(nmeas-1)
    λ = ones(nmeas-1)
    μλ = ones(nmeas-1)
    σλ = 0.4
    Fλ(x) = -0.5*sum(( ((log.(x) .- μλ)/σλ).^2)) - 2*log(σλ)
    Qλ(x) = exp.(log.(x) .+ 0.25*quantile.(Normal(),rand(nmeas-1)))
    proposal = (δI = Qδ,δθ = Qθ,gN = QN,gF = QF,λ=Qλ)
    logprior = (δI = Fδ,δθ = Fθ,gN = FN,gF = FF,λ=Fλ)
    x0 = [δI;δθ;gN;gF;λ]
    return x0,(pos=pos,Qp=proposal,P0=logprior)
end

#vlist = [:αc, :αH, :αA, :σH, :σC, :αWR,:αWR2, :αF,:wq,:gN,:β]

function ModPars(Γ = zeros(18),gF = 0.)
    pos = (αc = 1, αH = 2:9, αA = 10:17, σH = 18, σC = 19, αWR = 20:27, αWR2 = 28:35, αF = 36:43, wq = 44, gN = 45, β = 46)
    # αc
    Qc(x) = exp.(log.(x) .+ 0.005*quantile(Normal(),rand()))
    αc = 1.
    # gN
    gN = 0.
    QN(x) = x + 0.2*quantile(Normal(),rand())
    # αH
    αH = 1*ones(8)
    QH(x) = x .+ 0.05*quantile.(Normal(),rand(8))
    # αA
    QA(x) = x .+ 0.01*quantile.(Normal(),rand(8))
    αA = 1*ones(8)
    # αF
    QF(x) = x .+ 0.15*quantile.(Normal(),rand(8))
    αF = zeros(8)
    # αWR
    QWR(x) = x .+ 0.05*quantile.(Normal(),rand(8))
    αWR = zeros(8)
    # αWR2
    QWR2(x) = x .+ 0.05*quantile.(Normal(),rand(8))
    αWR2 = zeros(8)
    # σH
    QσH(x) = exp(log(x) + 0.05*quantile(Normal(),rand()))
    σH = 1.
    # σC
    QσC(x) = exp(log(x) + 0.05*quantile(Normal(),rand()))
    σC = 1.
    # β
    Qβ(x) = 0.4 + 0.6*rand()
    β = 0.8
    # wq
    Qwq(x) = exp(log(x) + 0.01*quantile(Normal(),rand())) #<- 0.01
    wq = 1.

    x0 = [αc;αH;αA;σH;σC;αWR;αWR2;αF;wq;gN;β]
    proposal = (αc = Qc,αH = QH, αA = QA, αF = QF, αWR = QWR, αWR2 = QWR2, gN = QN, σH = QσH, σC = QσC, β = Qβ, wq = Qwq)
    #logprior = (αc = Fc,αH = FH, αA = FA, αF = FF, αWR = FWR, αWR2 = FWR2, gN = FN, σH = FσH, σC = FσC, β = Fβ, wq = Fwq)

    return x0,(pos = pos, Qp = proposal,gF = gF,Γ=Γ)
end

function HyperPars()
    Q(x) = x + 2.5*quantile(Normal(),rand())
    F(x) = -0.5*sum((x/50).^2)
    Qσ(x) = exp(log(x) + 0.75*quantile(Normal(),rand()))
    Fσ(x) = -0.5*((log(x)-0)/50)^2
    proposal = (αH=Q,σαH = Qσ, αA=Q,σαA = Qσ, αF=Q, σαF = Qσ, αWR=Q, σαWR = Qσ, αWR2=Q, σαWR2 = Qσ)
    logprior = (αH=F,σαH = Fσ, αA=F,σαA = Fσ, αF=F, σαF = Fσ, αWR=F, σαWR = Fσ, αWR2=F, σαWR2 = Fσ)
    parmap = (αH = :αH, σαH = :αH, αA = :αA,σαA = :αA,αF = :αF,σαF = :αF,αWR = :αWR,σαWR = :αWR,αWR2 = :αWR2,σαWR2 = :αWR2)
    pos = (αH = 1, σαH = 2, αA = 3,σαA = 4,αF = 5,σαF = 6,αWR = 7,σαWR = 8,αWR2 = 9,σαWR2 = 10)
    x0 = [1.,1.,1.,1.,0.,1.,0.,1.,0.,1.] #<- too many parameters?
    return x0,(pos = pos, Qp = proposal, P0 = logprior, parmap = parmap)
end

function GetProdParameters(x,pars)
    return (δI = x[pars.pos.δI],δθ = x[pars.pos.δθ],gN = x[pars.pos.gN],gF = x[pars.pos.gF],λ = x[pars.pos.λ])
end

function GetModParameters(x,pars)
    return (αc=x[pars.pos.αc],gN=x[pars.pos.gN],αH=x[pars.pos.αH],αA=x[pars.pos.αA],β=x[pars.pos.β],σH=x[pars.pos.σH],σC=x[pars.pos.σC],wq=x[pars.pos.wq],αWR=x[pars.pos.αWR],αWR2=x[pars.pos.αWR2],αF=x[pars.pos.αF],gF = pars.gF,Γ=pars.Γ)
end

function GetSitePars(pars,i)
    return (αc = pars.αc,gN = pars.gN,gF = pars.gF,wq = pars.wq,σC = pars.σC,σH = pars.σH,αWR = pars.αWR[i],αWR2 = pars.αWR2[i],αF = pars.αF[i],αH = pars.αH[i],αA = pars.αA[i],Γ=pars.Γ,β=pars.β)
end


function LogLikeMod(x,pars,budget,moments,site_features)
    pmod = GetModParameters(x,pars)
    Qn = 0
    for i=1:length(site_features.site_list)
        yb = site_features.yb[i]
        T = site_features.T[i]
        years = (yb+1-1991):(yb-1991+T)
        pars_site = GetSitePars(pmod,i)
        sname = site_features.site_list[i]
        #println(sname)
        Y = getfield(budget,sname)
        moms = getfield(moments,sname)
        π0 = site_features.π0[i,:,:]
        year_meas = site_features.year_meas[i]
        Abar = size(Y)[1]
        price = site_features.prices[i,1]
        #moms_control = GetMoments(pars_site,Y[1,:,:,:,:],price,0,T,π0,year_meas)

        for a = 1:site_features.n_arms[i]
            WR = site_features.work_reqs[i,a]
            price = site_features.prices[i,a]
            abar = min(Abar,a)
            if site_features.time_limits[i,a]==1
                Y_I = budget[Symbol(sname,"_I")]
                TLlength = site_features.TLlength[i,a]
                moms_model = GetMomentsTimeLims(pars_site,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,π0,TLlength,year_meas)
            else
                moms_model = GetMoments(pars_site,Y[abar,:,:,:,:],price,WR,T,π0,year_meas)
            end
            #TE = moms[:,a] .- moms[:,1]
            Qn += -0.5*sum(((moms.moms[:,a] .- moms_model)./moms.se[:,a]).^2) - sum(log.(moms.se[:,a]))
        end
    end
    return Qn
end

function LogLikeMod(x,g,pars,budget,moments,site_features)
    f0 = LogLikeMod(x,pars,budget,moments,site_features) + LogPrior([],x,mpars,hpars,:wq)
    dL = ForwardDiff.gradient(x->LogLikeMod(x,pars,budget,moments,site_features),x)
    g[:] = dL
    return f0
end

function LogLikeProd(x,pars,ChoiceProbs,site_list,budget,moments,site_features)
    pars_prod = GetProdParameters(x,pars)
    wq = 1. # <- until we figure out a better way to estimate, leave at this
    Qn = 0
    NK=3 # might want to update this later
    λ = [1;pars_prod.λ]
    for i=1:length(site_list)
        sname = site_list[i]
        #println(sname)
        moms = getfield(moments,sname) #<- structure: age0,age1,1
        nmom = size(moms.TE)[1]
        π0 = site_features.π0[i,:,:]
        #wght = getfield(wghts,sname)
        year_meas = site_features.year_meas[i]
        Y = getfield(budget,sname)
        P = getfield(ChoiceProbs,sname)
        TH = zeros(Real,site_features.n_arms[i],NK,17)
        price = site_features.prices[i,1]
        TH[1,:,:] =  GetChildOutcomesStatic(year_meas,P.pA[1],P.pWork[1],P.pF[1],Y[1,:,:,:,:],pars_prod,price,wq)

        Abar = size(Y)[1]
        #moms_control = GetMoments(pars_site,Y[1,:,:,:,:],price,0,T,π0,year_meas)
        #Qn += sum(wght[:,1].*(moms[:,1] .- moms_control).^2)
        for a = 2:site_features.n_arms[i]
            WR = site_features.work_reqs[i,a]
            price = site_features.prices[i,a]
            abar = min(Abar,a)
            if site_features.time_limits[i,a]==1
                Y_I = budget[Symbol(sname,"_I")]
                TLlength = site_features.TLlength[i,a]
                TH[a,:,:] = GetChildOutcomesDynamic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],pars_prod,price,wq)
            else
                TH[a,:,:] = GetChildOutcomesStatic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[abar,:,:,:,:],pars_prod,price,wq)
            end
        end
        TE = TreatmenEffects(TH,π0,moms.a0,moms.a1,moms.arm)
        for k=1:size(moms.TE)[2]
            for m=1:size(moms.TE)[1]
                if moms.wght[m,k]>0
                    Qn += -0.5*((λ[k]*TE[m] .- moms.TE[m,k])/moms.SE[m,k])^2 - log(moms.SE[m,k])
                end
            end
        end
    end
    return Qn
end

# variable cases:
#vlist = [:αc, :αH, :αA, :σH, :σC, :αWR,:αWR2, :αF, :β,:wq,:gN]
function LogPrior(xh,xm,mpars,hpars,var)
    if var==:αH
        αH = xm[mpars.pos.αH]
        μα = xh[hpars.pos.αH]
        σα = xh[hpars.pos.σαH]
        return -0.5*sum(((αH .- μα)/σα).^2) - length(αH)*log(σα)
    elseif var==:αA
        α = xm[mpars.pos.αA]
        μα = xh[hpars.pos.αA]
        σα = xh[hpars.pos.σαA]
        return -0.5*sum(((α .- μα)/σα).^2) - length(α)*log(σα)
    elseif var==:αF
        α = xm[mpars.pos.αF]
        μα = xh[hpars.pos.αF]
        σα = xh[hpars.pos.σαF]
        return -0.5*sum(((α .- μα)/σα).^2) - length(α)*log(σα)
    elseif var==:αWR
        α = xm[mpars.pos.αWR]
        μα = xh[hpars.pos.αWR]
        σα = xh[hpars.pos.σαWR]
        return -0.5*sum(((α .- μα)/σα).^2) - length(α)*log(σα)
    elseif var==:αWR2
        α = xm[mpars.pos.αWR2]
        μα = xh[hpars.pos.αWR2]
        σα = xh[hpars.pos.σαWR2]
        return -0.5*sum(((α .- μα)/σα).^2) - length(α)*log(σα)
    elseif (var==:αc) | (var==:σH) | (var==:σC)
        x = xm[mpars.pos[var]]
        return -0.5*(x/20)^2 - log(20)
    elseif (var==:wq)
        x = xm[mpars.pos[var]]
        return -0.5*((log(x)-log(2.))/0.1)^2 - log(0.5)
        #return (x<=5)*0 + (x>5)*-Inf
    elseif (var==:β) #
        return 0.
    elseif var==:gN
        gN = xm[mpars.pos.gN]
        return -0.5*(gN/20)^2
    end
end

function LogLikeHyper(x,xmod,mpars,hpars)
    ll = 0
    for v in keys(mpars.pos)
        ll += LogPrior(x,xmod,mpars,hpars,v)
    end
    return ll
end

function LogLikeHyper(x,g,xmod,mpars,hpars)
    ll = LogLikeHyper(x,xmod,mpars,hpars)
    dL = ForwardDiff.gradient(x->LogLikeHyper(x,xmod,mpars,hpars),x)
    g[:] = dL
    return ll
end



function IterateMCMCHyper(xcurrent,xmod,mpars,hpars,vname)
    par_current = xcurrent[hpars.pos[vname]]
    par_propose = hpars.Qp[vname](par_current)
    x1 = copy(xcurrent)
    x1[hpars.pos[vname]] = par_propose
    vapply = hpars.parmap[vname]
    #pars.pos[vname]
    ll1 = LogPrior(x1,xmod,mpars,hpars,vapply)
    ll0 = LogPrior(xcurrent,xmod,mpars,hpars,vapply)
    p0 = hpars.P0[vname](par_current)
    p1 = hpars.P0[vname](par_propose)
    accept_prob = exp(ll1 + p1 - ll0 - p0)
    if rand()<accept_prob
        return x1,1
    else
        return xcurrent,0
    end
end

function IterateMCMCProd(xcurrent,vname,pars,ll0,ChoiceProbs,site_list,budget,moments,site_features)
    # current parameter values
    par_current = xcurrent[pars.pos[vname]]
    par_propose = pars.Qp[vname](par_current)
    x1 = copy(xcurrent)
    x1[pars.pos[vname]] = par_propose
    ll1 = LogLikeProd(x1,pars,ChoiceProbs,site_list,budget,moments,site_features)
    p0 = pars.P0[vname](par_current)
    p1 = pars.P0[vname](par_propose)
    accept_prob = exp(ll1 + p1 - ll0 - p0)
    if rand()<accept_prob
        return x1,ll1,1
    else
        return xcurrent,ll0,0
    end
end


function IterateMCMCMod(xcurrent,xhyper,vname,mpars,hpars,ll0,budget,moments,site_features)
    # current parameter values
    par_current = xcurrent[mpars.pos[vname]]
    par_propose = mpars.Qp[vname](par_current)
    x1 = copy(xcurrent)
    x1[mpars.pos[vname]] = par_propose
    ll1 = LogLikeMod(x1,mpars,budget,moments,site_features)
    p0 = LogPrior(xhyper,xcurrent,mpars,hpars,vname)
    p1 = LogPrior(xhyper,x1,mpars,hpars,vname)
    accept_prob = exp(ll1 + p1 - ll0 - p0)
    if rand()<accept_prob
        return x1,ll1,1
    else
        return xcurrent,ll0,0
    end
end

function GetChainMod(N,xm0,xh0,mlist,hlist,mpars,hpars,budget,moments,site_features)
    mnp = length(xm0)
    hnp = length(xh0)
    Ph = zeros(hnp,N)
    Pm = zeros(mnp,N)
    nm = length(mlist)
    nh = length(hlist)
    Amhist = zeros(nm,N)
    Ahhist = zeros(nh,N)
    Lhist = zeros(N)
    # calculate initial likelihood
    ll = LogLikeMod(xm0,mpars,budget,moments,site_features)
    Pm[:,1] = xm0
    Ph[:,1] = xh0
    for i=2:N
        # first part: lower level
        xm = Pm[:,i-1]
        xh = Ph[:,i-1]
        lp = 0
        for v=1:nm # (1:nm)
            vname = mlist[v]
            #println(vname)
            xm,ll,aresult = IterateMCMCMod(xm,xh,vname,mpars,hpars,ll,budget,moments,site_features)
            Amhist[v,i] = aresult
            # update parameters
            lp += LogPrior(xh,xm,mpars,hpars,vname)
        end
        Pm[:,i] = xm
        Lhist[i] = ll + lp
        # higher level: iterate over hyperparameters
        for v=1:nh
            vname = hlist[v]
            xh,aresult = IterateMCMCHyper(xh,xm,mpars,hpars,vname)
            Ahhist[v,i] = aresult
        end
        Ph[:,i] = xh
        if mod(i,100)==0
            println(i)
        end
    end
    return Pm,Ph,Amhist,Ahhist,Lhist
end


function GetChainProd(N,x0,vlist,pars,ChoiceProbs,site_list,budget,moments,site_features)
    np = length(x0)
    P = zeros(np,N)
    nv = length(vlist)
    Ahist = zeros(nv,N)
    Lhist = zeros(N)
    # calculate initial likelihood
    ll = LogLike(x0,pars,ChoiceProbs,site_list,budget,moments,site_features)
    P[:,1] = x0
    for i=2:N
        xcurrent = P[:,i-1]
        for v=1:nv
            vname = vlist[v]
            xcurrent,ll,aresult = IterateMCMCProd(xcurrent,vname,pars,ll,ChoiceProbs,site_list,budget,moments,site_features)
            Ahist[v,i] = aresult
            # update parameters
        end
        P[:,i] = xcurrent
        Lhist[i] = ll
        if mod(i,100)==0
            println(i)
        end
    end
    return P,Ahist,Lhist
end



function GetMoments(pars,Y,price,WR,T,π0,year_meas,wrapped=false)
    NK = size(π0)[1]
    pA,pWork,pF = GetStaticProbs(pars,Y,price,WR,T,NK)
    EA,EH,Inc = GetAggMoments(pA,pWork,π0,Y)
    Care = GetMeanCare(year_meas,0,9,pA,pWork,pF,π0)
    if wrapped
        return (Part = EA,LFP=EH,Inc = Inc,Care=Care)
    else
        return [EA; EH; Care]
    end
end

function GetMomentsTimeLims(pars,Y,Y_I,price,WR,T,π0,TLlength,year_meas,wrapped=false)
    NK = size(π0)[1]
    pA,pWork,pF = GetDynamicProbs(pars,Y,Y_I,price,WR,T,NK,TLlength)
    EA,EH,Care,Inc = GetDynamicMoments(pA,pWork,pF,π0,year_meas,0,9,Y,Y_I)
    #EA,EH = GetAggMoments(pA,pWork,π0)
    #Care = GetMeanCare(year_meas,0,9,pA,pWork,pF,π0)
    if wrapped
        return (Part = EA,LFP=EH,Inc = Inc,Care=Care)
    else
        return [EA; EH; Care]
    end
end
