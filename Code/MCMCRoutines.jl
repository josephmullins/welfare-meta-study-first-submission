
function ParsBayes(nmeas)
    pos = (δI = 1:2, δθ = 3, gN = 4:5,gF=6:7,λ = 8:(8+nmeas-2))
    δI = 0.1 *ones(2)
    δθ = 1.
    Qδ(x) = exp.(log.(x) .+ 0.3*quantile.(Normal(),rand(2)))
    μδ = log.([0.5,0.5])
    Fδ(x) = -0.5*sum(( ((log.(x) .- μδ)/0.5).^2)) - 2*log(0.5) #<- log prior
    Qθ(x) = rand() #exp.(log.(x) .+ 0.1*quantile.(Normal(),rand(2)))
    Fθ(x) = x
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

function GetProdParameters(x,pars)
    return (δI = x[pars.pos.δI],δθ = x[pars.pos.δθ],gN = x[pars.pos.gN],gF = x[pars.pos.gF],λ = x[pars.pos.λ])
end


function LogLike(x,pars,ChoiceProbs,site_list,budget,moments,site_features)
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

function IterateMCMC(xcurrent,vname,pars,ll0,ChoiceProbs,site_list,budget,moments,site_features)
    # current parameter values
    par_current = xcurrent[pars.pos[vname]]
    par_propose = pars.Qp[vname](par_current)
    x1 = copy(xcurrent)
    x1[pars.pos[vname]] = par_propose
    ll1 = LogLike(x1,pars,ChoiceProbs,site_list,budget,moments,site_features)
    p0 = pars.P0[vname](par_current)
    p1 = pars.P0[vname](par_propose)
    accept_prob = exp(ll1 + p1 - ll0 - p0)
    if rand()<accept_prob
        return x1,ll1,1
    else
        return xcurrent,ll0,0
    end
end

function GetChain(N,x0,vlist,pars,ChoiceProbs,site_list,budget,moments,site_features)
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
            xcurrent,ll,aresult = IterateMCMC(xcurrent,vname,pars,ll,ChoiceProbs,site_list,budget,moments,site_features)
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
