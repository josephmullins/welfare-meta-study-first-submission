# routines for estimating the baseline model
using ForwardDiff
using NLopt

# we want to go back to what we did before, but with adaptability?
# one idea: do a copy operation each time?

function parameters()
    np = (αc = 1, gN = 1, gF =1, αH = 35, αA = 8, σH = 8, σC = 1, wq = 1, αWR = 1, αWR2 = 1, αF = 1, αHT = 8)
    lb = (αc = 0., gN = -Inf, gF =-Inf, αH = -Inf*ones(35), αA = -Inf*ones(8), σH = zeros(8), σC = 0., wq = 0., αWR = -Inf, αWR2 = -Inf, αF = -Inf, αHT = -Inf*ones(8))
    ub = (αc = Inf, gN = Inf, gF =Inf, αH = Inf*ones(35), αA = Inf*ones(8), σH = Inf*ones(8), σC = Inf, wq = Inf, αWR = Inf, αWR2 = Inf, αF = Inf, αHT = Inf*ones(8))
    αc = 1.
    gN = 0.
    gF = 0.
    αH = zeros(35) #zeros(35)
    αA = zeros(8)
    β = 0.9
    σH = ones(8)
    σC = 1.
    wq = 2.
    αWR = 0.
    αWR2 = 0.
    αF = 0.
    αHT = zeros(8)
    Γ = zeros(18)
    return (np=np,lb=lb,ub=ub,αc=αc,gN=gN,gF=gF,αH=αH,αA=αA,β=β,σH=σH,σC=σC,wq=wq,αWR=αWR,αWR2=αWR2,αF=αF,αHT=αHT,Γ=Γ)
end

# Create a copy of the tuple, with the variables in vars updated
function UpdatePars(x,pars,vars::Array{Symbol,1})
    allvars = keys(pars)
    for v in keys(pars)
        eval(:($v=$(pars[v])))
    end
	pos = 1
	for s in vars
		np = pars.np[s]
		if typeof(pars.lb[s])==Float64
			#setfield!(Pars,s,x[pos])
            eval(:($s=$(x[pos])))
		else
			xrange = pos:(pos+np-1)
			#setfield!(Pars,s,x[xrange])
            eval(:($s=$(x[xrange])))
		end
		pos += np
	end
    vals = (eval(s) for s in allvars)
    return (;zip(allvars,vals)...)
end

function GetOptimization(pars,vars,site_list,budget,moments,wghts,site_features)
    np = 0;
	for s in vars
		np += pars.np[s]
	end
	x0 = zeros(Real,np)
	lb = zeros(np)
	ub = zeros(np)
    opt = Opt(:LD_LBFGS,np)
    pos = 1
	for s in vars
		np = pars.np[s]
		if typeof(pars.lb[s])==Float64
			x0[pos] = getfield(pars,s)
			lb[pos] = pars.lb[s]
			ub[pos] = pars.ub[s]
		else
			println(s)
			xrange = pos:(pos+np-1)
			x0[xrange] = getfield(pars,s)
			lb[xrange] = pars.lb[s]
			ub[xrange] = pars.ub[s]
		end
		pos += np
	end
	lower_bounds!(opt,lb)
	upper_bounds!(opt,ub)
	min_objective!(opt,(x,g)->Criterion(x,g,pars,vars,site_list,budget,moments,wghts,site_features))
    return opt,x0
end


function GetParameters(x,num_sites,Γ)
    # universal: αc,Γ,gN,gF,wq,σC,σH,αWR
    # site-specific: αH,αA
    αc,gN,gF,wq,σC,σH,αWR,αF = x[1:8]
    αH = x[9:(8+num_sites)]
    αA = x[(9+num_sites):(8+2*num_sites)]
    ndummies = 1999 - 1991 #<- 9 distinct years, 8 dummies
    αHT = [0;x[(9+2*num_sites):(8+ndummies+2*num_sites)]]
    return (αc=αc,gN = ones(2)*gN, gF = ones(2)*gF,wq = wq,σC = σC,σH = σH,αWR = αWR,αF = αF, αA = αA,αH=αH,Γ=Γ,αHT=αHT)
end

# site features contains: T,n_arms,work_requirement,π0,prices,year_meas
function Criterion(x,pars,vars,site_list,budget,moments,wghts,site_features)
    #num_sites = length(site_list)
    #pars = GetParameters(x,num_sites,Γ)
    pars = UpdatePars(x,pars,vars)
    Qn = CriterionP(pars,site_list,budget,moments,wghts,site_features)
    return Qn
end


# an extra two parameters, let's see how we do
function CriterionP(pars,site_list,budget,moments,wghts,site_features)
    Qn = 0
    αHT = [0;pars.αHT]
    for i=1:length(site_list)
        yb = site_features.yb[i]
        T = site_features.T[i]
        years = (yb+1-1991):(yb-1991+T)
        pos = sum(site_features.T[1:i-1])
        αH = pars.αH[(pos+1):(pos+site_features.T[i])]
        pars_site = (αc = pars.αc,gN = ones(2)*pars.gN,gF = ones(2)*pars.gF,wq = pars.wq,σC = pars.σC,σH = pars.σH[i],αWR = pars.αWR,αWR2 = pars.αWR2,αF = pars.αF,αH = αH,αA = pars.αA[i],Γ=pars.Γ,β=0.)
        sname = site_list[i]
        Y = getfield(budget,sname)
        moms = getfield(moments,sname)
        π0 = site_features.π0[i,:,:]
        wght = getfield(wghts,sname)
        year_meas = site_features.year_meas[i]
        Abar = size(Y)[1]
        for a = 1:site_features.n_arms[i]
            WR = site_features.work_reqs[i,a]
            price = site_features.prices[i,a]
            abar = min(Abar,a)
            moms_model = GetMoments(pars_site,Y[abar,:,:,:,:],price,WR,T,π0,year_meas)
            Qn += sum(wght[:,a].*(moms[:,a] .- moms_model).^2)
        end
    end
    return Qn
end

function GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)
    moms_collect = []
    αHT = [0;pars.αHT]
    for i=1:length(site_list)
        yb = site_features.yb[i]
        T = site_features.T[i]
        years = (yb+1-1991):(yb-1991+T)
        pos = sum(site_features.T[1:i-1])
        αH = pars.αH[(pos+1):(pos+site_features.T[i])]
        pars_site = (αc = pars.αc,gN = pars.gN,gF = pars.gF,wq = pars.wq,σC = pars.σC,σH = pars.σH[i],αWR = pars.αWR,αWR2 = pars.αWR2,αF = pars.αF,αH = αH,αA = pars.αA[i],Γ=pars.Γ,β=0.)
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
            moms_model[:,a] = GetMoments(pars_site,Y[abar,:,:,:,:],price,WR,T,π0,year_meas)
        end
        append!(moms_collect,[moms_model])
    end
    return (;zip(site_list,moms_collect)...)
end


function Criterion(x,g,pars,vars,site_list,budget,moments,wghts,site_features)
    f0 = Criterion(x,pars,vars,site_list,budget,moments,wghts,site_features)
    dF = ForwardDiff.gradient(x->Criterion(x,pars,vars,site_list,budget,moments,wghts,site_features),x)
    g[:] = dF
    return f0
end

function GetMoments(pars,Y,price,WR,T,π0,year_meas)
    NK = size(π0)[1]
    pA,pWork,pF = GetStaticProbs(pars,Y,price,WR,T,NK)
    EA,EH = GetAggMoments(pA,pWork,π0)
    Care = GetMeanCare(year_meas,0,9,pA,pWork,pF,π0)
    return [EA; EH; Care]
end

function InspectModelFit(model_moments,data_moments,site_features,site_list)
    colors = ["red","green","blue"]
    for i=1:8
        T = site_features.T[i]
        moms0 = getfield(model_moments,site_list[i])
        moms1 = getfield(data_moments,site_list[i])
        for a=1:site_features.n_arms[i]
            figure("AFDC")
            subplot(2,4,i)
            title(String(site_list[i]))
            plot(moms0[1:T,a],color=colors[a])
            plot(moms1[1:T,a],color=colors[a],linestyle="--")
            figure("LFP")
            subplot(2,4,i)
            title(String(site_list[i]))
            plot(moms0[T+1:2*T,a],color=colors[a])
            plot(moms1[T+1:2*T,a],color=colors[a],linestyle="--")
        end
    end
end
