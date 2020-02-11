# routines for estimating the baseline model
using ForwardDiff
using NLopt
using PyPlot

# we want to go back to what we did before, but with adaptability?
# one idea: do a copy operation each time?

function parameters(Γ = zeros(18),gF = 0.)
    np = (αc = 8, gN = 1, gF =1, αH = 8, αA = 8, σH = 1, σC = 1, wq = 1, αWR = 1, αWR2 = 1, αF = 8, αHT = 8, β = 1)
    lb = (αc = 0. *ones(8), gN = -Inf, gF =-Inf, αH = -Inf*ones(8), αA = -Inf*ones(8), σH = 0., σC = 0., wq = 0. , αWR = -Inf, αWR2 = -Inf, αF = -Inf*ones(8), αHT = -Inf*ones(8), β = 0.)
    ub = (αc = Inf*ones(8), gN = Inf, gF =Inf, αH = Inf*ones(8), αA = Inf*ones(8), σH = Inf, σC = Inf, wq = 5., αWR = Inf, αWR2 = Inf, αF = Inf*ones(8), αHT = Inf*ones(8), β = 1.)
    αc = 1. *ones(8)
    gN = 0.
    αH = 0.1 .+ zeros(8) #zeros(35) #zeros(35)
    αA = zeros(8)
    β = 0.9
    σH = 1. #ones(8)
    σC = 1.
    wq = 1. #* ones(8)
    αWR = 0. #ones(8)
    αWR2 = 0. #ones(8)
    αF = 0. * ones(8)
    αHT = zeros(8)

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

function GetSitePars(pars,i)
    return (αc = pars.αc[i],gN = pars.gN,gF = pars.gF,wq = pars.wq,σC = pars.σC,σH = pars.σH,αWR = pars.αWR,αWR2 = pars.αWR2,αF = pars.αF[i],αH = pars.αH[i],αA = pars.αA[i],Γ=pars.Γ,β=pars.β)
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
        #αH = pars.αH[(pos+1):(pos+site_features.T[i])]
        pars_site = GetSitePars(pars,i)
        sname = site_list[i]
        #println(sname)
        Y = getfield(budget,sname)
        moms = getfield(moments,sname)
        π0 = site_features.π0[i,:,:]
        wght = getfield(wghts,sname)
        year_meas = site_features.year_meas[i]
        Abar = size(Y)[1]
        price = site_features.prices[i,1]
        moms_control = GetMoments(pars_site,Y[1,:,:,:,:],price,0,T,π0,year_meas)
        Qn += sum(wght[:,1].*(moms[:,1] .- moms_control).^2)

        for a = 2:site_features.n_arms[i]
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
            TE = moms[:,a] .- moms[:,1]
            TE_mod = moms_model .- moms_control
            Qn += sum(wght[:,a].*(TE .- TE_mod).^2)
            #Qn += sum(wght[:,a].*(moms_model .- moms[:,a]).^2)
        end

        # for a = 1:site_features.n_arms[i]
        #     WR = site_features.work_reqs[i,a]
        #     price = site_features.prices[i,a]
        #     abar = min(Abar,a)
        #     moms_model = GetMoments(pars_site,Y[abar,:,:,:,:],price,WR,T,π0,year_meas)
        #     Qn += sum(wght[:,a].*(moms[:,a] .- moms_model).^2)
        # end
    end
    return Qn
end

function CriterionSite(x,pars,vars,site_list,budget,moments,wghts,site_features,i)
    pars = UpdatePars(x,pars,vars)
    Qn = Real(0.)
    yb = site_features.yb[i]
    T = site_features.T[i]
    years = (yb+1-1991):(yb-1991+T)
    sname = site_list[i]
    Y = getfield(budget,sname)
    moms = getfield(moments,sname)
    π0 = site_features.π0[i,:,:]
    wght = getfield(wghts,sname)
    year_meas = site_features.year_meas[i]
    Abar = size(Y)[1]

    # do control group
    price = site_features.prices[i,1]
    moms_control = GetMoments(pars,Y[1,:,:,:,:],price,0,T,π0,year_meas)
    Qn += sum(wght[:,1] .* (moms[:,1] .- moms_control).^2)
    for a = 2:site_features.n_arms[i]
        WR = site_features.work_reqs[i,a]
        price = site_features.prices[i,a]
        abar = min(Abar,a)
        if site_features.time_limits[i,a]==1
            Y_I = budget[Symbol(sname,"_I")]
            TLlength = site_features.TLlength[i,a]
            moms_model = GetMomentsTimeLims(pars,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,π0,TLlength,year_meas)
        else
            moms_model = GetMoments(pars,Y[abar,:,:,:,:],price,WR,T,π0,year_meas)
        end
        TE = moms[:,a] .- moms[:,1]
        TE_mod = moms_model .- moms_control
        Qn += sum(wght[:,a].*(TE .- TE_mod).^2)
    end
    return Qn
end

function CriterionSite(x,g,pars,vars,site_list,budget,moments,wghts,site_features,i)
    f0 = CriterionSite(x,pars,vars,site_list,budget,moments,wghts,site_features,i)
    dF = ForwardDiff.gradient(x->CriterionSite(x,pars,vars,site_list,budget,moments,wghts,site_features,i),x)
    # don't have a better way to do this currently
    for j in findall(isnan.(dF))
        Δ = 1e-6
        x1 = copy(x)
        x1[j] += Δ
        f1 = CriterionSite(x1,pars,vars,site_list,budget,moments,wghts,site_features,i)
        dF[j] = (f1-f0)/Δ
    end
    g[:] = dF
    return f0
end

function FitSite(vars,site_list,budget,moments,wghts,site_features,i)
    T = site_features.T[i]
    np = (αc = 1, gN = 1, gF =1, αH = 1, αA = 1, σH = 1, σC = 1, wq = 1, αWR = 1, αWR2 = 1, αF = 1, αHT = 8,β=1)
    lb = (αc = 0., gN = -Inf, gF =-Inf, αH = -Inf, αA = -Inf, σH = 0., σC = 0., wq = 0., αWR = -Inf, αWR2 = -Inf, αF = -Inf,β=0.)
    ub = (αc = Inf, gN = Inf, gF =Inf, αH = Inf, αA = Inf, σH = Inf, σC = Inf, wq = Inf, αWR = Inf, αWR2 = Inf, αF = Inf,β=1.)
    αc = 1.
    gN = 0.
    gF = 0.
    αH = 0. #zeros(T) #zeros(35)
    αA = 0.
    β = 0.5
    σH = 1.
    σC = 1.
    wq = 2.
    αWR = 0.
    αWR2 = 0.
    αF = 0.
    αHT = zeros(8)
    Γ = zeros(18)
    pars = (np=np,lb=lb,ub=ub,αc=αc,gN=gN,gF=gF,αH=αH,αA=αA,β=β,σH=σH,σC=σC,wq=wq,αWR=αWR,αWR2=αWR2,αF=αF,Γ=Γ)
    opt,x0 = GetOptimization(pars,vars,site_list,budget,moments,wghts,site_features)
    min_objective!(opt,(x,g)->CriterionSite(x,g,pars,vars,site_list,budget,moments,wghts,site_features,i))
    res = optimize(opt,x0)
    pars = UpdatePars(res[2],pars,vars)
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
            moms_model[:,a] = GetMomentsTimeLims(pars,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,π0,TLlength,year_meas)
        else
            moms_model[:,a] = GetMoments(pars,Y[abar,:,:,:,:],price,WR,T,π0,year_meas)
        end

    end
    println(res[1])
    println(res[3])
    return pars,moms_model,moms
end

function FitSite(pars,vars,site_list,budget,moments,wghts,site_features,i)
    T = site_features.T[i]
    opt,x0 = GetOptimization(pars,vars,site_list,budget,moments,wghts,site_features)
    min_objective!(opt,(x,g)->CriterionSite(x,g,pars,vars,site_list,budget,moments,wghts,site_features,i))
    res = optimize(opt,x0)
    pars = UpdatePars(res[2],pars,vars)
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
        #moms_model[:,a] = GetMoments(pars,Y[abar,:,:,:,:],price,WR,T,π0,year_meas)
        if site_features.time_limits[i,a]==1
            Y_I = budget[Symbol(sname,"_I")]
            TLlength = site_features.TLlength[i,a]
            moms_model[:,a] = GetMomentsTimeLims(pars,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,π0,TLlength,year_meas)
        else
            moms_model[:,a] = GetMoments(pars,Y[abar,:,:,:,:],price,WR,T,π0,year_meas)
        end

    end
    println(res[1])
    println(res[3])
    return pars,moms_model,moms
end


function GetMomentsAll(pars,site_list,budget,moments,wghts,site_features)
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
        moms_model = []
        for a = 1:site_features.n_arms[i]
            WR = site_features.work_reqs[i,a]
            price = site_features.prices[i,a]
            abar = min(Abar,a)
            if site_features.time_limits[i,a]==1
                Y_I = budget[Symbol(sname,"_I")]
                TLlength = site_features.TLlength[i,a]
                append!(moms_model,[GetMomentsTimeLims(pars_site,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,π0,TLlength,year_meas,true)])
            else
                append!(moms_model,[GetMoments(pars_site,Y[abar,:,:,:,:],price,WR,T,π0,year_meas,true)])
            end
        end
        append!(moms_collect,[moms_model])
    end
    return (;zip(site_list,moms_collect)...)
end

function GetChoiceProbsAll(pars,site_list,budget,site_features)
    probs_collect = []
    NK = size(site_features.π0)[2]
    for i=1:length(site_list)
        yb = site_features.yb[i]
        T = site_features.T[i]
        #years = (yb+1-1991):(yb-1991+T)
        #pos = sum(site_features.T[1:i-1])
        #αH = pars.αH[(pos+1):(pos+site_features.T[i])]
        pars_site = GetSitePars(pars,i)
        sname = site_list[i]
        Y = getfield(budget,sname)
        Abar = size(Y)[1]
        PA = []
        PW = []
        PF = []
        for a = 1:site_features.n_arms[i]
            WR = site_features.work_reqs[i,a]
            price = site_features.prices[i,a]
            abar = min(Abar,a)
            if site_features.time_limits[i,a]==1
                Y_I = budget[Symbol(sname,"_I")]
                TLlength = site_features.TLlength[i,a]
                pA,pWork,pF = GetDynamicProbs(pars_site,Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],price,WR,T,NK,TLlength)
            else
                pA,pWork,pF = GetStaticProbs(pars_site,Y[abar,:,:,:,:],price,WR,T,NK)
            end
            append!(PA,[pA])
            append!(PW,[pWork])
            append!(PF,[pF])
        end
        P = (pA = PA, pWork = PW, pF = PF)
        append!(probs_collect,[P])
    end
    return (;zip(site_list,probs_collect)...)
end


function Criterion(x,g,pars,vars,site_list,budget,moments,wghts,site_features)
    f0 = Criterion(x,pars,vars,site_list,budget,moments,wghts,site_features)
    dF = ForwardDiff.gradient(x->Criterion(x,pars,vars,site_list,budget,moments,wghts,site_features),x)
    g[:] = dF
    return f0
end

function GetMoments(pars,Y,price,WR,T,π0,year_meas,wrapped=false)
    NK = size(π0)[1]
    pA,pWork,pF = GetStaticProbs(pars,Y,price,WR,T,NK)
    EA,EH,Inc = GetAggMoments(pA,pWork,π0,Y)
    Care = GetMeanCare(year_meas,0,9,pA,pWork,pF,π0)
    if wrapped
        return (Part = EA,LFP=EH,Inc = Inc,Care=Care)
    else
        return [EA; EH; Care; Inc]
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
        return [EA; EH; Care; Inc]
    end
end

#
# function InspectTreatFit(model_moments,data_moments,site_features,site_list)
#     colors = ["red","green","blue"]
#     for i=1:8
#         T = site_features.T[i]
#         moms0 = getfield(model_moments,site_list[i])
#         moms1 = getfield(data_moments,site_list[i])
#         for a=2:site_features.n_arms[i]
#             figure("AFDC")
#             subplot(2,4,i)
#             title(String(site_list[i]))
#             plot(moms0[1:T,a] .- moms0[1:T,1],color=colors[a])
#             plot(moms1[1:T,a] .- moms1[1:T,1],color=colors[a],linestyle="--")
#             figure("LFP")
#             subplot(2,4,i)
#             title(String(site_list[i]))
#             plot(moms0[T+1:2*T,a] .- moms0[T+1:2*T,1],color=colors[a])
#             plot(moms1[T+1:2*T,a] .- moms1[T+1:2*T,1],color=colors[a],linestyle="--")
#             figure("Scatter")
#             scatter(moms0[:,a] .- moms0[:,1],moms1[:,a] .- moms1[:,1],color="blue")
#         end
#     end
# end

# function InspectModelFit(model_moments,data_moments,site_features,site_list)
#     colors = ["red","green","blue"]
#     for i=1:8
#         T = site_features.T[i]
#         moms0 = getfield(model_moments,site_list[i])
#         moms1 = getfield(data_moments,site_list[i])
#         for a=1:site_features.n_arms[i]
#             figure("AFDC")
#             subplot(2,4,i)
#             title(String(site_list[i]))
#             plot(moms0[1:T,a],color=colors[a])
#             plot(moms1[1:T,a],color=colors[a],linestyle="--")
#             figure("LFP")
#             subplot(2,4,i)
#             title(String(site_list[i]))
#             plot(moms0[T+1:2*T,a],color=colors[a])
#             plot(moms1[T+1:2*T,a],color=colors[a],linestyle="--")
#         end
#     end
# end

function InspectModelFit(model_moments,data_moments,site_features,site_list)
    colors = ["blue","green","red"]
    for i=1:8
        T = site_features.T[i]
        sname = site_list[i]
        moms0 = getfield(data_moments,sname)
        moms1 = getfield(model_moments,sname)
        for a=1:site_features.n_arms[i]
            for v in [:Part,:LFP,:Inc]
                figure(String(v))
                subplot(2,4,i)
                title(String(sname))
                plot(moms0[a][v],color=colors[a])
                plot(moms1[a][v],color=colors[a],linestyle="--")
            end
        end
    end
end

function InspectTreatFit(model_moments,data_moments,site_features,site_list)
    colors = ["blue","green","red"]
    for i=1:8
        T = site_features.T[i]
        sname = site_list[i]
        moms0 = getfield(data_moments,sname)
        moms1 = getfield(model_moments,sname)
        for a=2:site_features.n_arms[i]
            for v in [:Part,:LFP,:Inc]
                figure(String(v))
                subplot(2,4,i)
                title(String(sname))
                plot(moms0[a][v] .- moms0[1][v],color=colors[a])
                plot(moms1[a][v] .- moms1[1][v],color=colors[a],linestyle="--")
                if v==:Inc
                    color_="pink"
                else
                    color_="blue"
                end
                figure("Scatter")
                scatter(moms0[a][v] .- moms0[1][v],moms1[a][v] .- moms1[1][v],color=color_)
            end
        end
    end
end
