
# we need to think about how to do this a little better
function ProdPars(nmeas)
    np = (δI = 2,gN = 2, gF = 2, δθ=1, λ=nmeas-1)
    lb = (δI = -Inf*ones(2),gN = -Inf*ones(2),gF = -Inf*ones(2), δθ = 0., λ = -Inf*ones(nmeas-1))
    ub = (δI = Inf*ones(2),gN = Inf*ones(2),gF = Inf*ones(2), δθ = 1., λ = Inf*ones(nmeas-1))
    δI = 10. *ones(2)
    gN = zeros(2) #log.([1.5,1.5])
    gF = zeros(2) #log.([0.7,0.7])
    δθ = 0.5
    #λ = zeros(nmeas-1)
    λ = ones(nmeas-1)
    #λ = [11.,15.,0.]
    #λ = [1.,-1.,1.,-1.,1.,-1.,-1.]
    #λ = [2.3,-2.48,2.85,-2.9,1.29]
    return (np=np,lb=lb,ub=ub,δI=δI,gN=gN,gF=gF,δθ=δθ,λ=λ)
end

function TreatmenEffects(TH,π0,a0,a1,arms)
    nm = length(a0)
    means = zeros(Real,nm)
    for m=1:nm
        #println(m)
        arm = arms[m]
        ii = (a0[m]+1):(min(a1[m],16)+1)
        means[m] = sum(π0[:,ii].*(TH[arm,:,ii] .- TH[1,:,ii]))/sum(π0[:,ii])
    end
    return means
end

function GetOptimization(vars,pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features)
    np = 0;
    #vars = [:gN,:gF,:δI,:δθ,:λ]
	for s in vars
		np += pars_prod.np[s]
	end
	x0 = zeros(Real,np)
	lb = zeros(np)
	ub = zeros(np)
    opt = Opt(:LD_LBFGS,np)
    pos = 1
	for s in vars
		np = pars_prod.np[s]
		if typeof(pars_prod.lb[s])==Float64
			x0[pos] = getfield(pars_prod,s)
			lb[pos] = pars_prod.lb[s]
			ub[pos] = pars_prod.ub[s]
		else
			println(s)
			xrange = pos:(pos+np-1)
			x0[xrange] = getfield(pars_prod,s)
			lb[xrange] = pars_prod.lb[s]
			ub[xrange] = pars_prod.ub[s]
		end
		pos += np
	end
	lower_bounds!(opt,lb)
	upper_bounds!(opt,ub)
	min_objective!(opt,(x,g)->ProductionCriterion(x,g,vars,pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features))
    return opt,x0
end

function ProductionCriterion(x,vars,pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features)
    #vlist = [:gN,:gF,:δI,:δθ,:λ]
    pars_prod = UpdatePars(x,pars_prod,vars)
    return ProductionCriterion(pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features)
end

function ProductionCriterion(x,g,vars,pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features)
    f0 = ProductionCriterion(x,vars,pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features)
    dF = ForwardDiff.gradient(x->ProductionCriterion(x,vars,pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features),x)
    g[:] = dF
    return f0
end


function ProductionCriterion(pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features)
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
        TH[1,:,:] =  GetChildOutcomesStatic(year_meas,P.pA[1],P.pWork[1],P.pF[1],Y[1,:,:,:,:],pars_prod,price,pars.wq)

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
                TH[a,:,:] = GetChildOutcomesDynamic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],pars_prod,price,pars.wq)
            else
                TH[a,:,:] = GetChildOutcomesStatic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[abar,:,:,:,:],pars_prod,price,pars.wq)
            end
        end
        TE = TreatmenEffects(TH,π0,moms.a0,moms.a1,moms.arm)
        for k=1:size(moms.TE)[2]
            Qn += sum(moms.wght[:,k].*(λ[k]*TE .- moms.TE[:,k]).^2)
        end
    end
    return Qn
end

function GetTreatmentEffects(pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features)
    Qn = 0
    NK=3 # might want to update this later
    λ = [1;pars_prod.λ]
    colors = ["blue","red","green","pink","purple","orange","black"]
    TE_collect = []
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
        TH[1,:,:] =  GetChildOutcomesStatic(year_meas,P.pA[1],P.pWork[1],P.pF[1],Y[1,:,:,:,:],pars_prod,price,pars.wq)

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
                TH[a,:,:] = GetChildOutcomesDynamic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],pars_prod,price,pars.wq)
            else
                TH[a,:,:] = GetChildOutcomesStatic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[abar,:,:,:,:],pars_prod,price,pars.wq)
            end
        end
        TE = TreatmenEffects(TH,π0,moms.a0,moms.a1,moms.arm)
        append!(TE_collect,[TE])
    end
    return (;zip(site_list,TE_collect)...)
end

function InspectTreatFitProduction!(pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features)
    Qn = 0
    NK=3 # might want to update this later
    λ = [1;pars_prod.λ]
    colors = ["blue","red","green","pink","purple","orange","black","magenta","grey"]
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
        TH[1,:,:] =  GetChildOutcomesStatic(year_meas,P.pA[1],P.pWork[1],P.pF[1],Y[1,:,:,:,:],pars_prod,price,pars.wq)

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
                TH[a,:,:] = GetChildOutcomesDynamic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[abar,:,:,:,:],Y_I[abar,:,:,:,:],pars_prod,price,pars.wq)
            else
                TH[a,:,:] = GetChildOutcomesStatic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[abar,:,:,:,:],pars_prod,price,pars.wq)
            end
        end
        TE = TreatmenEffects(TH,π0,moms.a0,moms.a1,moms.arm)
        #plt.scatter(TE,moms.TE[:,1],color="blue")
        for k=1:size(moms.TE)[2]
            #Qn += sum(moms.wght[:,k].*(λ[k]*TE .- moms.TE[:,k]).^2)
            ii = moms.wght[:,k].>0
            plt.scatter(λ[k]*TE[ii].*sqrt.(moms.wght[ii,k]), moms.TE[ii,k],color=colors[k])
        end
    end
end
