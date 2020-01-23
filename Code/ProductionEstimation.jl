
# we need to think about how to do this a little better
function ProdPars(nmeas)
    np = (δI = 2,gN = 2, gF = 2, δθ=1, λ=nmeas-1
    lb = (δI = zeros(2),gN = -Inf*ones(2),gF = -Inf*ones(2), δθ = 0., λ = -Inf*ones(nmeas-1))
    lb = (δI = Inf*ones(2),gN = Inf*ones(2),gF = Inf*ones(2), δθ = Inf, λ = Inf*ones(nmeas-1))
    δI = 0.1
    gN = zeros(2)
    gF = zeros(2)
    δθ = 0.9
    λ = zeros(nmeas-1)
    return (np=np,lb=lb,δI=δI,gN=gN,gF=gF,δθ=δθ,λ=λ)
end

function TreatmenEffects(TH,π0,a0,a1,arms)
    nm = length(a0)
    means = zeros(Real,2)
    for m=1:nm
        arm = arms[m]
        ii = (a0[m]+1):(a1[m]+1)
        means[m] = sum(π0[:,ii].*(TH[arm,:,ii] .- TH[1,:,ii]))/sum(π0[:,ii])
    end
    return means
end


function ProductionCriterion(pars_prod,pars,ChoiceProbs,site_list,budget,moments,site_features)
    Qn = 0
    NK=3 # might want to update this later
    λ = [1;pars_prod.λ]
    for i=1:length(site_list)
        sname = site_list[i]
        moms = getfield(moments,sname) #<- structure: age0,age1,1
        nmom = size(moments)[1]
        π0 = site_features.π0[i,:,:]
        #wght = getfield(wghts,sname)
        year_meas = site_features.year_meas[i]
        Y = getfield(budget,sname)
        P = getfield(budget,sname)
        TH = zeros(Real,site_features.n_arms[i],NK)

        TH[1,:,:] =  GetChildOutcomesStatic(year_meas,P.pA[1],P.pWork[1],P.pF[1],Y[1,:,:,:,:],pars_prod,pars)

        Abar = size(Y)[1]
        price = site_features.prices[i,1]
        #moms_control = GetMoments(pars_site,Y[1,:,:,:,:],price,0,T,π0,year_meas)
        #Qn += sum(wght[:,1].*(moms[:,1] .- moms_control).^2)
        for a = 2:site_features.n_arms[i]
            WR = site_features.work_reqs[i,a]
            price = site_features.prices[i,a]
            abar = min(Abar,a)
            if site_features.time_limits[i,a]==1
                Y_I = budget[Symbol(sname,"_I")]
                TLlength = site_features.TLlength[i,a]
                TH[a,:,:] = GetChildOutcomesDynamic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[a,:,:,:,:],Y_I[a,:,:,:,:],pars_prod,pars)
            else
                TH[a,:,:] = GetChildOutcomesStatic(year_meas,P.pA[a],P.pWork[a],P.pF[a],Y[a,:,:,:,:],pars_prod,pars)
            end
        end
        TE = TreatmenEffects(TH,π0,moms.a0,moms.a1,moms.arms)
        for k=1:size(moms.TE)[2]
            Qn += sum(moms.wght[:,k].*(λ[k]*TE .- moms.TE[:,k]).^2)
        end
    end
    return Qn
end
