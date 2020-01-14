# routines for estimating the baseline model
using ForwardDiff

function GetParameters(x,num_sites,Γ)
    # universal: αc,Γ,gN,gF,wq,σC,σH,αWR
    # site-specific: αH,αA
    αc,gN,gF,wq,σC,σH,αWR,αF = x[1:8]
    αH = x[9:(8+num_sites)]
    αA = x[(9+num_sites):(8+2*num_sites)]
    return (αc=αc,gN = ones(2)*gN, gF = ones(2)*gF,wq = wq,σC = σC,σH = σH,αWR = αWR,αF = αF, αA = αA,αH=αH,Γ=Γ)
end

# site features contains: T,n_arms,work_requirement,π0,prices,year_meas
function Criterion(x,site_list,budget,moments,wghts,site_features,Γ)
    num_sites = length(site_list)
    pars = GetParameters(x,num_sites,Γ)
    Qn = CriterionP(pars,site_list,budget,moments,wghts,site_features)
    return Qn
end

function CriterionP(pars::NamedTuple,site_list,budget,moments,wghts,site_features)
    Qn = 0
    for i=1:length(site_list)
        pars_site = (αc = pars.αc,gN = pars.gN,gF = pars.gF,wq = pars.wq,σC = pars.σC,σH = pars.σH,αWR = pars.αWR,αF = pars.αF,αH = pars.αH[i],αA = pars.αA[i],Γ=pars.Γ,β=0.)
        sname = site_list[i]
        Y = getfield(budget,sname)
        moms = getfield(moments,sname)
        π0 = site_features.π0[i,:,:]
        T = site_features.T[i]
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

function GetMomentsAll(pars::NamedTuple,site_list,budget,moments,wghts,site_features)
    moms_collect = []
    for i=1:length(site_list)
        pars_site = (αc = pars.αc,gN = pars.gN,gF = pars.gF,wq = pars.wq,σC = pars.σC,σH = pars.σH,αWR = pars.αWR,αF = pars.αF,αH = pars.αH[i],αA = pars.αA[i],Γ=pars.Γ,β=0.)
        sname = site_list[i]
        Y = getfield(budget,sname)
        moms = getfield(moments,sname)
        π0 = site_features.π0[i,:,:]
        T = site_features.T[i]
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


function Criterion(x,g,site_list,budget,moments,wghts,site_features,Γ)
    f0 = Criterion(x,site_list,budget,moments,wghts,site_features,Γ)
    dF = ForwardDiff.gradient(x->Criterion(x,site_list,budget,moments,wghts,site_features,Γ),x)
    g[:] = dF
    return f0
end

# this function is wrong! We have to increase the age of each child as they grow (though π0?)
function GetMoments(pars,Y,price,WR,T,π0,year_meas)
    NK = size(π0)[1]
    pA,pWork,pF = GetStaticProbs(pars,Y,price,WR,T,NK)
    EA,EH = GetAggMoments(pA,pWork,π0)
    Care = GetMeanCare(year_meas,0,9,pA,pWork,pF,π0)
    return [EA; EH; Care]
end
