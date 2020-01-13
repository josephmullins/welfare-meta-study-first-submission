# routines for estimating the baseline model
using ForwardDiff

function GetParameters(x,num_sites)
    # universal: αc,Γ,gN,gF,wq,σC,σH,αWR
    # site-specific: αH,αA
    αc,gN,gF,wq,σC,σH,αWR,αF = x[1:8]
    αH = x[9:(8+num_sites)]
    αA = x[(9+num_sites):(8+2*num_sites)]
    return (αc=αc,gN = ones(2)*gN, gF = ones(2)*gF,wq = wq,σC = σC,σH = σH,αWR = αWR,αF = αF, αA = αA,αH=αH)
end

# site features contains: T,n_arms,work_requirement,π0,prices
function Criterion(x,site_list,budget,moments,wghts,site_features)
    num_sites = length(site_list)
    pars = GetParameters(x,num_sites)
    Qn = 0
    for i=1:site_list
        pars_site = (αc = pars.αc,gN = pars.gN,gF = pars.gF,wq = pars.wq,σC = pars.σC,σH = pars.σH,αWR = pars.αWR,αF = pars.αF,αH = pars.αH[i],αA = pars.αA[i])
        sname = site_list[i]
        Y = getfield(budget,sname)
        moms = getfield(moments,sname)
        π0 = site_features.π0[i,:,:]
        T = site_features.T[i]
        wght = getfield(wghts,sname)
        for a = 1:site_features.n_arms[i]
            WR = site_features.work_reqs[i,a]
            price = site_features.prices[i,a]
            moms_model = GetMoments(pars_site,Y,price,WR,T,π0)
            Qn += sum(wght[:,a].*(moms[:,a] .- moms_model))
        end
    end
    return Qn
end

function Criterion(x,g,site_list,budget,moments,wghts,site_features)
    f0 = Criterion(x,site_list,budget,moments,wghts,site_features)
    dF = ForwardDiff.gradient(x->Criterion(x,site_list,budget,work_reqs,site_years,price,moments,wghts),x)
    g[:] = dF
    return f0
end

function GetMoments(pars,Y,price,WR,T,π0)
    NK = size(π0)[1]
    pA,pWork,pF = GetStaticProbs(pars,Y,price,WR,T,NK)
    EA,EH = GetAggMoments(pA,pWork,π0)
    Care = GetMeanCare(T,0,9,pA,pWork,π0)
    return [EA; EH; Care]
end
