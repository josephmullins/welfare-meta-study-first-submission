
function GetEffects(pars,Y,Y_fs,price,T,NK,π0)
    EffectA = zeros(2,T)
    EffectH = zeros(2,T)
    ElastA = zeros(2,T)
    ElastH = zeros(2,T)
    ElastF = zeros(2,T)

    pA,pWork,pF = GetStaticProbs(pars,Y,price,0,T,NK)
    EA0,EH0 = GetAggMoments(pA,pWork,π0)
    # four counterfactuals to work through

    # 1 - Work requirement
    pA1,pWork1,pF1 = GetStaticProbs(pars,Y,price,1,T,NK)
    EA1,EH1 = GetAggMoments(pA1,pWork1,π0)
    EffectA[1,:] = EA1 .- EA0
    EffectH[1,:] = EH1 .- EH0

    # 2 - time limit
    pA1,pWork1,pF1 = GetDynamicProbs(pars,Y,Y_fs,price,0,T,NK,5)
    EA1,EH1,c = GetDynamicMoments(pA1,pWork1,pF1,π0,1,0,9)
    EffectA[2,:] = EA1 .- EA0
    EffectH[2,:] = EH1 .- EH0

    # 3 - 10% increase in post-transfer wages
    Yalt = copy(Y)
    Yalt[:,:,:,2] += 0.1*(Y[:,:,:,2] .- Y[:,:,:,1]) #<- 10% increase in post-tax wage
    pA1,pWork1,pF1 = GetStaticProbs(pars,Yalt,price,0,T,NK)
    #EA1,EH1,c = GetDynamicMoments(pA1,pWork1,pF1,π0,3,0,9)
    ElastA[1,:],ElastH[1,:],ElastF[1,:] = GetElasticities((pA,pA1),(pWork,pWork1),(pF,pF1),π0)

    # 4 - 10% increase in price of care
    pA1,pWork1,pF1 = GetStaticProbs(pars,Y,1.1*price,0,T,NK)
    EA1,EH1 = GetAggMoments(pA1,pWork1,π0)
    ElastA[2,:],ElastH[2,:],ElastF[2,:] = GetElasticities((pA,pA1),(pWork,pWork1),(pF,pF1),π0)

    return EffectA,EffectH,ElastA,ElastH,ElastF
end

function GetChildEffects(pars_prod,pars,Y,Y_fs,price,T,NK,π0)
    pA,pWork,pF = GetStaticProbs(pars,Y,price,0,T,NK)
    # four counterfactuals to work through
    TH = zeros(T,3,NK,17)
    for t=1:T
        TH[t,1,:,:] = GetChildOutcomesStatic(t,pA,pWork,pF,Y,pars_prod,price,pars.wq)
    end

    # 1 - Work requirement
    pA1,pWork1,pF1 = GetStaticProbs(pars,Y,price,1,T,NK)
    for t=1:T
        TH[t,2,:,:] = GetChildOutcomesStatic(t,pA1,pWork1,pF1,Y,pars_prod,price,pars.wq)
    end

    # 2 - Extra thousand dollars per year
    Yalt = copy(Y);
    Yalt .+= 1000/52
    pA1,pWork1,pF1 = GetStaticProbs(pars,Yalt,price,1,T,NK)
    for t=1:T
        TH[t,3,:,:] = GetChildOutcomesStatic(t,pA1,pWork1,pF1,Y,pars_prod,price,pars.wq)
    end

    # 2 - time limit
    #pA1,pWork1,pF1 = GetDynamicProbs(pars,Y,Y_fs,price,0,T,NK,5)
    #EA1,EH1,c = GetDynamicMoments(pA1,pWork1,pF1,π0,1,0,9)
    #EffectA[2,:] = EA1 .- EA0
    #EffectH[2,:] = EH1 .- EH0
    TE = zeros(T,2,2)
    for t=1:T
        for j=2:3
            ii = 1:6
            TE[t,1,j-1] = sum(π0[:,ii].*(TH[t,j,:,ii] .- TH[t,1,:,ii]))/sum(π0[:,ii])
            ii = 7:13
            TE[t,2,j-1] = sum(π0[:,ii].*(TH[t,j,:,ii] .- TH[t,1,:,ii]))/sum(π0[:,ii])
        end
    end
    return TE
end
# plan for calculating elasticity of formal care! When age<=5? Most useful.
# need to make a decision about it.
function GetElasticities(pA,pWork,pF,π0)
	T = size(pA[1])[1]
	NK = size(pA[1])[2]
	EA = zeros(Real,T)
	EH = zeros(Real,T)
    EF = zeros(Real,T)
	for t=1:T
        denom = 0
        numerator = 0
		for a0 = 0:16
			for nk=1:NK
				EA[t] += π0[nk,a0+1]*(pA[2][t,nk,a0+1]-pA[1][t,nk,a0+1])/pA[1][t,nk,a0+1]
                pW0 = pWork[1][t,nk,a0+1,1] + pA[1][t,nk,a0+1]*(pWork[1][t,nk,a0+1,2]-pWork[1][t,nk,a0+1,1])
                pW1 = pWork[2][t,nk,a0+1,1] + pA[2][t,nk,a0+1]*(pWork[2][t,nk,a0+1,2]-pWork[2][t,nk,a0+1,1])
				EH[t] += π0[nk,a0+1]*(pW1-pW0)/pW0
                mass1 = pA[1][t,nk,1+a0]*pWork[1][t,nk,1+a0,2]
                mass0 = (1-pA[1][t,nk,1+a0])*pWork[1][t,nk,1+a0,1]
                e1 = (pF[2][t,nk,1+a0,2]-pF[1][t,nk,1+a0,2])/pF[1][t,nk,1+a0,2]
                e0 = (pF[2][t,nk,1+a0,1]-pF[1][t,nk,1+a0,1])/pF[1][t,nk,1+a0,1]
                if a0<=9
                    numerator += π0[nk,a0+1]*(mass1*e1+mass0*e0)
                    denom += π0[nk,a0+1]*(mass0+mass1)
                end
			end
		end
        EF[t] = numerator/denom
	end
	return EA,EH,EF
end


function GetEffects(x,pars,budget,site_features)
    pmod = GetModParameters(x,pars)
    nsite = length(site_features.site_list)
    EffectA = zeros(nsite,2,5)
    EffectH = zeros(nsite,2,5)
    ElastA = zeros(nsite,2,5)
    ElastH = zeros(nsite,2,5)
    ElastF = zeros(nsite,2,5)
    for i=1:nsite
        yb = site_features.yb[i]
        T = site_features.T[i]
        pars_site = GetSitePars(pmod,i)
        sname = site_features.site_list[i]
        #println(sname)
        Y = getfield(budget,sname)
        Y_I = budget[Symbol(sname,"_I")]
        moms = getfield(moments,sname)
        π0 = site_features.π0[i,:,:]
        year_meas = site_features.year_meas[i]
        Abar = size(Y)[1]
        price = site_features.prices[i,1]
        #moms_control = GetMoments(pars_site,Y[1,:,:,:,:],price,0,T,π0,year_meas)
        # GetEffects(pars,Y,Y_fs,price,T,NK)
        NK = size(π0)[1]
        EffectA[i,:,1:T],EffectH[i,:,1:T],ElastA[i,:,1:T],ElastH[i,:,1:T],ElastF[i,:,1:T] = GetEffects(pars_site,Y[1,:,:,:,:],Y_I[1,:,:,:,:],price,T,NK,π0)
    end
    return EffectA,EffectH,ElastA,ElastH,ElastF
end


function GetChildEffects(x,pmod,budget,site_features)
    #prod = GetProdParameters(x,prod_pars)
    δI = x[1:2]
    δθ = x[3]
    gN = x[4:5]
    gF = x[6:7]
    prod_pars = (δI = x[1:2],δθ=x[3],gN = x[4:5],gF = x[6:7])

    #pmod = GetModParameters(x,pars)
    nsite = length(site_features.site_list)
    TE = zeros(nsite,5,2,2)
    for i=1:nsite
        yb = site_features.yb[i]
        T = site_features.T[i]
        pars_site = GetSitePars(pmod,i)
        sname = site_features.site_list[i]
        #println(sname)
        Y = getfield(budget,sname)
        Y_I = budget[Symbol(sname,"_I")]
        moms = getfield(moments,sname)
        π0 = site_features.π0[i,:,:]
        year_meas = site_features.year_meas[i]
        Abar = size(Y)[1]
        price = site_features.prices[i,1]
        #moms_control = GetMoments(pars_site,Y[1,:,:,:,:],price,0,T,π0,year_meas)
        # GetEffects(pars,Y,Y_fs,price,T,NK)
        NK = size(π0)[1]
        TE[i,1:T,:,:] = GetChildEffects(prod_pars,pars_site,Y[1,:,:,:,:],Y_I[1,:,:,:,:],price,T,NK,π0)
    end
    return TE
end
