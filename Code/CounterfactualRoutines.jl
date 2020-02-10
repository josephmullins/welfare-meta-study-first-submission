
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
                numerator += π0[nk,a0+1]*(mass1*e1+mass0*e0)
                denom += π0[nk,a0+1]*(mass0+mass1)
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
        #Y_I = budget[Symbol(sname,"_I")]
        moms = getfield(moments,sname)
        π0 = site_features.π0[i,:,:]
        year_meas = site_features.year_meas[i]
        Abar = size(Y)[1]
        price = site_features.prices[i,1]
        #moms_control = GetMoments(pars_site,Y[1,:,:,:,:],price,0,T,π0,year_meas)
        # GetEffects(pars,Y,Y_fs,price,T,NK)
        NK = size(π0)[1]
        EffectA[i,:,1:T],EffectH[i,:,1:T],ElastA[i,:,1:T],ElastH[i,:,1:T],ElastF[i,:,1:T] = GetEffects(pars_site,Y[1,:,:,:,:],Y[1,:,:,:,:],price,T,NK,π0)
    end
    return EffectA,EffectH,ElastA,ElastH,ElastF
end
