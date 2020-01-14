# for now we're going to assume annual frequency (for simplicity)

# this function calculates the choice probabilities, and produces the Emax function
function ChoiceProb(pars,Y,price,nk,age,CV,WR)
	pWork = zeros(Real,2) #<- work choice probabilities
	pF = zeros(Real,2) #<- probability of formal care given work and welfare choice
	vA = zeros(Real,2) #<- value from participation choice
	for p=0:1
		if age<=17 #<- assume quarterly
			if age<=5 #<- assume quarterly
				gN,gF = pars.Γ[age+1]*pars.gN[1],pars.Γ[age+1]*pars.gF[1]
			else
				gN,gF = pars.Γ[age+1]*pars.gN[1],pars.Γ[age+1]*pars.gF[1]
			end
			uF = (pars.αc+pars.Γ[age+1])*log(Y[nk+1,1+p,2]+(112-30)*pars.wq-price) + pars.αF - gF
			uN = (pars.αc+pars.Γ[age+1])*log(Y[nk+1,1+p,2]+(112-30)*pars.wq) - gN
			pF[p+1] = 1/(1+exp((uN-uF)/pars.σC)) #<- probability of formal care
			vW1 = pars.σC*log(exp(uN/pars.σC)+exp(uF/pars.σC)) - pars.αH
			vW0 = (pars.αc+pars.Γ[age+1])*log(Y[nk+1,1+p,1] + 112*pars.wq) - pars.αWR*WR*p
			pWork[p+1] = 1/(1+exp((vW0-vW1)/pars.σH))
			vA[p+1] = pars.σH*log(exp(vW0/pars.σH)+exp(vW1/pars.σH)) + pars.β*CV[p+1] - pars.αA*p
		else
			vW1 = pars.αc*log(Y[1,1+p,2]+(112-30)*pars.wq) - pars.αH
			vW0 = pars.αc*log(Y[1,1+p,1]+(112)*pars.wq)
			pWork[p+1] = 1/(1+exp((vW0-vW1)/pars.σH))
			vA[p+1] = pars.σH*log(exp(vW0/pars.σH)+exp(vW1/pars.σH)) + pars.β*CV[p+1] - pars.αA*p
		end
	end
	pA = 1/(1+exp(vA[1]-vA[2]))
	V = log(exp(vA[1])+exp(vA[2]))
	return pA,pWork,pF,V
end

# NK is the max number of kids
function GetStaticProbs(pars,Y,price,WR,T,NK)
	pWork = zeros(Real,T,NK,17,2)
	Tbar = size(Y)[1]
	pF = zeros(Real,T,NK,17,2)
	pA = zeros(Real,T,NK,17)
	CV = zeros(2) #<- static, don't bother with this.
	for t=1:T
		for k=1:NK
			for a0 = 0:16
				a = a0+t-1
				tbar = min(Tbar,t)
				pA[t,k,a0+1],pWork[t,k,a0+1,:],pF[t,k,a0+1,:],V = ChoiceProb(pars,Y[tbar,:,:,:],price,k,a,CV,WR)
			end
		end
	end
	return pA,pWork,pF
end

function GetAggMoments(pA,pWork,π0)
	T = size(pA)[1]
	NK = size(pA)[2]
	EA = zeros(Real,T)
	EH = zeros(Real,T)
	for t=1:T
		for a0 = 0:16
			for nk=1:NK
				EA[t] += π0[nk,a0+1]*pA[t,nk,a0+1]
				EH[t] += π0[nk,a0+1]*(pWork[t,nk,a0+1,1] + pA[t,nk,a0+1]*(pWork[t,nk,a0+1,2]-pWork[t,nk,a0+1,1]))
			end
		end
	end
	return EA,EH
end

# what's available?
# this will be the mean fraction among those working
# a0 is the lower bound and a1 is the upper bound at entry
function GetMeanCare(year,a0,a1,pA,pWork,pF,π0)
	NK = size(pA)[2]
	denom = 0
	numerator = 0
	for a=a0:a1
		for nk=1:NK
			mass1 = pA[year,nk,1+a]*pWork[year,nk,1+a,2]
			mass0 = (1-pA[year,nk,1+a])*pWork[year,nk,1+a,1]
			numerator += π0[nk,a+1]*(mass1*pF[year,nk,1+a,2]+mass0*pF[year,nk,1+a,1])
			denom += π0[nk,a+1]*(mass0+mass1)
		end
	end
	return numerator/denom
end
