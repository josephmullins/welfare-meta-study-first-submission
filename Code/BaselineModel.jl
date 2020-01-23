# for now we're going to assume annual frequency (for simplicity)

function IncVal(v0,v1,σ)
	if v1>v0
		return v1 + σ*log(1+exp((v0-v1)/σ))
	else
		return v0 + σ*log(1+exp((v1-v0)/σ))
	end
end

# this function calculates the choice probabilities, and produces the Emax function
function ChoiceProb(pars,Y,price,t,nk,age,CV,WR)
	pWork = zeros(Real,2) #<- work choice probabilities
	pF = zeros(Real,2) #<- probability of formal care given work and welfare choice
	vA = zeros(Real,2) #<- value from participation choice
	for p=0:1
		if age<=17 #<-
			if age<=5 #<-
				gN,gF = pars.Γ[age+1]*pars.gN[1],pars.Γ[age+1]*pars.gF[1]
			else
				gN,gF = pars.Γ[age+1]*pars.gN[1],pars.Γ[age+1]*pars.gF[1]
			end
			uF = (pars.αc+pars.Γ[age+1])*log(Y[t,nk+1,1+p,2]+(112-30)*pars.wq-price) + pars.αF - gF
			uN = (pars.αc+pars.Γ[age+1])*log(Y[t,nk+1,1+p,2]+(112-30)*pars.wq) - gN
			pF[p+1] = 1/(1+exp((uN-uF)/pars.σC)) #<- probability of formal care
			vW1 = IncVal(uN,uF,pars.σC) - pars.αH - WR*p*pars.αWR2
			vW0 = (pars.αc+pars.Γ[age+1])*log(Y[t,nk+1,1+p,1] + 112*pars.wq) - pars.αWR*WR*p
			pWork[p+1] = 1/(1+exp((vW0-vW1)/pars.σH))
			vA[p+1] = IncVal(vW0,vW1,pars.σH) + pars.β*CV[p+1] - pars.αA*p
		else
			vW1 = pars.αc*log(Y[t,1,1+p,2]+(112-30)*pars.wq) - pars.αH
			vW0 = pars.αc*log(Y[t,1,1+p,1]+(112)*pars.wq)
			pWork[p+1] = 1/(1+exp((vW0-vW1)/pars.σH))
			vA[p+1] = IncVal(vW0,vW1,pars.σH) + pars.β*CV[p+1] - pars.αA*p
		end
	end
	pA = 1/(1+exp(vA[1]-vA[2]))
	V = IncVal(vA[1],vA[2],1.)
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
				pA[t,k,a0+1],pWork[t,k,a0+1,:],pF[t,k,a0+1,:],V = ChoiceProb(pars,Y,price,tbar,k,a,CV,WR)
			end
		end
	end
	return pA,pWork,pF
end

# one thing is missing: the budget when someone is ineligible (just food stamps)
# can we recover this?
# problem ends when child turns 18

function GetDynamicProbs(pars,Y,Y_fs,price,WR,T,NK,TLlength)
	Tbar = size(Y)[1]
	pWork = zeros(Real,T,NK,17,TLlength+1,2)
	pF = zeros(Real,T,NK,17,TLlength+1,2)
	pA = zeros(Real,T,NK,17,TLlength+1)
	for a0=0:16
		for nk=1:NK
			Tprob = max(18-a0,T)
			V = zeros(Real,Tprob+1,TLlength+1)
			for t=Tprob:-1:1
				age = a0+t-1
				tbar = min(Tbar,t)
				pA_ = zeros(Real,TLlength+1)
				pW_ = zeros(Real,TLlength+1,2)
				pF_ = zeros(Real,TLlength+1,2)
				for w=1:TLlength
					cv = V[t+1,[w,w+1]]
					pA_[w],pW_[w,:],pF_[w,:],V[t,w] = ChoiceProb(pars,Y,price,tbar,nk,age,cv,WR)
				end
				w_ = TLlength+1
				cv = [V[t+1,w_],V[t+1,w_]]
				pA_[w_],pW_[w_,:],pF_[w_,:],V[t,w_] = ChoiceProb(pars,Y_fs,price,tbar,nk,age,cv,WR)
				if t<=T
					pA[t,nk,a0+1,:] .= pA_
					pWork[t,nk,a0+1,:,:] .= pW_
					pF[t,nk,a0+1,:,:] .= pF_
				end
			end
		end
	end
	return pA,pWork,pF
end

function GetDynamicMoments(pA,pWork,pF,π0,year,a0,a1)
	# π0 is a NK x 17 probability distribution
	T = size(pA)[1]
	NK = size(pA)[2]
	TLlength = size(pA)[4]-1
	EA = zeros(Real,T)
	EH = zeros(Real,T)
	numerator = 0
	denom = 0
	for a = 0:16
		for nk=1:NK
			π1 = zeros(Real,T,TLlength+1)
			π1[1,1] = π0[nk,a+1]
			for t=1:T
				for w=1:TLlength+1
					EA[t] += π1[t,w]*pA[t,nk,a+1,w]
					EH[t] += π1[t,w]*(pWork[t,nk,a+1,w,1] + pA[t,nk,a+1,w]*(pWork[t,nk,a+1,w,2]-pWork[t,nk,a+1,w,1]))
					if t<T
						if w<=TLlength
							π1[t+1,w] += π1[t,w]*(1-pA[t,nk,a+1,w])
							π1[t+1,w+1] += π1[t,w]*pA[t,nk,a+1,w]
						else
							π1[t+1,w] += π1[t,w]
						end
					end
				end
			end
			if (a>=a0) & (a<=a1)
				for w=1:(TLlength+1)
					mass1 = pA[year,nk,1+a,w]*pWork[year,nk,1+a,w,2]
					mass0 = (1-pA[year,nk,1+a,w])*pWork[year,nk,1+a,w,1]
					numerator += π1[year,w]*(mass1*pF[year,nk,1+a,w,2]+mass0*pF[year,nk,1+a,w,1])
					denom += π1[year,w]*(mass0+mass1)
				end
			end
		end
	end
	return EA,EH,numerator/denom
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

function GetChildOutcomesStatic(year_meas,pA,pWork,pF,Y,pars_prod,pars)
	NK = size(pA)[2]
	TH = zeros(Real,NK,17)
	for nk=1:NK
		for a0 = 0:16
			th = 0
			for t=1:year_meas
				age = a0+t-1
				I = zeros(Real,2)
				if age<=17
					if age<=5
						gN,gF,δI = pars_prod.gN[1],pars_prod.gN[1],pars_prod.δI[1]
					else
						gN,gF,δI = pars_prod.gN[2],pars_prod.gN[2],pars_prod.δI[1]
					end
					for p=0:1
						If = log(Y[t,nk+1,1+p,2]+(112-30)*pars.wq-price) - gF
						In = log(Y[t,nk+1,1+p,2]+(112-30)*pars.wq) - gN
						Im = log(Y[t,nk+1,1+p,1] + 112*pars.wq)
						I[p+1] += Im + pWork[t,nk,age+1,p+1]*(In + pF[t,nk,age+1,p+1]*(If-In)-Im)
					end
					th = δI*(I[1] + pA[t,nk,age+1]*(I[2]-I[1])) + pars_prod.δθ*th
				end
			end
			TH[nk,a0] = th
		end
	end
	return TH
end

# NEXT: replace with wq and price (then don't need)
function GetChildOutcomesDynamic(year_meas,pA,pWork,pF,Y,Y_I,pars_prod,pars)
	NK = size(pA)[2]
	TH = zeros(Real,NK,17)
	TLlength = size(pA)[4]-1
	for nk=1:NK
		for a0 = 0:16
			th = 0
			π1 = zeros(Real,year_meas+1,TLlength+1)
			π1[1,1] = 1.
			for t=1:year_meas
				age = a0+t-1
				Imean = 0
				if age<=17
					if age<=5
						gN,gF,δI = pars_prod.gN[1],pars_prod.gN[1],pars_prod.δI[1]
					else
						gN,gF,δI = pars_prod.gN[2],pars_prod.gN[2],pars_prod.δI[2]
					end
					# case: time limit not reached yet
					for w=1:TLlength
						I = zeros(Real,2)
						for p=0:1
							If = log(Y[t,nk+1,1+p,2]+(112-30)*pars.wq-price) - gF
							In = log(Y[t,nk+1,1+p,2]+(112-30)*pars.wq) - gN
							Im = log(Y[t,nk+1,1+p,1] + 112*pars.wq)
							I[p+1] += Im + pWork[t,nk,age+1,w,p+1]*(In + pF[t,nk,age+1,w,p+1]*(If-In)-Im)
						end
						Imean += π1[t,w]*(I[1]+ pA[t,nk,w,age+1]*(I[2]-I[1]))
						π1[t+1,w] += π1[t,w]*(1-pA[t,nk,age+1,w])
						π1[t+1,w+1] += π1[t,w]*pA[t,nk,age+1,w]
					end
					# case: time limit reached, update
					w_ = TLlength+1
					I = zeros(Real,2)
					for p=0:1
						If = log(Y_I[t,nk+1,1+p,2]+(112-30)*pars.wq-price) - gF
						In = log(Y_I[t,nk+1,1+p,2]+(112-30)*pars.wq) - gN
						Im = log(Y_I[t,nk+1,1+p,1] + 112*pars.wq)
						I[p+1] += Im + pWork[t,nk,age+1,w_,p+1]*(In + pF[t,nk,age+1,w_,p+1]*(If-In)-Im)
					end
					Imean += π1[t,w_]*(I[1]+ pA[t,nk,w_,age+1]*(I[2]-I[1]))
					π1[t+1,w_] += π1[t,w_]
					th = δI*Imean + pars_prod.δθ*th
				end
			end
			TH[nk,a0] = th
		end
	end

	return TH
end
