using Random
using Distributions
# what structure here? Maybe make so we can simulate across sites??

# NS - number of sites
# NT - number of treatments (3)
mutable struct Model
	# preferences
	αc::Float64
	αθ::Float64
	αH::Array{Float64,1} #<- one per site
	αA::Array{Float64,1} #<- one per site
	β::Float64 #<- discounting

	# production
	δI::Array{Float64,1} #<- one per age of child
	δθ::Array{Float64,1} #<- one per age of child
	ϵ::Array{Float64,1} #<- one per age of child

	# prices
	pc::Array{Float64,1} #<- one per age of child (relative price)
	earnings::Array{Float64,2} #<- number of sites x number of quarters
	wq::Float64 #<- one for now, may need more!

	# policies (maybe don't need to store these all here? Only matter for budget))
	αWR::Float64 #<- cost imposed by work requirement
	TL::Array{Bool,2} #<- time limit indicator, NS x NT
	TLmax::Array{Int64,2} #<- length of time limit NS x NT
	B::Array{Float64,4} #<- max benefit(Quarter,NumChild,Site,Arm) (Q x 3 x NS x NT)
	D::Array{Float64,3} #<- Dollar disregard(NumChild,Site,Arm) (3 x NS x NT)
	R::Array{Float64,2} #<- % disregard(Site,Arm) NS x NT
	Bf::Array{Float64,4} #<- max food stamp benefit(Quarter,NumChild) (Q x 3)
	τ::Array{Float64,2} #<- NS x NT, childcare subsidy.

	# primitives
	πk::Array{Float64,2} #<- distribution of age youngest NS x 3 (<=2,3-5,6+)
	πK::Array{Float64,2} #<- distribution of number of children NS x 3 (1,2,3+)

	# endogenous objects
	Γδ::Array{Float64,1} #<- recursive coefficient for investment value
	budget::Array{Float64,7} # budget(site,arm,NumChild,age0,quarter,program,work) (NS x NT x 3 x 16 x Q x 2 x 2)
	utility::Array{Float64,7} # utility(site,arm,NumChild,age0,quarter,program,work) (NS x NT x 3 x 16 x Q x 2 x 2)
	utility_welf::Array{Float64,6} # = log(exp(u[p,1]) + exp(u[p,2])) for p=1,2
	V::Array{Float64,6} # value(site,arm,NumChild,age0,quarter,usage) (NS x NT x 3 x 16 x Q x TLmax+1)
	welf_prob::Array{Float64,6} # (site,arm,NumChild,age0,quarter,usage)
	work_prob::Array{Float64,6} # (site,arm,NumChild,age0,quarter,welfare_choice) (NS x NT x 3 x 16 x Q x 2)

end

function GetRecursiveCoefficient!(M)
	Q = 17*4+1
	M.Γδ[Q] = 1/(1-M.β)
	for q=Q-1:-1:1
		M.Γδ[q] = 1 + M.β*M.δθ[q]*M.Γδ[q+1]
	end
end

# this is already a bit fussy, I think we need to code up an individual file
function CalculateBudget!(M)
	Q = 17*4+1
	for s=1:3,tr=1:3,nk=1:3,age0=1:16,q=1:Q
		qT = (17-age0)*4
		earn = M.earnings[s,q]
		if M.TL[s,tr] & (q<=qT)
			for wu = 1:M.TLmax[s,tr]
				B0 = M.Bf[q,nk] + M.B[s,tr,nk,q]
				netinc_fs = max(0.8*earn - 134,0)
				netinc_a = max(earn - M.D[s,tr,nk,q])
				M.budget[s,tr,nk,age0,q,wu,1,1] = 0
				M.budget[s,tr,nk,age0,q,wu,1,2] = earn #<- no taxes yet
				M.budget[s,tr,nk,age0,q,wu,2,1] = B0
				M.budget[s,tr,nk,age0,q,wu,2,1] = B0 - 0.3*netinc_fs - (1-M.R[s,tr])*netinc_a
				for p=1:2

				end
			end
			# special case here for wu=TLmax+1
		else #<- if time limits don't apply, don't need to solve over dimension of welfare use
			for p=1:2
				u0 = M.utility[s,tr,nk,age0,q,1,p,1]
				u1 = M.utility[s,tr,nk,age0,q,1,p,2]
				M.work_prob[s,tr,nk,age0,q,1,p] = 1/exp(u0-u1)
				M.utility_welf[s,tr,nk,age0,q,1,p] = log(exp(u0)+exp(u1))
			end
		end
	end
end

#<- thought: only impact of # kids is on benefits, but should we maybe include something in alphaV?
function CalculateUtilities!(M)
	# first get recursive coefficient
	Q = 17*4+1
	for s=1:3,tr=1:3,nk=1:3,age0=1:16,q=1:Q
		age = age0*4+q
		if age<Q
			αV = M.αθ*M.δI[age]*M.β*M.Γδ[age+1]
			if M.TL[s,tr]
				for wu=1:M.TLmax[s,tr]
					for p=1:2,h=1:2
						Y = M.budget[s,tr,nk,age0,q,wu,p,h]
						hr = (h-1)*30
						U = (M.αC+αV)*log(Y+M.wq*(112-hr)) - αV/(1-M.ϵ[age])*log(112-hrs+hrs*M.pc[age]^(1-M.ϵ[age])) - M.αH[s]*(h-1) - M.αA[s]*(p-1) - M.αWR*(p-1)*(2-h)
						M.utility[s,tr,nk,age0,q,wu,p,h] = U
					end
					# if reach max time limit, only get food stamps
				end
				wu = M.TLmax[s,tr]+1
				for p=1:2,h=1:2
					Y = M.budget[s,tr,nk,age0,q,wu,p,h]
					hr = (h-1)*30
					U = (M.αC+αV)*log(Y+M.wq*(112-hr)) - αV/(1-M.ϵ[age])*log(112-hrs+hrs*M.pc[age]^(1-M.ϵ[age])) - M.αH[s]*(h-1) - M.αA[s]*(p-1)
					M.utility[s,tr,nk,age0,q,wu,p,h] = U
				end
			else
				for p=1:2,h=1:2
					Y = M.budget[s,tr,nk,age0,q,1,p,h]
					hr = (h-1)*30
					U = (M.αC+αV)*log(Y+M.wq*(112-hr)) - αV/(1-M.ϵ[age])*log(112-hrs+hrs*M.pc[age]^(1-M.ϵ[age])) - M.αH[s]*(h-1) - M.αA[s]*(p-1) - M.αWR*(p-1)*(2-h)
					M.utility[s,tr,nk,age0,q,1,p,h] = U
				end
			end
		else #<- no investment stuff
			for p=1:2,h=1:2
				Y = M.budget[s,tr,nk,age0,q,1,p,h]
				hr = (h-1)*30
				U = M.αC*log(Y+M.wq*(112-hr)) - M.αH[s]*(h-1) - M.αA[s]*(p-1)
				M.utility[s,tr,nk,age0,q,1,p,h] = U
			end
		end
	end
end

# this function takes utility calculations and calculates work choice probabilities, as well as expected utilities for program participation choices
function SolveWorkProb!(M::Model)
	for s=1:3,tr=1:3,nk=1:3,age0=1:16,q
		qT = (17-age0)*4
		if M.TL[s,tr] & (q<=qT)
			for wu = 1:M.TLmax[s,tr]+1
				for p=1:2
					u0 = M.utility[s,tr,nk,age0,q,wu,p,1]
					u1 = M.utility[s,tr,nk,age0,q,wu,p,2]
					M.work_prob[s,tr,nk,age0,q,wu,p] = 1/exp(u0-u1)
					M.utility_welf[s,tr,nk,age0,q,wu,p] = log(exp(u0)+exp(u1))
				end
			end
		else #<- if time limits don't apply, don't need to solve over dimension of welfare use
			for p=1:2
				u0 = M.utility[s,tr,nk,age0,q,1,p,1]
				u1 = M.utility[s,tr,nk,age0,q,1,p,2]
				M.work_prob[s,tr,nk,age0,q,1,p] = 1/exp(u0-u1)
				M.utility_welf[s,tr,nk,age0,q,1,p] = log(exp(u0)+exp(u1))
			end
		end
	end
end


# let's assume we don't need values further than for probabilities
# this model solves for welfare choice probabilities
function SolveModel!(M::Model,s,tr,nk,age0)
	Q = 17*4 + 1
	# first solve terminal period values (improve this for choice probs, etc)
	#M.V[s,tr,nk,age0,Q,:] .= log.(sum(exp.(M.utility[s,tr,nk,age0,Q,:,:])))
	# terminal period for eligibility:
	qT = (17-age0)*4
	for q=Q-1:-1:1
		age = age0 + floor(q-1,4)
		if M.TL[s,tr] & (q<=qT)
			for wu = 1:M.TLmax[s,tr]+1
				# an efficiency gain here, potentially
				v0 = M.utility_welf[s,tr,nk,age0,q,wu,1] + M.β*M.V[s,tr,nk,age0,q+1,wu]
				v1 = M.utility_welf[s,tr,nk,age0,q,wu,2] + M.β*M.V[s,tr,nk,age0,q+1,max(wu+1,M.TLmax[s,tr]+1)]
				M.V[s,tr,nk,age0,q,wu] = log(exp(v0)+exp(v1))
				M.welf_prob[s,tr,nk,age0,q,wu] = 1/(1+exp(v0-v1))
			end
		else
			# no time limits
			# an efficiency gain here, potentially
			v0 = M.utility_welf[s,tr,nk,age0,q,1,1] + M.β*M.V[s,tr,nk,age0,q+1,1]
			v1 = M.utility_welf[s,tr,nk,age0,q,1,2] + M.β*M.V[s,tr,nk,age0,q+1,1]
			M.welf_prob[s,tr,nk,age0,q,1] = 1/(1+exp(v0-v1))
		end
	end
end



# initial variables: age of youngest, # of kids, wage

# parameters needed:
# - δI,δθ,wq,pc*τ,ϵ

# if we want the government's expenditure, it's Xc*τ?
function Simulate(M::Model,Q,s,tr,nk,age0)
	A = zeros(Q) #<- participation
	A2 = zeros(Q) #<- receipt
	Y = zeros(Q) #<- earnings
	L = zeros(Q) #<- work
	θ = zeros(Q+1)
	Xc = zeros(Q)
	age = age0*4 #<- convert age to quarterly number
	wu = 0
	for q=1:Q
		welf = rand()<M.welf_prob[s,tr,nk,age0,q,1+wu] #<- no heterogeneity here
		L[q] = rand()<M.work_prob[s,tr,nk,age0,q,1+welf]
		Y[q] = L[q]*M.earnings[s,q] #<-
		if welf
			payment = M.budget[s,tr,nk,age0,q,2,1+L[q]]-M.earnings[s,q]
			A2[q] = payment
		end
		inc = Y[q] + A2[q]
		h = 30*L[q]
		if age<=17*4
			pc = M.pc[1+age]*(1-M.τ[s,tr])
			ϵ = M.ϵ[1+age]
			Xc[q] = h*pc^(1-ϵ)/(112 -h + h*pc^(1-ϵ))*inc
			θq[q+1] = M.δI[1+age]*log(inc + M.wq*(112-h)) - 1/(1-ϵ)*log(112-h + h*pc^(1-ϵ)) + M.δθ[1+age]*θ[q]
		end
		age += 1
		if M.TL[s,tr]
			A[q] = (age<=17*4)*(A2[q]>0)*(wu<M.TLmax[s,tr]+1)
			wu = max(wu+welf,M.TLmax[s,tr]+1)
		else
			A[q] = (age<=17*4)*(A2[q]>0)
		end
	end
	return A,A2,Y,L,θ,Xc
end

# simulate R panels of length Q for a site x treatment combination
function Simulate(M::Model,R,Q,s,tr)
	# ay_cat = rand(Multinomial(M.πk[s,tr]))
	# cats = [0:2,3:5,6:16]
	# age0 =
	A = zeros(Q,R) #<- participation
	A2 = zeros(Q,R) #<- receipt
	Y = zeros(Q,R) #<- earnings
	L = zeros(Q,R) #<- work
	θ = zeros(Q+1,R)
	Xc = zeros(Q,R)
	AGE = zeros(Q,R) #<- convert age to quarterly number
	age_dist = Multinomial(1,M.πk[s,:])
	age_cats = [0:2,3:5,6:16]
	nk_dist = Multinomial(1,M.πK[s,:])
	for r=1:R
		ac = rand(age_dist)
		age0 = rand(age_cats[ac])
		nk = rand(nk_dist)
		A[:,r],A2[:,r],Y[:,r],L[:,r],θ[:,r],Xc[:,r] = Simulate(M,Q,s,tr,nk,age0)
		AGE[:,r] = age0 + floor((0:Q-1)/4)
	end
	return AGE,A,A2,Y,L,θ,Xc
end
