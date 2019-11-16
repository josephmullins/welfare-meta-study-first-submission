using Random
using Distributions
using Revise
using Profile
using CSV
using GLM
using DataFrames
#cd("/Users/FilipB/github/welfare-meta-study/Code")
#include("Budget_Function_Code.jl") # will use budget1, Earnings, Budget_Ageout

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
	δI::Array{Float64,1} #<- one per age of child in years
	δθ::Array{Float64,1} #<- one per age of child in years
	ϵ::Array{Float64,1} #<- one per age of child in years

	# prices
	pc::Array{Float64,1} #<- one per age of child (relative price)
	earnings::Array{Float64,2} #<- number of sites x number of quarters
	SNAP::Array{Float64,8}
	wq::Float64 #<- one for now, may need more!

	# policies (maybe don't need to store these all here? Only matter for budget))
	αWR_p::Float64 #<- cost imposed by work requirement
	αWR::Array{Float64,2} # only some arms have work reqs
	work_reqs::Array{Float64,2}
	TL::Array{Bool,2} #<- time limit indicator, NS x NT
	TLmax::Array{Int64,2} #<- length of time limit NS x NT

	τ::Array{Float64,2} #<- NS x NT, childcare subsidy.

	# primitives
	πk::Array{Float64,2} #<- distribution of age youngest NS x 3 (<=2,3-5,6+)
	πK::Array{Float64,2} #<- distribution of number of children NS x 3 (1,2,3+)

	# endogenous objects
	Γδ::Array{Float64,1} #<- recursive coefficient for investment value
	budget::Array{Float64,8} #  budget(site,arm,NumChild,age0,quarter,eligible,program,work)
	budget_ageout::Array{Float64,4} # (site,  quarter, program, work)
	utility::Array{Float64,8} # utility(site,arm,NumChild,age0,quarter,program,work) (NS x NT x 3 x 18 x Q x 2 x 2)
	utility_welf::Array{Float64,7} # = log(exp(u[p,1]) + exp(u[p,2])) for p=1,2
	V::Array{Float64,6} # value(site,arm,NumChild,age0,quarter,usage) (NS x NT x 3 x 18 x Q x TLmax+1)
	welf_prob::Array{Float64,6} # (site,arm,NumChild,age0,quarter,usage)
	work_prob::Array{Float64,7} # (site,arm,NumChild,age0,quarter,welfare_choice) (NS x NT x 3 x 18 x Q x 2)

end


function Model(N_sites, N_arms, N_age, budget,budget_ageout,qs, Earn, TimeLimit_Ind,TimeLimits,τ, Work_Reqs, Foodstamps_receipt; πk=ones(N_sites,N_arms).*(1/3),πK=ones(N_sites,N_arms).*(1/3))
	αc=0.5
	αθ=0.5
	αH=ones(N_sites)*0.5 #<- one per site
	αA=ones(N_sites)*0.5#<- one per site
	β=0.95 # seems like an ok guess

	# production
	δI=ones(N_age)*0.5 #<- one per age of child in years
	δθ=ones(N_age)*0.5 #<- one per age of child in years
	ϵ=ones(N_age)*0.5 #<- one per age of child in years

	# prices
	pc=ones(N_age)*0.5 #<- one per age of child (relative price)
	earnings=Earn #<- number of sites x number of quarters
	wq=50.0 #<- one for now, may need more!

	# policies (maybe don't need to store these all here? Only matter for budget))
	αWR_p=10.0#<- cost imposed by work requirement
	αWR=Work_Reqs.*αWR_p
	work_reqs=Work_Reqs
	TL=TimeLimit_Ind #<- time limit indicator, NS x NT
	TLmax=TimeLimits #<- length of time limit NS x NT
		htl=findmax(TLmax)[1]+1

	τ=τ#::Array{Float64,2} #<- NS x NT, childcare subsidy.

	# primitives
	πk=πk#::Array{Float64,2} #<- distribution of age youngest NS x 3 (<=2,3-5,6+)
	πK=πK#::Array{Float64,2} #<- distribution of number of children NS x 3 (1,2,3+)

	# endogenous objects
	Γδ=zeros(qs+1) #<- recursive coefficient for investment value, number of possible quarters with kid plus 1
	budget=budget # budget(site,arm,NumChild,age0,quarter,eligible,program,work)
	budget_ageout=budget_ageout # (site, arm, quarter, program, work)
	utility=zeros(N_sites,N_arms,3,N_age,qs+1,htl,2,2)#::Array{Float64,8} # utility(site,arm,NumChild,age0,quarter,usage,program,work) (NS x NT x 3 x 16 x Q x 2 x 2)
	utility_welf=zeros(N_sites,N_arms,3,N_age,qs+1,htl,2)#::Array{Float64,7} # = log(exp(u[p,1]) + exp(u[p,2])) for p=1,2
	V=zeros(N_sites,N_arms,3,N_age,qs+1,htl+1)#::Array{Float64,6} # value(site,arm,NumChild,age0,quarter,usage) (NS x NT x 3 x 16 x Q x TLmax+1)
	welf_prob=zeros(N_sites,N_arms,3,N_age,qs+1,htl)#::Array{Float64,6} # (site,arm,NumChild,age0,quarter,usage)
	work_prob=zeros(N_sites,N_arms,3,N_age,qs+1,htl,2)#::Array{Float64,7} # (site,arm,NumChild,age0,quarter,usage,welfare_choice) (NS x NT x 3 x 16 x Q x 2)

	return Model(αc,αθ,αH,αA,β,δI,δθ,ϵ,pc, earnings,Foodstamps_receipt,wq,αWR_p,αWR,Work_Reqs,TL,TLmax,τ,πk,πK,Γδ,budget,budget_ageout,utility,utility_welf,V,welf_prob,work_prob)
end



function GetRecursiveCoefficient!(M)
	Q = 18*4+1
	M.Γδ[18+1] = 1/(1-M.β)
	for q=Q-1:-1:1
		Age_Year=convert(Int,ceil(q/4))
		M.Γδ[q] = 1 + M.β*M.δθ[Age_Year]*M.Γδ[Age_Year+1]
	end
end


	#<- thought: only impact of # kids is on benefits, but should we maybe include something in alphaV?
function CalculateUtilities!(M)
	# first get recursive coefficient
	Q = 18*4+1
	for s=1:length(M.αH),tr=1:3,nk=1:3,age0=1:18,q=1:Q
		age = age0*4+q
		if age<Q
			Age_Year=convert(Int,ceil(q/4))
			αV = M.αθ*M.δI[Age_Year]*M.β*M.Γδ[Age_Year+1]
			if M.TL[s,tr]
				for wu=1:M.TLmax[s,tr]
					for p=1:2,h=1:2
						Y = M.budget[s,tr,nk,age0,q,2,p,h] # I change eligibility to be a 0-1 variable, 2 indexing eligible
						hr = (h-1)*30 # the old draft had hrs below not hr--I assume it was a typo?
						U = (M.αc+αV)*log(Y+M.wq*(112-hr)) - αV/(1-M.ϵ[Age_Year])*log(112-hr+hr*M.pc[Age_Year]^(1-M.ϵ[Age_Year])) - M.αH[s]*(h-1) - M.αA[s]*(p-1) - M.αWR[s,tr]*(p-1)*(2-h)
						M.utility[s,tr,nk,age0,q,wu,p,h] = U
					end
					# if reach max time limit, only get food stamps
				end
				wu = M.TLmax[s,tr]+1
				for p=1:2,h=1:2
					pc = M.pc[Age_Year]*(1-M.τ[s,tr])
					Y = M.budget[s,tr,nk,age0,q,1,p,h] # if past time limit, then ineligible
					hr = (h-1)*30 # the old draft had hrs below not hr--I assume it was a typo?
					U = (M.αc+αV)*log(Y+M.wq*(112-hr)) - αV/(1-M.ϵ[Age_Year])*log(112-hr+hr*pc^(1-M.ϵ[Age_Year])) - M.αH[s]*(h-1) - M.αA[s]*(p-1)
					M.utility[s,tr,nk,age0,q,wu,p,h] = U
				end
			else
				for p=1:2,h=1:2
					pc = M.pc[Age_Year]*(1-M.τ[s,tr])
					Y = M.budget[s,tr,nk,age0,q,2,p,h]
					hr = (h-1)*30 # the old draft had hrs below not hr--I assume it was a typo?
					U = (M.αc+αV)*log(Y+M.wq*(112-hr)) - αV/(1-M.ϵ[Age_Year])*log(112-hr+hr*pc^(1-M.ϵ[Age_Year])) - M.αH[s]*(h-1) - M.αA[s]*(p-1) - M.αWR[s,tr]*(p-1)*(2-h)
					M.utility[s,tr,nk,age0,q,1,p,h] = U
				end
			end
		else #<- no investment stuff
			for p=1:2,h=1:2
				Y = M.budget_ageout[s,q,p,h] # this is after the kids have aged out
				hr = (h-1)*30
				U = M.αc*log(Y+M.wq*(112-hr)) - M.αH[s]*(h-1) - M.αA[s]*(p-1)
				M.utility[s,tr,nk,age0,q,1,p,h] = U
			end
		end
	end
end
# utility(site,arm,NumChild,age0,quarter,program,work) (NS x NT x 3 x 16 x Q x 2 x 2)


# this function takes utility calculations and calculates work choice probabilities, as well as expected utilities for program participation choices
function SolveWorkProb!(M::Model)
	for s=1:4,tr=1:3,nk=1:3,age0=1:18,q=1:(length(M.δI)+1)
		qT = (18-age0)*4
		if M.TL[s,tr] & (q<=qT)
			for wu = 1:M.TLmax[s,tr]+1
				for p=1:2
					u0 = M.utility[s,tr,nk,age0,q,wu,p,1]
					u1 = M.utility[s,tr,nk,age0,q,wu,p,2]
					M.work_prob[s,tr,nk,age0,q,wu,p] = 1/(1+exp(u0-u1) )
					M.utility_welf[s,tr,nk,age0,q,wu,p] = log(exp(u0)+exp(u1))
				end
			end
		else #<- if time limits don't apply, don't need to solve over dimension of welfare use
			for p=1:2
				u0 = M.utility[s,tr,nk,age0,q,1,p,1]
				u1 = M.utility[s,tr,nk,age0,q,1,p,2]
				M.work_prob[s,tr,nk,age0,q,1,p] = 1/(1+exp(u0-u1) )# this didn't have the +1 in the denominator before
				M.utility_welf[s,tr,nk,age0,q,1,p] = log(exp(u0)+exp(u1))
			end
		end
	end
end


# let's assume we don't need values further than for probabilities
# this model solves for welfare choice probabilities
function SolveModel!(M::Model,s,tr,nk,age0)
	Q = 18*4 + 1
	# first solve terminal period values (improve this for choice probs, etc)
	#M.V[s,tr,nk,age0,Q,:] .= log.(sum(exp.(M.utility[s,tr,nk,age0,Q,:,:])))
	# terminal period for eligibility:
	qT = (18-age0)*4
	for q=Q-1:-1:1
		age = age0 + min(q-1,4)# this was floor--should it not be min as both are ints?
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




function initialize_model()
	Mod1=Model(4,3,18,budget1,Budget_Ageout,18*4,Earnings,TimeLimit_Ind,TimeLimits,τ, Work_Reqs_Ind,Foodstamps_receipt)
	GetRecursiveCoefficient!(Mod1)
	CalculateUtilities!(Mod1)
	SolveWorkProb!(Mod1)
	for s in 1:4
		for tr in 1:3
			for nk in 1:3
				for age0 in 1:18
					SolveModel!(Mod1,s,tr,nk,age0)
				end
			end
		end
	end
	return Mod1
end

function SolveModel!(Mod1::Model)
	GetRecursiveCoefficient!(Mod1)
	CalculateUtilities!(Mod1)
	SolveWorkProb!(Mod1)
	n_arms = [2,2,3,3]
	for s in 1:4
		for tr in 1:n_arms[s]
			for nk in 1:3
				for age0 in 1:18
					SolveModel!(Mod1,s,tr,nk,age0)
				end
			end
		end
	end
end


function UpdateSpecificParams!(M::Model;
	αc=M.αc, αθ=M.αθ, αH=M.αH, αA=M.αA,αWR=M.αWR,
	δI=M.δI,δθ=M.δθ,
	ϵ=M.ϵ, wq=M.wq,
	pc=M.pc)
	M.αc=αc
	M.αθ=αθ
	M.αH=αH
	M.αA=αA
	M.αWR=αWR
	M.δI=δI
	M.δθ=δθ
	M.ϵ=ϵ
	M.pc=pc
	M.wq=M.wq

end


# initial variables: age of youngest, # of kids, wage

# parameters needed:
# - δI,δθ,wq,pc*τ,ϵ

# if we want the government's expenditure, it's Xc*τ/(1-τ)
function Simulate(M::Model,Q,s,tr,nk,age0)
	A = zeros(Q) #<- participation
	A2 = zeros(Q) #<- receipt
	Y = zeros(Q) #<- earnings
	L = zeros(Q) #<- work
	θ = zeros(Q+1)
	Xc = zeros(Q)
	Foodstamps=zeros(Q)
	age = age0*4 #<- convert age to quarterly number
	if age==0
		age=1 # I don't want to worry about 1-indexing
	end
	wu = 1 # changed it to 1 here
	for q=1:Q
		Age_Year=convert(Int,ceil(q/4))
		welf = rand()<M.welf_prob[s,tr,nk,age0+1,q,wu] #<- no heterogeneity here
		L[q] = rand()<M.work_prob[s,tr,nk,age0+1,q,wu,1+welf] #pretty sure this needs another dimension
		Y[q] = L[q]*M.earnings[s,q] #<-
		Elig=1 # start off eligible
		if M.TL[s,tr] && wu==M.TLmax[s,tr]+1 # switch off only in places with time limits
			Elig=0
		end
		if welf
			payment=0
			if age<=18*4 && Elig==1
				payment = M.budget[s,tr,nk,age0+1,q,2,2,1+convert(Int,L[q])]-(M.earnings[s,q]*L[q])-
									M.SNAP[s,tr,nk,age0+1,q,2,2,1+convert(Int,L[q])]
				Foodstamps[q]=M.SNAP[s,tr,nk,age0+1,q,2,2,1+convert(Int,L[q])]
			elseif age<=18*4 && Elig==0 # switch index to 1 if no longer eligible
				payment = M.budget[s,tr,nk,age0+1,q,1,2,1+convert(Int,L[q])]-(M.earnings[s,q]*L[q])-
									M.SNAP[s,tr,nk,age0+1,q,1,2,1+convert(Int,L[q])]
				Foodstamps[q]=M.SNAP[s,tr,nk,age0+1,q,1,2,1+convert(Int,L[q])]
			else
				if s==1 || s==2 # count this as snap for MN, earnings otherwise
					payment = M.budget_ageout[s,q,2,1+convert(Int,L[q])]-(M.earnings[s,q]*L[q])
				else
					Foodstamps[q]=M.budget_ageout[s,q,2,1+convert(Int,L[q])]-(M.earnings[s,q]*L[q])
				end
			end

			A2[q] = payment
		end

		inc = Y[q] + A2[q]
		h = 30*L[q]
		if age<=18*4
			pc = M.pc[Age_Year]*(1-M.τ[s,tr])
			ϵ = M.ϵ[Age_Year]
			Xc[q] = h*pc^(1-ϵ)/(112 -h + h*pc^(1-ϵ))*inc
			θ[q+1] = M.δI[Age_Year]*log(inc + M.wq*(112-h)) - 1/(1-ϵ)* M.δI[Age_Year]*log(112-h + h*pc^(1-ϵ)) + M.δθ[Age_Year]*θ[q]
		else
			θ[q+1] = θ[q]
		end



		age += 1
		if M.TL[s,tr]
			A[q] = (age<=18*4)*(A2[q]>0)*(wu<M.TLmax[s,tr]+1)
			wu = min(wu+welf,M.TLmax[s,tr]+1) # I think this should be min not max?
		else
			A[q] = (age<=18*4)*(A2[q]>0)
		end

	end
	return A,A2,Y,L,θ,Xc, Foodstamps
end

# simulate R panels of length Q for a site x treatment combination
function Simulate(M::Model,R,Q,s,tr)
	# ay_cat = rand(Multinomial(M.πk[s,tr]))
	# cats = [0:2,3:5,6:16]
	# age0 =
	Random.seed!(123)
	A = zeros(Q,R) #<- participation
	A2 = zeros(Q,R) #<- receipt
	Y = zeros(Q,R) #<- earnings
	L = zeros(Q,R) #<- work
	θ = zeros(Q+1,R)
	Xc = zeros(Q,R)
	Foodstamps=zeros(Q,R)
	AGE = zeros(Q,R) #<- convert age to quarterly number
	age_dist = Multinomial(1,M.πk[s,:])
	age_cats = [0:2,3:5,6:16]
	nk_dist = Multinomial(1,M.πK[s,:])
	for r=1:R
		ac = rand(age_dist)
		ac2=findmax(ac)[2]
		age0 = rand(age_cats[ac2])
		nk = rand(nk_dist)
		nk2=findmax(nk)[2]
		A[:,r],A2[:,r],Y[:,r],L[:,r],θ[:,r],Xc[:,r],Foodstamps[:,R] = Simulate(M,Q,s,tr,nk2,age0)
		for i in 1:Q
			AGE[i,r] = age0 + floor((i-1)/4) # ages start at 0, had issue with 0-indexinf
		end
	end
	return (AGE=AGE,Participation=A,Benefit_Receipt=A2,
	Earned_Income=Y,LFP=L,Skills=θ,Childcare=Xc,Foodstamps=Foodstamps)
end

function MomentsBaseline(M::Model,R,lengths,TE_index)
	# forget about making general for now

	# -- set up some primitives
	sample_size = [4803,1405+1410,3208,6009]
	simsize = R*sample_size
	N1 = 2*(lengths[1] + lengths[2]) + 3*(lengths[3] + lengths[4])
	E = zeros(N1)
	A = zeros(N1)
	A2 = zeros(N1)
	n_arms = [2,2,3,3]
	years_childcare = [3,3,3,3]
	year_meas = [3,4,3,3]
	θ = -1*ones(4,3,maximum(simsize))
	AGE = -1*ones(4,3,maximum(simsize))
	XG = zeros(4,2)
	# --- Run the simulation ----- #

	curr_pos = 0
	for i=1:4
		for j=1:n_arms[i]
			dat = Simulate(M,simsize[i],lengths[i],i,j)
			curr_slice = (curr_pos+1):(curr_pos+lengths[i])
			E[curr_slice] = mean(dat.LFP,dims=2)
			A[curr_slice] = mean(dat.Participation,dims=2)
			A2[curr_slice] = mean(dat.Benefit_Receipt,dims=2)
			θ[i,j,1:simsize[i]] = dat.Skills[year_meas[i]*4,:]
			AGE[i,j,1:simsize[i]] = dat.AGE[year_meas[i]*4,:]
			if j<=2
				τ = Mod1.τ[i,j]
				XG[i,j] = τ/(1-τ)*mean(dat.Childcare[1:years_childcare[i]*4,:])*4 #<- annual
			end
			curr_pos += lengths[i]
		end
	end
	# --- Child Care Moments ------ #
	XG_moms = [(XG[i,2]-XG[i,1])/XG[i,1] for i=1:3]
	XG_mom2 = XG[3,1]/XG[1,1]

	# --- Skill Outcome Moments ------ #
	# assume TE_index stores: site,treatment,agemin,agemax
	skill_moms = zeros(size(TE_index)[1])
	for m=1:size(TE_index)[1]
		s,a,amin,amax = TE_index[m,:]
		slice_control = (AGE[s,1,:].>=amin) .& (AGE[s,1,:].<=amax)
		slice_treat = (AGE[s,1+a,:].>=amin) .& (AGE[s,1+a,:].<=amax)
		skill_moms[m] = mean(θ[s,1+a,slice_treat]) - mean(θ[s,1,slice_control])
	end
	return E,A,A2,[XG_moms;XG_mom2],skill_moms
end
