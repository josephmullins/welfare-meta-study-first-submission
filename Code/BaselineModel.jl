using Random
using Distributions
using Revise
using Profile
using CSV
using GLM
using DataFrames
cd("/Users/FilipB/github/welfare-meta-study/Code")
include("Budget_Function_Code.jl") # will use budget1, Earnings, Budget_Ageout

budget1
findmax(budget1[:,3,:,:,:,:,:,:])

Budget_Ageout
Earnings
TimeLimit_Ind
TimeLimits
τ=ones(3,3).*0.5
findmax(TimeLimits)[1]
1+1
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

	τ::Array{Float64,2} #<- NS x NT, childcare subsidy.

	# primitives
	πk::Array{Float64,2} #<- distribution of age youngest NS x 3 (<=2,3-5,6+)
	πK::Array{Float64,2} #<- distribution of number of children NS x 3 (1,2,3+)

	# endogenous objects
	Γδ::Array{Float64,1} #<- recursive coefficient for investment value
	budget::Array{Float64,8} #  budget(site,arm,NumChild,age0,quarter,eligible,program,work)
	budget_ageout::Array{Float64,4} # (site,  quarter, program, work)
	utility::Array{Float64,8} # utility(site,arm,NumChild,age0,quarter,program,work) (NS x NT x 3 x 16 x Q x 2 x 2)
	utility_welf::Array{Float64,7} # = log(exp(u[p,1]) + exp(u[p,2])) for p=1,2
	V::Array{Float64,6} # value(site,arm,NumChild,age0,quarter,usage) (NS x NT x 3 x 16 x Q x TLmax+1)
	welf_prob::Array{Float64,6} # (site,arm,NumChild,age0,quarter,usage)
	work_prob::Array{Float64,7} # (site,arm,NumChild,age0,quarter,welfare_choice) (NS x NT x 3 x 16 x Q x 2)

end
j=ones(6,6)*0.5





function Model(N_sites, N_arms, N_age, budget,budget_ageout,qs, Earn, TimeLimit_Ind,TimeLimits,τ; πk=ones(3,3).*0.33,πK=ones(3,3).*0.33)
	αc=0.5
	αθ=0.5
	αH=ones(3)*0.5 #<- one per site
	αA=ones(N_sites)*0.5#<- one per site
	β=0.95 # seems like an ok guess

	# production
	δI=ones(qs)*0.5 #<- one per age of child
	δθ=ones(qs)*0.5 #<- one per age of child, in quarters for compatibility later
	ϵ=ones(qs)*0.5 #<- one per age of child

	# prices
	pc=ones(qs)*0.5 #<- one per age of child (relative price)
	earnings=Earn #<- number of sites x number of quarters
	wq=50.0 #<- one for now, may need more!

	# policies (maybe don't need to store these all here? Only matter for budget))
	αWR=10.0#<- cost imposed by work requirement
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

	return Model(αc,αθ,αH,αA,β,δI,δθ,ϵ,pc, earnings,wq,αWR,TL,TLmax,τ,πk,πK,Γδ,budget,budget_ageout,utility,utility_welf,V,welf_prob,work_prob)
end



function GetRecursiveCoefficient!(M)
	Q = 17*4+1
	M.Γδ[Q] = 1/(1-M.β)
	for q=Q-1:-1:1
		M.Γδ[q] = 1 + M.β*M.δθ[q]*M.Γδ[q+1]
	end
end



	#<- thought: only impact of # kids is on benefits, but should we maybe include something in alphaV?
function CalculateUtilities!(M)
	# first get recursive coefficient
	Q = 17*4+1
	for s=1:3,tr=1:3,nk=1:3,age0=1:17,q=1:Q
		age = age0*4+q
		if age<Q
			αV = M.αθ*M.δI[age]*M.β*M.Γδ[age+1]
			if M.TL[s,tr]
				for wu=1:M.TLmax[s,tr]
					for p=1:2,h=1:2
						Y = M.budget[s,tr,nk,age0,q,2,p,h] # I change eligibility to be a 0-1 variable, 2 indexing eligible
						hr = (h-1)*30 # the old draft had hrs below not hr--I assume it was a typo?
						U = (M.αc+αV)*log(Y+M.wq*(112-hr)) - αV/(1-M.ϵ[age])*log(112-hr+hr*M.pc[age]^(1-M.ϵ[age])) - M.αH[s]*(h-1) - M.αA[s]*(p-1) - M.αWR*(p-1)*(2-h)
						M.utility[s,tr,nk,age0,q,wu,p,h] = U
					end
					# if reach max time limit, only get food stamps
				end
				wu = M.TLmax[s,tr]+1
				for p=1:2,h=1:2
					Y = M.budget[s,tr,nk,age0,q,1,p,h] # if past time limit, then ineligible
					hr = (h-1)*30 # the old draft had hrs below not hr--I assume it was a typo?
					U = (M.αc+αV)*log(Y+M.wq*(112-hr)) - αV/(1-M.ϵ[age])*log(112-hr+hr*M.pc[age]^(1-M.ϵ[age])) - M.αH[s]*(h-1) - M.αA[s]*(p-1)
					M.utility[s,tr,nk,age0,q,wu,p,h] = U
				end
			else
				for p=1:2,h=1:2
					Y = M.budget[s,tr,nk,age0,q,1,p,h]
					hr = (h-1)*30 # the old draft had hrs below not hr--I assume it was a typo?
					U = (M.αc+αV)*log(Y+M.wq*(112-hr)) - αV/(1-M.ϵ[age])*log(112-hr+hr*M.pc[age]^(1-M.ϵ[age])) - M.αH[s]*(h-1) - M.αA[s]*(p-1) - M.αWR*(p-1)*(2-h)
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
	for s=1:3,tr=1:3,nk=1:3,age0=1:17,q=1:(length(M.δI)+1)
		qT = (17-age0)*4
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
	Q = 17*4 + 1
	# first solve terminal period values (improve this for choice probs, etc)
	#M.V[s,tr,nk,age0,Q,:] .= log.(sum(exp.(M.utility[s,tr,nk,age0,Q,:,:])))
	# terminal period for eligibility:
	qT = (17-age0)*4
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
	Mod1=Model(3,3,17,budget1,Budget_Ageout,17*4,Earnings,TimeLimit_Ind,TimeLimits,τ)
	GetRecursiveCoefficient!(Mod1)
	CalculateUtilities!(Mod1)
	SolveWorkProb!(Mod1)
	for s in 1:3
		for tr in 1:3
			for nk in 1:3
				for age0 in 1:17
					SolveModel!(Mod1,s,tr,nk,age0)
				end
			end
		end
	end
	return Mod1
end

@time Mod1=initialize_model()



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
Mod1.αc
UpdateSpecificParams!(Mod1; αc=0.3)
Mod1.αc
UpdateSpecificParams!(Mod1)
Mod1.αc
1+1


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
	if age==0
		age=1 # I don't want to worry about 1-indexing
	end
	wu = 1 # changed it to 1 here
	for q=1:Q

		welf = rand()<M.welf_prob[s,tr,nk,age0+1,q,wu] #<- no heterogeneity here
		L[q] = rand()<M.work_prob[s,tr,nk,age0+1,q,wu,1+welf] #pretty sure this needs another dimension
		Y[q] = L[q]*M.earnings[s,q] #<-

		if welf
			payment=0
			if age<=17*4
				payment = M.budget[s,tr,nk,age0+1,q,2,2,1+convert(Int,L[q])]-(M.earnings[s,q]*L[q])
				# only subtract earnings from budget if working...I think???
			else
				payment = M.budget_ageout[s,q,2,1+convert(Int,L[q])]-(M.earnings[s,q]*L[q])


			end

			A2[q] = payment
		end

		inc = Y[q] + A2[q]
		h = 30*L[q]
		if age<=17*4
			pc = M.pc[age]*(1-M.τ[s,tr])
			ϵ = M.ϵ[age]
			Xc[q] = h*pc^(1-ϵ)/(112 -h + h*pc^(1-ϵ))*inc
			θ[q+1] = M.δI[age]*log(inc + M.wq*(112-h)) - 1/(1-ϵ)* M.δI[age]*log(112-h + h*pc^(1-ϵ)) + M.δθ[age]*θ[q]

		end



		age += 1
		if M.TL[s,tr]
			A[q] = (age<=17*4)*(A2[q]>0)*(wu<M.TLmax[s,tr]+1)
			wu = min(wu+welf,M.TLmax[s,tr]+1) # I think this should be min not max?
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
		ac2=findmax(ac)[2]
		age0 = rand(age_cats[ac2])
		nk = rand(nk_dist)
		nk2=findmax(nk)[2]
		A[:,r],A2[:,r],Y[:,r],L[:,r],θ[:,r],Xc[:,r] = Simulate(M,Q,s,tr,nk2,age0)
		for i in 1:Q
		AGE[i,r] = age0 + floor((i)/4) # ages start at 0, had issue with 0-indexinf
		end
	end
	return AGE,A,A2,Y,L,θ,Xc
end





@time Mod2=initialize_model()








#=


Below are a few sanity checks

=#



@time S1_c=Simulate(Mod2,10000,30,1,1)
@time S2_c=Simulate(Mod2,10000,30,2,1)
@time S3_c=Simulate(Mod2,10000,30,3,1)


mean(S1_c[2])
mean(S2_c[2])
mean(S3_c[2])


mean(S1_c[3])
mean(S2_c[3])
mean(S3_c[3])

@time S1_t=Simulate(Mod2,100,30,1,3)
@time S2_t=Simulate(Mod2,100,30,2,3)
@time S3_t=Simulate(Mod2,100,30,3,3)

budget1[3,3,1,3,25,2,2,2]
Earnings[3,25]

mean(S1_t[2])
mean(S2_t[2])
mean(S3_t[2])


Mod2.welf_prob[1,1,:,:,:,:]

Mod2.welf_prob[1,2,:,:,:,:]

mean(S1_t[3])
mean(S2_t[3])
mean(S3_t[3])

mean(S1_c[4])
mean(S1_t[4])

mean(S1_c[5])
mean(S1_t[5])

mean(S1_c[6])
mean(S1_t[6])

mean(S1_c[7])
mean(S1_t[7])
findmax(budget1[:,3,:,:,:,:,:,:])

UpdateSpecificParams!(Mod2; pc=ones(length(Mod2.pc))*0.1,ϵ=ones(length(Mod2.ϵ))*0.1)

@time S1_t=Simulate(Mod2,100,30,1,3)

mean(S1_t[6])

S1_t[6]
