# goal: get nsims, model, set of pars; feed into function that generates moment fit
using NLopt

mutable struct Parameters
	np::NamedTuple{(:αc, :αθ, :αH, :αA, :β, :δI, :δθ, :ϵ, :τ, :pc, :wq, :αWR),NTuple{12,Int64}}
	lb::NamedTuple{(:αc, :αθ, :αH, :αA, :β, :δI, :δθ, :ϵ, :τ, :pc, :wq, :αWR),Tuple{Float64, Float64, Array{Float64,1}, Array{Float64,1}, Float64, Array{Float64,1},Float64,Float64,Array{Float64,1},Array{Float64,1},Float64,Float64 }}
	ub::NamedTuple{(:αc, :αθ, :αH, :αA, :β, :δI, :δθ, :ϵ, :τ, :pc, :wq, :αWR),Tuple{Float64, Float64, Array{Float64,1}, Array{Float64,1}, Float64, Array{Float64,1},Float64,Float64,Array{Float64,1},Array{Float64,1},Float64,Float64 }}


    # preferences
    αc::Float64
    αθ::Float64
    αH::Array{Float64,1} #<- one per site, 4
    αA::Array{Float64,1} #<- one per site,4
    β::Float64 #<- discounting

    # production
    δI::Array{Float64,1} #<- one per age of child in years, 2 params
    δθ::Float64 #<- in theory one per age of child in years, in practice 1 number
    ϵ::Float64 #<- one per age of child in years in theory, 1 in practice


    τ::Array{Float64,1} #<- NS x NT, childcare subsidy, 2 in practice
		# will only update 2 (1 treatment, 1 for both controls for both arms in MN)

    pc::Array{Float64,1} #<- one per age of child (relative price), 2 params

    wq::Float64 #<- one for now, may need more! 1

	# policies (maybe don't need to store these all here? Only matter for budget))
	αWR::Float64 # only some arms have work reqs 1
end


function UpdateModel!(M::Model, Pars::Parameters)

	M.αc=Pars.αc
	M.αθ=Pars.αθ
	M.αH=Pars.αH
	M.β=Pars.β
	for i in 1:length(M.δI)
		M.δI[i]=exp(Pars.δI[1]+i*Pars.δI[2])
		M.δθ[i]=Pars.δθ
		M.ϵ[i]=Pars.ϵ
	end
	for i in 3:4
		M.τ[i,1]=Pars.τ[1]
		for j in 2:3
			M.τ[i,j]=Pars.τ[2]
		end
	end
	for i in 1:length(M.pc)
		M.pc[i]=exp(Pars.pc[1]+i*Pars.pc[2])

	end

	M.wq=Pars.wq
	M.αWR=Pars.αWR.*M.work_reqs
end



function UpdatePars!(x::Array{Float64,1},
	Pars::Parameters,
	vars::Array{Symbol,1})

	pos = 1
	for s in vars
		np = Pars.np[s]
		if np>1
			xrange = pos:(pos+np-1)
			setfield!(Pars,s,x[xrange])
		else
			setfield!(Pars,s,x[pos])
		end
		pos += np
	end
end



# this function assumed a dataframe. This will still be useful but let's work on it later.
# function CheckFit!(M::Model,moms,wghts,R,lengths,TE_index)
# 		SolveModel!(M)
# 		mom_sim=BaselineMoments(M,R,lengths,TE_index)
# 		moms[:mom_sim] = mom_sim
# 		moms[:Q] = moms[:wght].*(moms[:mom] - moms[:mom_sim]).^2
# end
#
# function CheckFit!(Pars::Parameters, M::Model, moms::DataFrame)
# 	UpdateModel!(M, Pars)
# 	SolveModel!(M)
# 	CheckFit!(moms,M Nsims=10000)
# end


function Criterion(x::Array{Float64,1},
		Pars::Parameters,
		M::Model,
		vars::Array{Symbol,1},moms,wghts,R,lengths,TE_index)

	UpdatePars!(x,Pars,vars)
	UpdateModel!(M, Pars)
	SolveModel!(M)
	E,A,A2,XG,skill = MomentsBaseline(M,R,lengths,TE_index)
	mom_sim = [E;A;A2;XG;skill]
	Q = sum(wghts.* (mom_sim .- moms).^2)
	println(Q)
	return Q
end


# another version with an approximated gradient
# might be useful for different optimizers
function Criterion!(x,g,Pars,M,moms,vars)
	np = length(x)
	Δ = 0.01
	F0=Criterion(x,Pars,M,moms, vars)
	for i in 1:np
		x1 = copy(x)
		x1[i] += Δ
		F1 = Criterion(x,Pars,M,moms, vars)
		g[i] = (F1-F0)/Δ
	end

end

function GetOptimization(Pars::Parameters,M::Model,vars,moms,wghts,R,lengths,TE_index;
	Global=0,
	SBPLX=0,
	maxevals=1000,
	tightness=1e-3)

	UpdateModel!(M, Pars)
	SolveModel!(M)
	np = 0;
	for s in vars
		np += Pars.np[s]
	end
	x0 = zeros(np)
	lb = zeros(np)
	ub = zeros(np)
	if Global==1
		opt = Opt(:G_MLSL_LDS,np)
		#opt = Opt(:GN_CRS2_LM,np)

		if SBPLX==1
			local_opt = Opt(:LN_SBPLX,np)
		else
			local_opt = Opt(:LN_NELDERMEAD,np)
		end
		ftol_abs!(local_opt,1e-1) #<-!?
		xtol_abs!(local_opt,1e-1)
		opt.local_optimizer = local_opt::Opt
		opt.population=50
	elseif Global==0 && SBPLX==1
		opt = Opt(:LN_SBPLX,np)
	else
		opt = Opt(:LN_NELDERMEAD,np)
	end
	pos = 1
	for s in vars
		np = Pars.np[s]
		if np>1
			xrange = pos:(pos+np-1)
			x0[xrange] = getfield(Pars,s)
			lb[xrange] = pars.lb[s]
			ub[xrange] = pars.ub[s]
		else
			x0[pos] = getfield(Pars,s)
			lb[pos] = pars.lb[s]
			ub[pos] = pars.ub[s]
		end
		pos += np
	end
	ftol_abs!(opt,tightness) #<-!?
	xtol_abs!(opt,tightness)
	lower_bounds!(opt,lb)
	upper_bounds!(opt,ub)
		# below it was (x,g) in Jo's old code
		# don't know what the g was for--test this and see if it breaks (yes it breaks)
	min_objective!(opt,(x,g)->Criterion(x,Pars,M,vars,moms,wghts,R,lengths,TE_index))
	maxeval!(opt,maxevals)

	return opt,x0

end
