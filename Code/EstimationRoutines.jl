# goal: get nsims, model, set of pars; feed into function that generates moment fit
using NLopt

mutable struct Parameters
	np::NamedTuple{(:αc, :αθ, :αH, :αA, :β, :σA, :δI, :δθ, :ϵ, :τ, :pc, :wq, :αWR),NTuple{13,Int64}}
	lb::NamedTuple{(:αc, :αθ, :αH, :αA, :β, :σA, :δI, :δθ, :ϵ, :τ, :pc, :wq, :αWR),Tuple{Array{Float64,1}, Float64, Array{Float64,1}, Array{Float64,1}, Float64, Float64, Array{Float64,1},Float64,Float64,Array{Float64,1},Array{Float64,1},Float64,Float64 }}
	ub::NamedTuple{(:αc, :αθ, :αH, :αA, :β, :σA, :δI, :δθ, :ϵ, :τ, :pc, :wq, :αWR),Tuple{Array{Float64,1}, Float64, Array{Float64,1}, Array{Float64,1}, Float64, Float64, Array{Float64,1},Float64,Float64,Array{Float64,1},Array{Float64,1},Float64,Float64 }}


    # preferences
    αc::Array{Float64,1}
    αθ::Float64
    αH::Array{Float64,1} #<- one per site, 4
    αA::Array{Float64,1} #<- one per site,4
    β::Float64 #<- discounting
	σA::Float64

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
	M.αA = Pars.αA
	M.β=Pars.β
	M.σA = Pars.σA
	for i in 1:length(M.δI)
		M.δI[i]=exp(Pars.δI[1]+i*Pars.δI[2])
		M.δθ[i]=Pars.δθ
		M.ϵ[i]=Pars.ϵ
		M.pc[i] = Pars.pc[1]
	end
	for i in 3:4
		M.τ[i,1]=Pars.τ[1]
		for j in 2:3
			M.τ[i,j]=Pars.τ[2]
		end
	end
	# for i in 1:length(M.pc)
	# 	M.pc[i]=exp(Pars.pc[1]+i*Pars.pc[2])
	#
	# end

	M.wq=Pars.wq
	M.αWR=Pars.αWR.*M.work_reqs
end



function UpdatePars!(x::Array{Float64,1},
	Pars::Parameters,
	vars::Array{Symbol,1})

	pos = 1
	for s in vars
		np = Pars.np[s]
		if typeof(Pars.lb[s])==Float64
			setfield!(Pars,s,x[pos])
		else
			xrange = pos:(pos+np-1)
			setfield!(Pars,s,x[xrange])
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


function Criterion(x,Pars,M,vars,moms,wghts,R,lengths,TE_index,solve,verbose)
	UpdatePars!(x,Pars,vars)
	UpdateModel!(M, Pars)
	if solve
		SolveModel!(M)
	end
	Random.seed!(151119)
	E,A,A2,XG,skill,Y = MomentsBaseline(M,R,lengths,TE_index)
	mom_sim = [E;A;A2;XG;skill;Y]
	Q = sum(wghts.* (mom_sim .- moms).^2)
	if verbose
		println(Q)
	end
	return Q
end


# another version with an approximated gradient
# might be useful for different optimizers
function Criterion(x,g,Pars,M,vars,moms,wghts,R,lengths,TE_index,solve,verbose)
	np = length(x)
	Δ = 0.001
	F0=Criterion(x,Pars,M,vars,moms,wghts,R,lengths,TE_index,solve,verbose)
	for i in 1:np
		x1 = copy(x)
		x1[i] += Δ
		F1 = Criterion(x1,Pars,M,vars,moms,wghts,R,lengths,TE_index,solve,false)
		g[i] = (F1-F0)/Δ
	end
	return F0
end

function GetOptimization(Pars,M,vars,moms,wghts,R,lengths,TE_index;
	Global=0,
	SBPLX=0,
	LBFGS=0,
	maxevals=1000,
	tightness=1e-3,
	solve = true)

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
		elseif LBFGS==1
			local_opt = Opt(:LD_LBFGS,np)
		else
			local_opt = Opt(:LN_NELDERMEAD,np)
		end
		ftol_abs!(local_opt,1e-1) #<-!?
		xtol_abs!(local_opt,1e-1)
		opt.local_optimizer = local_opt::Opt
		opt.population=50
	elseif Global==0 && SBPLX==1
		opt = Opt(:LN_SBPLX,np)
	elseif Global==0 && LBFGS==1
		opt = Opt(:LD_LBFGS,np)
	else
		opt = Opt(:LN_NELDERMEAD,np)
	end
	pos = 1
	for s in vars
		np = Pars.np[s]
		if typeof(Pars.lb[s])==Float64
			x0[pos] = getfield(Pars,s)
			lb[pos] = pars.lb[s]
			ub[pos] = pars.ub[s]
		else
			println(s)
			xrange = pos:(pos+np-1)
			x0[xrange] = getfield(Pars,s)
			lb[xrange] = pars.lb[s]
			ub[xrange] = pars.ub[s]
		end
		pos += np
	end
	ftol_abs!(opt,tightness) #
	xtol_abs!(opt,tightness)
	lower_bounds!(opt,lb)
	upper_bounds!(opt,ub)
	if LBFGS==0
		min_objective!(opt,(x,g)->Criterion(x,Pars,M,vars,moms,wghts,R,lengths,TE_index,solve,true))
	else
		min_objective!(opt,(x,g)->Criterion(x,g,Pars,M,vars,moms,wghts,R,lengths,TE_index,solve,true))
	end
	maxeval!(opt,maxevals)

	return opt,x0
end

function Fit(M)
	SolveModel!(M)
	E,A,A2,XG,th,Y = MomentsBaseline(M,5,lengths,TE_index)
	# CTJF
	subplot(2,4,1)
	title("Emp. CTJF")
	plot(E[1:16],color="blue",linestyle="dashed")
	plot(E_mom[1:16],color="blue")
	plot(E[17:32],color="red",linestyle="dashed")
	plot(E_mom[17:32],color="red")
	subplot(2,4,5)
	title("Part. CTJF")
	plot(A[1:16],color="blue",linestyle="dashed")
	plot(A_mom[1:16],color="blue")
	plot(A[17:32],color="red",linestyle="dashed")
	plot(A_mom[17:32],color="red")

	# FTP
	subplot(2,4,2)
	title("Emp. FTP")
	plot(E[33:50],color="blue",linestyle="dashed")
	plot(E_mom[33:50],color="blue")
	plot(E[51:68],color="red",linestyle="dashed")
	plot(E_mom[51:68],color="red")
	subplot(2,4,6)
	title("Part. FTP")
	plot(A[33:50],color="blue",linestyle="dashed")
	plot(A_mom[33:50],color="blue")
	plot(A[51:68],color="red",linestyle="dashed")
	plot(A_mom[51:68],color="red")

	# MFIP-LR
	subplot(2,4,3)
	title("Emp. MFIP-LR")
	plot(E[69:80],color="blue",linestyle="dashed")
	plot(E_mom[69:80],color="blue")
	plot(E[81:92],color="red",linestyle="dashed")
	plot(E_mom[81:92],color="red")
	plot(E[93:104],color="green",linestyle="dashed")
	plot(E_mom[93:104],color="green")
	subplot(2,4,7)
	title("Part. MFIP-LR")
	plot(A[69:80],color="blue",linestyle="dashed")
	plot(A_mom[69:80],color="blue")
	plot(A[81:92],color="red",linestyle="dashed")
	plot(A_mom[81:92],color="red")
	plot(A[93:104],color="green",linestyle="dashed")
	plot(A_mom[93:104],color="green")

	# MFIP-RA
	subplot(2,4,4)
	title("Emp. MFIP-RA")
	plot(E[105:116],color="blue",linestyle="dashed")
	plot(E_mom[105:116],color="blue")
	plot(E[117:128],color="red",linestyle="dashed")
	plot(E_mom[117:128],color="red")
	plot(E[129:140],color="green",linestyle="dashed")
	plot(E_mom[129:140],color="green")
	subplot(2,4,8)
	title("Part. MFIP-RA")
	plot(A[105:116],color="blue",linestyle="dashed")
	plot(A_mom[105:116],color="blue")
	plot(A[117:128],color="red",linestyle="dashed")
	plot(A_mom[117:128],color="red")
	plot(A[129:140],color="green",linestyle="dashed")
	plot(A_mom[129:140],color="green")

end
function FitTE(M)
	SolveModel!(M)
	E,A,A2,XG,th,Y = MomentsBaseline(M,5,lengths,TE_index)
	# CTJF
	subplot(3,4,1)
	plot(E[17:32] .- E[1:16],color="blue",linestyle="dashed")
	plot(E_mom[17:32] .- E_mom[1:16],color="blue")
	subplot(3,4,5)
	plot(A[17:32] .- A[1:16],color="blue",linestyle="dashed")
	plot(A_mom[17:32] .- A_mom[1:16],color="blue")
	subplot(3,4,9)
	plot(Y[17:32] .- Y[1:16],color="blue",linestyle="dashed")
	plot(Y_mom[17:32] .- Y_mom[1:16],color="blue")

	# FTP
	subplot(3,4,2)
	plot(E[51:68] .- E[33:50],color="blue",linestyle="dashed")
	plot(E_mom[51:68] .- E_mom[33:50],color="blue")
	subplot(3,4,6)
	plot(A[51:68] .- A[33:50],color="blue",linestyle="dashed")
	plot(A_mom[51:68] .- A_mom[33:50],color="blue")
	subplot(3,4,10)
	plot(Y[51:68] .- Y[33:50],color="blue",linestyle="dashed")
	plot(Y_mom[51:68] .- Y_mom[33:50],color="blue")

	# MFIP-LR
	subplot(3,4,3)
	plot(E[81:92] .- E[69:80],color="blue",linestyle="dashed")
	plot(E_mom[81:92] .- E_mom[69:80],color="blue")
	plot(E[93:104] .- E[69:80],color="green",linestyle="dashed")
	plot(E_mom[93:104] .- E_mom[69:80],color="green")
	subplot(3,4,7)
	plot(A[81:92] .- A[69:80],color="blue",linestyle="dashed")
	plot(A_mom[81:92] .- A_mom[69:80],color="blue")
	plot(A[93:104] .- A[69:80],color="green",linestyle="dashed")
	plot(A_mom[93:104] .- A_mom[69:80],color="green")
	subplot(3,4,11)
	plot(Y[81:92] .- Y[69:80],color="blue",linestyle="dashed")
	plot(Y_mom[81:92] .- Y_mom[69:80],color="blue")
	plot(Y[93:104] .- Y[69:80],color="green",linestyle="dashed")
	plot(Y_mom[93:104] .- Y_mom[69:80],color="green")

	# MFIP-LR
	subplot(3,4,4)
	plot(E[117:128] .- E[105:116],color="blue",linestyle="dashed")
	plot(E_mom[117:128] .- E_mom[105:116],color="blue")
	plot(E[129:140] .- E[105:116],color="green",linestyle="dashed")
	plot(E_mom[129:140] .- E_mom[105:116],color="green")
	subplot(3,4,8)
	plot(A[117:128] .- A[105:116],color="blue",linestyle="dashed")
	plot(A_mom[117:128] .- A_mom[105:116],color="blue")
	plot(A[129:140] .- A[105:116],color="green",linestyle="dashed")
	plot(A_mom[129:140] .- A_mom[105:116],color="green")
	subplot(3,4,12)
	plot(Y[117:128] .- Y[105:116],color="blue",linestyle="dashed")
	plot(Y_mom[117:128] .- Y_mom[105:116],color="blue")
	plot(Y[129:140] .- Y[105:116],color="green",linestyle="dashed")
	plot(Y_mom[129:140] .- Y_mom[105:116],color="green")

end


function testfunc(x::Array{Float64,1},p,verbose)
	Q = (x[1]-p[1])^2 + (x[2]-p[2])^2 + (x[3]-p[3])^2
	if verbose
		println(Q)
	end
	return Q
end

function testfunc(x,g,p,verbose)
	F0 = testfunc(x,verbose)
	Δ = 0.001
	for i=1:3
		x1 = copy(x)
		x1[i] += Δ
		F1 = testfunc(x1,p,false)
		g[i] = (F1-F0)/Δ
	end
	pause(5)
	println(F0)
	return F0
end

opttest = Opt(:LD_LBFGS,3)
min_objective!(opttest,(x,g)->testfunc(x,g,[2.,3.,2.],true))
