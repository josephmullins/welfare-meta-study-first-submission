# this assumes we have run the script "SetupBaseline.jl"
include("ProductionEstimation.jl")
measures = [:Achievement,:AchieveBelowAverage,:Repeat]
D = CSV.read("../Data/ChildTreatmentEffects.csv")

# create a moments object that stores all the info we want
moms_collect = []
for s in site_list
    m = D[D.Site.==String(s),:]
    TE_ = m[:,measures]
    wght = ones(size(TE_))
    Idrop = [ismissing.(TE_[v]) for v in measures] #convert(Array{Bool,2},ismissing.(TE))
    TE = zeros(size(TE_)) #convert(Array{Float64,2},coalesce.(TE,0.))
    for i=1:length(measures)
        TE[:,i] = convert(Array{Float64,1},coalesce.(TE_[measures[i]],0.))
        wght[Idrop[i],i] .= 0.
    end
    arm = convert(Array{Int64,1},m.Treatment) .+ 1
    a0 = convert(Array{Int64,1},m.AgeMin)
    a1 = convert(Array{Int64,1},m.AgeMax)
    append!(moms_collect,[(TE=TE,wght=wght,arm=arm,a0=a0,a1=a1)])
end
moments = (;zip(site_list,moms_collect)...)

# Initialize the parameters
pars_prod = ProdPars(length(measures))

CP = GetChoiceProbsAll(pars2,site_list,budget,site_features);

ProductionCriterion(pars_prod,pars2,CP,site_list[1:5],budget,moments,site_features)

opt,x0 = GetOptimization(pars_prod,pars2,CP,site_list[1:5],budget,moments,site_features)

res = optimize(opt,x0)
vlist = [:gN,:gF,:δI,:δθ,:λ]
pars_prod = UpdatePars(res[2],pars_prod,vlist)
