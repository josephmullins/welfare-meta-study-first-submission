# this assumes we have run the script "SetupBaseline.jl"

measures = [:Achievement,:AchieveBelowAverage,:Repeat]
D = CSV.read("../Data/ChildTreatmentEffects.csv")

# load and read in moments
# how do we want to do this
# ISSUE: moments are not balanced across treatment arms, so this is a problem
TEs = []
arms = []
a0s = []
a1s = []
wghts = []

moms_collect = []
for s in site_list
    m = D[D.Site.==String(s),:]
    TE = D[:,measures]
    wght = ones(size(TE))
    Idrop = convert(Array{Bool,2},ismissing.(TE))
    TE = convert(Array{Float64,2},coalesce.(TE,0.))
    wght[Idrop] .= 0.
    arm = convert(Array{Int64,1},D.Treatment) .+ 1
    #append!(arms,[arm])
    #append!(wghts,[wght])
    #append!(TEs,[TE])
    append!(moms_collect,[(TE=TE,wght=wght,arm=arm)])
end
moments = (;zip(site_list,moms_collect)...)
