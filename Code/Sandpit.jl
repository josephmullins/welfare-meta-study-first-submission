include("BaselineModel.jl")
using DelimitedFiles


earnings = reshape(readdlm("earnings"),3,72)
budget = reshape(readdlm("budget"),3,3,3,17,17*4,2,2,2)
budget_ageout = reshape(readdlm("budget_ageout"),3,72,2,2)
#findmax(budget1[:,3,:,:,:,:,:,:])

TimeLimit_Ind=[false true true; false true true; false false false]

TimeLimits=[0 7 7; 0 8 8; 0 0 0]

τ=ones(3,3).*0.5
