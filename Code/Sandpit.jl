include("BaselineModel.jl")
using DelimitedFiles


earnings = reshape(readdlm("earnings"),3,72)
budget = reshape(readdlm("budget"),3,3,3,17,17*4,2,2,2)
budget_ageout = reshape(readdlm("budget_ageout"),3,72,2,2)
#findmax(budget1[:,3,:,:,:,:,:,:])

TimeLimit_Ind=[false true true; false true true; false false false]

TimeLimits=[0 7 7; 0 8 8; 0 0 0]

Work_Reqs_Ind=[false true true; false true true; false false true];

τ=ones(3,3).*0.5

@time Mod1=initialize_model()

Mod1.αc
UpdateSpecificParams!(Mod1; αc=0.3)
Mod1.αc
UpdateSpecificParams!(Mod1)
Mod1.αc



@time Mod2=initialize_model()



Mod2.αWR




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

@time S1_t=Simulate(Mod2,10000,30,1,3)
@time S2_t=Simulate(Mod2,10000,30,2,3)
@time S3_t=Simulate(Mod2,10000,30,3,3)
@time S3_t2=Simulate(Mod2,10000,30,3,2)

budget1[3,3,1,3,25,2,2,2]
Earnings[3,25]

mean(S1_t[2])
mean(S2_t[2])
mean(S3_t[2])
mean(S3_t2[2])

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
