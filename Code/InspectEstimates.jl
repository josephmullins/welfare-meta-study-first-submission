include("BaselineModel.jl")
include("MCMCRoutines.jl")
include("CounterfactualRoutines.jl")
include("SetupBaseline.jl")
using Printf
using Random
Γ = readdlm("FirstStageSIPP/Gamma_est")[:]
gF = readdlm("FirstStageSIPP/gF")[:] #<- we need to figure this out and decide what we're doing!!!

xm0,mpars = ModPars(Γ,gF)
xh0,hpars = HyperPars()

Pm = readdlm("BaselineChain")
Ph = readdlm("BaselineChainHyper")

post_est = [mode(Pm[i,20000:end]) for i=1:size(Pm)[1]]
writedlm("posterior_ests",post_est)

# test
nsite = length(site_features.site_list)
nb = 100
EffectA = zeros(nb,nsite,2,5)
EffectH = zeros(nb,nsite,2,5)
ElastA = zeros(nb,nsite,2,5)
ElastH = zeros(nb,nsite,2,5)
ElastF = zeros(nb,nsite,2,5)


Random.seed!(202020)
#ii = rand(10000:50000,nb)
for j=1:nb
    i = rand(10000:50000)
    EffectA[j,:,:,:],EffectH[j,:,:,:],ElastA[j,:,:,:],ElastH[j,:,:,:],ElastF[j,:,:,:] = GetEffects(Pm[:,i],mpars,budget,site_features)
end

# save these results to file
D = DataFrame()
E = DataFrame()
for i=1:8
    T = site_features.T[i]

    d1=DataFrame(year=1:T,var="Participation",Site=site_str[i],case="Work Requirement",Effect=mean(EffectA[:,i,1,1:T],dims=1)[:],lb=[quantile(EffectA[:,i,1,t],0.025) for t=1:T],ub=[quantile(EffectA[:,i,1,t],0.975) for t=1:T])
    d2=DataFrame(year=1:T,var="Participation",Site=site_str[i],case="Time Limit",Effect=mean(EffectA[:,i,2,1:T],dims=1)[:],lb=[quantile(EffectA[:,i,2,t],0.025) for t=1:T],ub=[quantile(EffectA[:,i,2,t],0.975) for t=1:T])
    append!(D,[d1;d2])
    d1=DataFrame(year=1:T,var="LFP",Site=site_str[i],case="Work Requirement",Effect=mean(EffectH[:,i,1,1:T],dims=1)[:],lb=[quantile(EffectH[:,i,1,t],0.025) for t=1:T],ub=[quantile(EffectH[:,i,1,t],0.975) for t=1:T])
    d2=DataFrame(year=1:T,var="LFP",Site=site_str[i],case="Time Limit",Effect=mean(EffectH[:,i,2,1:T],dims=1)[:],lb=[quantile(EffectH[:,i,2,t],0.025) for t=1:T],ub=[quantile(EffectH[:,i,2,t],0.975) for t=1:T])
    append!(D,[d1;d2])

    e1=DataFrame(year=1:T,var="Formal Care",Site=site_str[i],case="Wage Change",Elasticity=mean(ElastF[:,i,1,1:T],dims=1)[:]/0.1,lb=[quantile(ElastF[:,i,1,t],0.025)/0.1 for t=1:T],ub=[quantile(ElastF[:,i,1,t],0.975)/0.1 for t=1:T])
    e2=DataFrame(year=1:T,var="Formal Care",Site=site_str[i],case="Price Change",Elasticity=mean(ElastF[:,i,2,1:T],dims=1)[:]/0.1,lb=[quantile(ElastF[:,i,2,t],0.025)/0.1 for t=1:T],ub=[quantile(ElastF[:,i,2,t],0.975)/0.1 for t=1:T])
    append!(E,[e1;e2])
    e1=DataFrame(year=1:T,var="LFP",Site=site_str[i],case="Wage Change",Elasticity=mean(ElastH[:,i,1,1:T],dims=1)[:]/0.1,lb=[quantile(ElastH[:,i,1,t],0.025)/0.1 for t=1:T],ub=[quantile(ElastH[:,i,1,t],0.975)/0.1 for t=1:T])
    e2=DataFrame(year=1:T,var="LFP",Site=site_str[i],case="Price Change",Elasticity=mean(ElastH[:,i,2,1:T],dims=1)[:]/0.1,lb=[quantile(ElastH[:,i,2,t],0.025)/0.1 for t=1:T],ub=[quantile(ElastH[:,i,2,t],0.975)/0.1 for t=1:T])
    append!(E,[e1;e2])


end
CSV.write("/Users/joseph/Dropbox/Research Projects/WelfareMetaAnalysis/Figures/Effects.csv",D)
CSV.write("/Users/joseph/Dropbox/Research Projects/WelfareMetaAnalysis/Figures/Elasticities.csv",E)

break
# Write Gobal Parameters to file
# αC, σH,σC,σF, β,w_q
vlist1 = [:αc,:σH,:σC,:wq,:β]
vstring1 = ["\$\\alpha_C\$","\$\\sigma_H\$","\$\\sigma_F\$","\$w_q\$","\$\\beta\$"]
vlist2 = [:αH, :σαH, :αA, :σαA, :αF, :σαF, :αWR, :σαWR, :αWR2, :σαWR2]
vstring2 = ["\$\\alpha_H\$","\$\\sigma_H\$","\$\\alpha_A\$","\$\\sigma_A\$","\$\\alpha_F\$","\$\\sigma_F\$","\$\\alpha_{WR}\$","\$\\sigma_{WR}\$","\$\\alpha_{WR2}\$","\$\\sigma_{WR2}\$"]

io = open("/home/joseph/Dropbox/Research Projects/WelfareMetaAnalysis/Tables/EstsTable.tex", "w");

for i=1:5
    s = vstring1[i]
    ii = mpars.pos[vlist1[i]]
    ests = [quantile(Pm[ii,20000:end],0.025),quantile(Pm[ii,20000:end],0.25),mode(Pm[ii,20000:end]),quantile(Pm[ii,20000:end],0.75),quantile(Pm[ii,20000:end],0.975)]
    #println(ests)
    for j=1:5
        s = string(s," & ",@sprintf("%0.2f",ests[j]))
    end
    s = string(s,"\\\\")
    write(io,s)
end

for i=1:10
    s = vstring2[i]
    ii = hpars.pos[vlist2[i]]
    ests = [quantile(Ph[ii,20000:end],0.025),quantile(Ph[ii,20000:end],0.25),mode(Ph[ii,20000:end]),quantile(Ph[ii,20000:end],0.75),quantile(Ph[ii,20000:end],0.975)]
    #println(ests)
    for j=1:5
        s = string(s," & ",@sprintf("%0.2f",ests[j]))
    end
    s = string(s,"\\\\")
    write(io,s)
end
close(io)


# αA,σA,αH,σH,α_F,σ_F,αR,σR,αR2,σR2
# write to table site-specific parameters
vlist = [:αH,:αA,:αF,:αWR,:αWR2]
vstring2 = ["\$\\alpha_H\$","\$\\alpha_A\$","\$\\alpha_F\$","\$\\alpha_{WR}\$","\$\\alpha_{WR2}\$"]
io = open("/home/joseph/Dropbox/Research Projects/WelfareMetaAnalysis/Tables/SiteEstsTable.tex", "w");
for i=1:8
    s = site_str[i]
    for j=1:5
        ii = mpars.pos[vlist[j]][i]
        est = mode(Pm[ii,20000:end])
        s = string(s," & ",@sprintf("%0.2f",est))
    end
    s = string(s,"\\\\")
    write(io,s)
    s = ""
    for j=1:5
        ii = mpars.pos[vlist[j]][i]
        sdev = std(Pm[ii,20000:end])
        s = string(s," & ","(",@sprintf("%0.2f",sdev),")")
    end
    s = string(s,"\\\\")
    write(io,s)
end
close(io)
