include("BaselineModel.jl")
include("MCMCRoutines.jl")
include("CounterfactualRoutines.jl")
include("SetupBaseline.jl")
using Printf
Γ = readdlm("FirstStageSIPP/Gamma_est")[:]
gF = readdlm("FirstStageSIPP/gF")[:] #<- we need to figure this out and decide what we're doing!!!

xm0,mpars = ModPars(Γ,gF)
xh0,hpars = HyperPars()

Pm = readdlm("BaselineChain")
Ph = readdlm("BaselineChainHyper")

# test
nsite = length(site_features.site_list)
nb = 10
EffectA = zeros(nb,nsite,2,5)
EffectH = zeros(nb,nsite,2,5)
ElastA = zeros(nb,nsite,2,5)
ElastH = zeros(nb,nsite,2,5)
ElastF = zeros(nb,nsite,2,5)

for i=1:nb
    EffectA[i,:,:,:],EffectH[i,:,:,:],ElastA[i,:,:,:],ElastH[i,:,:,:],ElastF[i,:,:,:] = GetEffects(Pm[:,i],mpars,budget,site_features)
end


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
    ests = [quantile(Pm[ii,5000:end],0.025),quantile(Pm[ii,5000:end],0.25),mode(Pm[ii,5000:end]),quantile(Pm[ii,5000:end],0.75),quantile(Pm[ii,5000:end],0.975)]
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
    ests = [quantile(Ph[ii,5000:end],0.025),quantile(Ph[ii,5000:end],0.25),mode(Ph[ii,5000:end]),quantile(Ph[ii,5000:end],0.75),quantile(Ph[ii,5000:end],0.975)]
    #println(ests)
    for j=1:5
        s = string(s," & ",@sprintf("%0.2f",ests[j]))
    end
    s = string(s,"\\\\")
    write(io,s)
end


close(io)


# αA,σA,αH,σH,α_F,σ_F,αR,σR,αR2,σR2
