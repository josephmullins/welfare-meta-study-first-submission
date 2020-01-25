using CSV
using DataFrames
using NLopt
using Distributions

D = CSV.read("../Data/ChildTreatmentEffects.csv")

measures = [:Achievement,:AchieveBelowAverage,:BPI,:PB,:Repeat]
J = length(measures)
V = zeros(J,J)

for i=1:J,j=1:J
    Ik = ismissing.(D[measures[i]]) .| ismissing.(D[measures[j]])
    V[i,j] = cov(D[.~Ik,measures[i]],D[.~Ik,measures[j]])
end

function DistanceFunc(x,V,J)
    v = x[1]
    λ = [1;x[2:J]]
    σ = x[J+1:end]
    Q = 0
    for i=1:J
        Q += (λ[i]^2*v + σ[i] - V[i,i])^2
        for j=1:i-1
            Q += (λ[i]*λ[j]*v - V[i,j])^2
        end
    end
    return Q
end

function DistanceFunc(x,g,V,J)
    f0 = DistanceFunc(x,V,J)
    Δ = 1e-7
    for i=1:length(x)
        x1 = copy(x)
        x1[i] += Δ
        f1 = DistanceFunc(x1,V,J)
        g[i] = (f1-f0)/Δ
    end
    return f0
end

x0 = ones(2*J)
DistanceFunc(x0,V,J)

#opt = Opt(:LN_NELDERMEAD,2*J)
opt = Opt(:LD_LBFGS,2*J)
lower_bounds!(opt,[0;-100*ones(J-1);zeros(J)])
upper_bounds!(opt,Inf*ones(2*J))
min_objective!(opt,(x,g)->DistanceFunc(x,g,V,J))

res = optimize(opt,x0)
res2 = optimize(opt,1.5*ones(2*J))
res3 = optimize(opt,2*ones(2*J))
# can confirm starting points wind up in same spot

λ = [1;res[2][2:J]]
σ = res[2][J+1:end]
wght = (λ.^2)./σ
D[:FacScore] = 0.

# don't want factor score!!!
for i=1:size(D)[1]
    denom = 0
    for j=1:J
        if ~ismissing(D[i,measures[j]])
            D[i,:FacScore] += D[i,measures[j]]/λ[j]
            denom += 1
        end
    end
    D[i,:FacScore] /= denom
end

CSV.write("../Data/ChildTreatmentEffectsScore.csv",D)
