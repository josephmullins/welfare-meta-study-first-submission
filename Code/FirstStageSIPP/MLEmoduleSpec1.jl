using ForwardDiff
using NLopt
using CSV
using DataFrames
using LinearAlgebra

# spec 1 has state-specific αc,μF,μH,wage,price of care


function LogLike(x,age,wage,price,work,pay,state)
    NS = 48
    P = reshape(x[1:5*NS],NS,5)
    αc = P[:,1]
    μH = P[:,2]
    μF = P[:,3]
    pF = P[:,4]
    mwage = P[:,5]
    wq = x[5*NS+1]
    σw = x[5*NS+2]
    σP = x[5*NS+3]
    σC = x[5*NS+4]
    Γ = x[5*NS .+ (5:7)] #<- 3 parameters, polynomial in age
    gN = x[5*NS .+ (8:9)]
    gF = x[5*NS .+ (10:11)]
    ll = 0
    for i=1:length(state)
        G = exp(Γ[1] + Γ[2]*age[i] + Γ[3]*age[i]^2)
        st = state[i]
        mwagei = mwage[st]
        if age[i]<=5
            gNi = gN[1]
            gFi = gF[1]
        else
            gNi = gN[2]
            gFi = gF[2]
        end
        v0 = (αc[st] + G)*log(wq*112)
        v1N = (αc[st] + G)*log(mwagei + (112-30)*wq) - μH[st] - G*gNi
        v1F = (αc[st] + G)*log(mwagei + (112-30)*wq - pF[st]) + μF[st] - μH[st] - G*gFi
        v1 = σC*log(exp(v1N/σC)+exp(v1F/σC))
        if work[i]==0
            ll += v0 - log(exp(v0)+exp(v1))
        else
            ll += v1 - log(exp(v0)+exp(v1))
            # wage contribution
            ll += -0.5*((wage[i]-log(mwagei))/σw)^2 - log(σw)
            if pay[i]==0
                ll += (v1N - v1)/σC
            else
                ll += (v1F - v1)/σC
                ll += -0.5*((price[i]-pF[st])/σP)^2 - log(σP)
            end
        end
    end
    return ll
end


function LogLike(x,g,age,wage,price,work,pay,state)
    LL = LogLike(x,age,wage,price,work,pay,state)
    g[:] = ForwardDiff.gradient(x->LogLike(x,age,wage,price,work,pay,state),x)
    return LL
end

function GetPars(x)
    NS = 48
    P = reshape(x[1:5*NS],NS,5)
    res = (αc = P[:,1],
    μH = P[:,2],
    μF = P[:,3],
    pF = P[:,4],
    mwage = P[:,5],
    wq = x[5*NS+1],
    σw = x[5*NS+2],
    σP = x[5*NS+3],
    σC = x[5*NS+4],
    Γ = x[5*NS .+ (5:7)],
    gN = x[5*NS .+ (8:9)],
    gF = x[5*NS .+ (10:11)])
    return res
end
