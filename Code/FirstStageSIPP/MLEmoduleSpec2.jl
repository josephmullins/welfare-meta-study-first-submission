using ForwardDiff
using NLopt
using CSV
using DataFrames
using LinearAlgebra

# spec 2 is same a spec 1, but ignores work decision so fewer parameters to estimate


function LogLike(x,age,wage,price,work,pay,state)
    NS = 48
    P = reshape(x[1:4*NS],NS,4)
    μH = P[:,1]
    μF = P[:,2]
    pF = P[:,3]
    mwage = P[:,4]
    wq = x[4*NS+1]
    σw = x[4*NS+2]
    σP = x[4*NS+3]
    Γ = x[4*NS .+ (4:6)] #<- 3 parameters, polynomial in age
    #gN = x[4*NS .+ (7:8)]
    gF = x[4*NS .+ (7:8)]
    ll = 0
    for i=1:length(state)
        G = exp(Γ[1] + Γ[2]*age[i] + Γ[3]*age[i]^2)
        st = state[i]
        mwagei = mwage[st]
        if age[i]<=5
            #gNi = gN[1]
            gFi = gF[1]
        else
            #gNi = gN[2]
            gFi = gF[2]
        end
        v1N = G*log(mwagei + (112-30)*wq) + μH[st]
        v1F = G*log(mwagei + (112-30)*wq - pF[st]) + μF[st] + G*gFi
        v1 = log(exp(v1N)+exp(v1F))
        # alternative: use actual income here
        if work[i]==1
            # wage contribution
            ll += -0.5*((wage[i]-log(mwagei))/σw)^2 - log(σw)
            if pay[i]==0
                ll += (v1N - v1)
            else
                ll += (v1F - v1)
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
    P = reshape(x[1:4*NS],NS,4)
    res = (μH = P[:,1],
    μF = P[:,2],
    pF = P[:,3],
    mwage = P[:,4],
    wq = x[4*NS+1],
    σw = x[4*NS+2],
    σP = x[4*NS+3],
    Γ = x[4*NS .+ (4:6)], #<- 3 parameters, polynomial in age
    #gN = x[4*NS .+ (7:8)],
    gF = x[4*NS .+ (7:8)])
    return res
end
