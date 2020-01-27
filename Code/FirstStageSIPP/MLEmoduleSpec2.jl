using ForwardDiff
using NLopt
using CSV
using DataFrames
using LinearAlgebra

function LogLike(x,age,wage,price,work,pay,state)
    NS = 48
    αc = x[1]
    μF = x[2]
    wq = x[3]
    σw = x[4]
    σP = x[5]
    σC = x[6]
    Γ = x[7:9] #<- 3 parameters, polynomial in age
    gN = x[10:11]
    gF = x[12:13]
    μH = x[14:(13+NS)]
    mwage = x[(13+NS+1):(13+2*NS)]
    pF = x[(13+2*NS+1):(13+3*NS)]
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
        v0 = (αc + G)*log(wq*112)
        v1N = (αc + G)*log(mwagei + (112-30)*wq) - μH[st] - G*gNi
        v1F = (αc + G)*log(mwagei + (112-30)*wq - pF[st]) + μF - μH[st] - G*gFi
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
    res = (αc = x[1],
    μF = x[2],
    wq = x[3],
    σw = x[4],
    σP = x[5],
    σC = x[6],
    Γ = x[7:9], #<- 3 parameters, polynomial in age
    gN = x[10:11],
    gF = x[12:13],
    μH = x[14:(13+NS)],
    mwage = x[(13+NS+1):(13+2*NS)],
    pF = x[(13+2*NS+1):(13+3*NS)])
    return res
end
