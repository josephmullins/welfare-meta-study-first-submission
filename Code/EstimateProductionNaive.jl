# this scrip naively estimates production parameters by plugging in means for labor supply and income.
# it's a useful proof of concept to see if we get reasonable numbers

using CSV
using DataFrames

TE_moms = CSV.read("../Data/ChildTreatmentEffectsScore.csv")[1:16,:]
TE_moms[:Site] = [ones(2);2*ones(2);3*ones(7);4*ones(5)]
TE_index = convert(Array{Int64,2},TE_moms[:,[:Site,:Treatment,:AgeMin,:AgeMax]])
TE_moms0 = convert(Array{Float64,1},TE_moms.FacScore)
N_control = convert(Array{Float64,1},TE_moms.N_control)

Q_moms = CSV.read("../Data/QuarterlyMoms.csv")
E_mom = convert(Array{Float64,1},Q_moms.LFP)/100
A_mom = convert(Array{Float64,1},Q_moms.Participation)/100
A2_mom = convert(Array{Float64,1},Q_moms.Receipt)
Y_mom = convert(Array{Float64,1},Q_moms.TotInc)


τ = 0.1*ones(4,3)
XG0 = [1082.75 1351.5; 94.33 200.25; 857.67 1089.67] #<- average subsidy expense for each program
XGm = XG0/12 #<- monthly
Xcm = [57.3 71;21. 20] #<- monthly payment from CTJF and FTP members
τ[1:2,1:2] = XGm[1:2,:]./(XGm[1:2,:] .+ Xcm)
XGmom = [(XG0[:,2] .- XG0[:,1])./XG0[:,1]; XG0[3,1]/XG0[1,1]]

moms0 = [XGmom;TE_moms0]
lengths = [16,18,12,12]
year_meas = [3,4,3,3]
n_arms = [2,2,3,3]

L = zeros(4,3,18)
Y = zeros(4,3,18)
pos = 0
for i=1:4
    global pos
    ll = lengths[i]
    for s=1:n_arms[i]
        curr_slice = (pos+1):(pos+ll)
        L[i,s,1:ll] = E_mom[curr_slice]
        Y[i,s,1:ll] = Y_mom[curr_slice]
        pos += ll
    end
end

x0 = [0,-0.1,0.9,2.,0.,1.2,36.,0.2,0.4]
#x0 = [0,0.,0.7,-0.1,0.2,1.2,36.,0.4,0.4]
lb = [-5,-0.3,0.,-3.,-0.3,0.,0.,0.01,0.01]
ub = [3.,0.1,3.,4.,0.3,100.,200.,0.99,0.99]
#τ = 0.01*ones(4,3)

function Update(x)
    δI = exp.(x[1] .+ x[2]*collect(1:18))
    δθ = x[3]*ones(18)
    pc = exp.(x[4] .+ x[5]*collect(1:18))
    ϵ = x[6]
    wq = x[7]
    τ[3:4,1] .= x[8]
    τ[3,2:3] .= x[9]
    τ[4,2:3] .= x[9]
    return (δI = δI,δθ=δθ,pc=pc,ϵ=ϵ,wq=wq,τ=τ)
end

function NaiveMoments(x)
    δI = exp.(x[1] .+ x[2]*collect(1:18))
    δθ = x[3]*ones(18)
    pc = exp.(x[4] .+ x[5]*collect(1:18))
    ϵ = x[6]
    wq = x[7]
    τ[3:4,1] .= x[8]
    τ[3,2:3] .= x[9]
    τ[4,2:3] .= x[9]

    # - get treatment effects
    TE = zeros(size(TE_index)[1])
    for i=1:size(TE_index)[1]
        s = TE_index[i,1]
        arm = TE_index[i,2]
        amin = TE_index[i,3] - year_meas[s]
        amax = TE_index[i,4] - year_meas[s]
        n_age = amax+1-amin
        th0 = zeros(n_age)
        th1 = zeros(n_age)
        for a=1:n_age
            a0 = amin+a-1
            θ0 = 0.
            θ1 = 0.
            for q=1:year_meas[s]*4
                inc0 = Y[s,1,q]
                h0 = 30*L[s,1,q]
                inc1 = Y[s,arm+1,q]
                h1 = 30*L[s,arm+1,q]
                Age_Year = a0+floor(Int64,q/4)
                    pc0 = (1-τ[s,1])*pc[Age_Year]
                    pc1 = (1-τ[s,1+arm])*pc[Age_Year]
                if Age_Year<18
                    θ0 = δI[Age_Year]*log(inc0 + wq*(112-h0)) - 1/(1-ϵ)* δI[Age_Year]*(L[s,1,q]*log(112-30 + 30*pc0^(1-ϵ)) + (1-L[s,1,q])*log(112)) + δθ[Age_Year]*θ0
                    θ1 = δI[Age_Year]*log(inc1 + wq*(112-h1)) - 1/(1-ϵ)* δI[Age_Year]*(L[s,arm+1,q]*log(112-30 + 30*pc1^(1-ϵ)) + (1-L[s,arm+1,q])*log(112)) + δθ[Age_Year]*θ1
                end
            end
            th0[a] = θ0
            th1[a] = θ1
        end
        TE[i] = mean(th1 .- th0)
    end
    XG = zeros(3,2)
    for s=1:3
        for arm=1:2
            Xc = zeros(12)
            for q=1:12
                inc1 = Y[s,arm,q]
                h1 = 30*L[s,arm,q]
                pc1 = (1-τ[s,arm])*mean(pc)
                Xc[q] = 30*pc1^(1-ϵ)/(112 -30 + 30*pc1^(1-ϵ))*inc1*L[s,arm,q]
            end
            XG[s,arm] = τ[s,arm]/(1-τ[s,arm])*mean(Xc)
        end
    end
    XG_moms = [(XG[i,2]-XG[i,1])/XG[i,1] for i=1:3]
    XG_mom2 = XG[3,1]/XG[1,1]
    return [XG_moms;XG_mom2;TE]
end
wghts = ones(20)
wghts[1:4] .= 0.
wghts = [100*ones(4);N_control]

function Criterion(x)
    moms = NaiveMoments(x)
    Q = sum(wghts .* (moms .- moms0).^2)
    println(Q)
    return Q
end

function Criterion(x,g)
    Δ = 1e-7
    F0 = Criterion(x)
    for i=1:length(x)
        x1 = copy(x)
        x1[i] += Δ
        F1 = Criterion(x1)
        g[i] = (F1-F0)/Δ
    end
    println(F0)
    return F0
end

Criterion(x0)
opt = Opt(:LD_LBFGS,9)
lower_bounds!(opt,lb)
upper_bounds!(opt,ub)
min_objective!(opt,(x,g)->Criterion(x,g))

res = optimize(opt,x0)

pars = Update(res[2])

break
pc = 10.
pc0 = (1-τ[1,1])*pc
pc1 = (1-τ[1,2])*pc
ϵ = 10.
inc0 = Y[1,1,1]
inc1 = Y[1,2,1]
h0 = L[1,1,1]*30
h1 = L[1,2,1]*30
th0 = log(inc0 + wq*(112-h0)) - 1/(1-ϵ)*log(112-h0 + h0*pc0^(1-ϵ))
th1 = log(inc1 + wq*(112-h1)) - 1/(1-ϵ)*log(112-h1 + h1*pc1^(1-ϵ))

function TestOut(E,wq,ϵ,pc,h1,h0)
    th0 = log(h0*E + wq*(112-h0)) - (1/(1-ϵ))*log(112-h0 + h0*pc^(1-ϵ))
    th1 = log(h1*E + wq*(112-h1))- (1/(1-ϵ))*log(112-h1 + h1*pc^(1-ϵ))
    return th1-th0
end
