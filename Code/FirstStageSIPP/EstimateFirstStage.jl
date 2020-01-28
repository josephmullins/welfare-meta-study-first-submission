include("MLEmoduleSpec1.jl")

D = CSV.read("SIPPdata.csv")
NS = 48

#age::Int64,wage,price,work,pay,state
age = convert(Array{Int64,1},D.age_youngest)
wage = log.(convert(Array{Float64,1},D.Wage))
price = convert(Array{Float64,1},D.ChCare)
work = convert(Array{Int64,1},(D.Work.=="TRUE") .& (D.Wage.>0))
pay = convert(Array{Int64,1},D.Pay.=="TRUE")
state = convert(Array{Int64,1},D.state)

# αc = x[1]
# μF = x[2]
# wq = x[3]
# σw = x[4]
# σP = x[5]
# σC = x[6]
# Γ = x[7:21] #<- 15 parameters (for now, fix later)
# gN = x[22:23]
# gF = x[24:25]
# μH = x[26:(25+NS)]
# mwage = x[(25+NS+1):(25+2*NS)]
# pF = x[(25+2*NS+1):(25+3*NS)]
dw = sort(by(D,:state,x->mean(log.(x.Wage[x.Wage.>0]))))
w0 = exp.(convert(Array{Float64,1},dw.x1))

x0 = [ones(NS);zeros(NS);zeros(NS);40*ones(NS);w0;
    [3.,1.,10.,1.,-1.,0.,0.,0.,0.,0.,0.]]

#x0 = [[1.,0.,3.,1.,10.,1.];
 #[-1.,0.,0.];zeros(2);zeros(2);zeros(NS);w0;40*ones(NS)
 #]

lb = [[-Inf*ones(NS);-Inf*ones(NS);-Inf*ones(NS);zeros(NS);4*30*ones(NS)];
    [0.,0.,0.,0.,-10.,-.1,0.]; -Inf*ones(4)]

ub = [Inf*ones(NS*5);[5,Inf,Inf,Inf,1.,1.,0.]; Inf*ones(4)]


#ub = [[Inf,Inf,5.,Inf,Inf,Inf];
        #[1.,1.,0.];Inf*ones(4);Inf*ones(NS);20*30*ones(NS);90*ones(NS)]

#x0 = GetParVec(mod)

g = zeros(length(x0))

#LL = LogLike(x0,g,mod,age,wage,price,work,pay,state)


LL = LogLike(x0,g,age,wage,price,work,pay,state)

#g = ForwardDiff.gradient(x->LogLike(x,age,wage,price,work,pay,state),x0)

opt = Opt(:LD_LBFGS,length(x0))
max_objective!(opt,(x,g)->LogLike(x,g,age,wage,price,work,pay,state))
lower_bounds!(opt,lb)
upper_bounds!(opt,ub)
res = optimize(opt,x0)

pars = GetPars(res[2]);
D[:pPay] = 0.
for i=1:length(state)
    G = exp(pars.Γ[1] + pars.Γ[2]*age[i] + pars.Γ[3]*age[i]^2)
    st = state[i]
    mwagei = pars.mwage[st]
    if age[i]<=5
        gNi = pars.gN[1]
        gFi = pars.gF[1]
    else
        gNi = pars.gN[2]
        gFi = pars.gF[2]
    end
    v0 = (pars.αc + G)*log(pars.wq*112)
    v1N = (pars.αc + G)*log(mwagei + (112-30)*pars.wq) - pars.μH[st] - G*gNi
    v1F = (pars.αc + G)*log(mwagei + (112-30)*pars.wq - pars.pF[st]) + pars.μF - pars.μH[st] - G*gFi
    #v1 = σC*log(exp(v1N/σC)+exp(v1F/σC))
    D[i,:pPay] = 1/(1+exp((v1N-v1F)/pars.σC))
end
