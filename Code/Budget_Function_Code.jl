# budget matrix
cd("/Users/FilipB/github/welfare-meta-study/Code")

using CSV
using GLM
using DataFrames
using DelimitedFiles
budget1=Array{Float64,8}# budget(site,arm,NumChild,age0,quarter,program,work, eligible) (NS x NT x 3 x 17 x Q x 2 x 2 x 2)

Quarterly_Data=CSV.read("../Data/QuarterlyData.csv")

BenStd=CSV.read("../Data/WelfareRules/BenStd.csv")

PovGuideline=CSV.read("../Data/WelfareRules/PovGuideline.csv")

SNAPRules=CSV.read("../Data/WelfareRules/SNAPRules.csv")

#=
I observe the following rules:

    - 1 indexes the control, 2 the treatment
    - 1 indexes Connecticut, 2 Florida, 3 Minneapolis
    - 1 indexes ineligible, 2 eligible

    There are 17*4=68 quarters

=#
println("Checkpoint 1")

Dev_Years=18

budget1=zeros(4,3,3,Dev_Years,Dev_Years*4,2,2,2)

site=["CTJF" "FTP" "MFIP-LR-F" "MFIP-R-F"]


# I store the first year I observed a program. This will come in handy later.
First_Year=zeros(Int64,4)
for i in 1:4
    Control = Quarterly_Data[Quarterly_Data[:, :Site] .== site[i], :]
    First_Year[i]=(minimum(Control.Year))
end

Earnings=zeros(4,Dev_Years*4+1) # I add some extra years for the parents whose kids will age out

for i in 1:4
    Control = Quarterly_Data[Quarterly_Data[:, :Site] .== site[i], :]
    Control = Control[Control[:, :Treatment] .== "Control", :]
    Control.Monthly_Earnings = Control.Earnings./(3*Control.LFP/100)
    Control.Period = 1:size(Control)[1]
    lm2 = lm(@formula(Monthly_Earnings ~ Period), Control)
    observed=length(Control.Monthly_Earnings)
    for j in 1:observed
        Earnings[i,j]=Control[j,:Monthly_Earnings]
    end
    for j in (observed+1):length(Earnings[1,:])
        Earnings[i,j]=coef(lm2)[1]+j*coef(lm2)[2]
    end
end

#=

SNAP payouts next

Some of this will be a little ugly:
CTJF started 2 years after the others. Hence, I need to shift some years by 2.

=#


SNAP=zeros(3,4,Dev_Years*4) # first index is the number of kids, second is site, third is quarter
for i in 1:3
    SNAPRulesC = SNAPRules[SNAPRules[:, :NumChild] .== i, :]
    for j in 1:4
        SNAPRulesC2 = SNAPRulesC[SNAPRulesC[:, :year] .>= First_Year[j], :]
            for k in 1:Dev_Years
                for z in 1:4
                SNAP[i,j,(4*k+z-4)]=SNAPRulesC2.MA[k]
                end
            end
    end
end


# need to consider parents whose kids aged out
SnapRules_0=SNAPRules[SNAPRules[:, :NumChild] .== 0, :]
SNAP0=zeros(4,Dev_Years*4+1)
for i in 1:4

    SnapRules_02 = SnapRules_0[SnapRules_0[:, :year] .>= First_Year[i], :]
            for k in 1:18
                for z in 1:4
                    SNAP0[i,(4*k+z-4)]=SnapRules_02.MA[k]
                end
            end
            SNAP0[i,73]=SnapRules_02.MA[19] # one more year for parents with aged-out kids
end


#=

Poverty Cutoff

=#

Poverty=zeros(3,4,Dev_Years*4+1)# first index is the number of kids, second is site, third is quarter for compatibility
PovGuideline
for j in 1:4
    #PovGuideline2=PovGuideline[:,(i+1)]
    PovGuideline2=PovGuideline[PovGuideline[:, :year] .>= First_Year[j], :] # but we do site first
    for i in 1:3
        for k in 1:Dev_Years
            for z in 1:4
            Poverty[i,j,(4*k+z-4)]=PovGuideline2[k,(i+2)]./12 # the index of this spreadsheet is family size, so 2=1+1kid
            end
        end
    end
end


#=

Benefits Guideline

=#
println("Checkpoint 2")
State=["Connecticut" "Florida" "Minnesota" "Minnesota"]
YearVec=[1996 1994 1994 1994]

BenStd
BS2=BenStd[BenStd[:, :NumChild] .== 1, :]
BS2=BS2[BS2[:, :state] .== State[1], :]
deletecols!(BS2, :NumChild)
BS3=melt(BS2,:state)
sort!(BS3, :variable)
BS3[!,:YearString]=string.(BS3.variable)
BS3[!,:Year]=parse.(Int64, BS3.YearString)
BS3=BS3[BS3[:, :Year] .>= 1996, :]
BS3.value[1]

Benefit=zeros(3,4,Dev_Years*4)
for i in 1:3
    BS2=BenStd[BenStd[:, :NumChild] .== i, :]
    for j in 1:4
        BS3=BS2[BS2[:, :state] .== State[j], :]
        deletecols!(BS3,:NumChild)
        BS3=melt(BS3,:state)
        sort!(BS3, :variable)
        BS3.YearString=string.(BS3.variable)
        BS3.Year=parse.(Int64, BS3.YearString)
        BS3=BS3[BS3[:, :Year] .>= 1994, :] # fix this for CT later on!!!
        if i!=1
        for k in 1:17
            for z in 1:4
                Benefit[i,j,(4*k+z-4)]=BS3.value[k] # I assume this was quarterly?
            end
        end
        else # this is ugly but I need to feed in something else for CT
            BS3=BS3[BS3[:, :Year] .>= 1996, :]
                for k in 1:(Dev_Years-2) # snip off the last 2 years
                    for z in 1:4
                        Benefit[i,j,(4*k+z-4)]=BS3.value[k] # I assume this was quarterly?
                    end
                end
                for k in (Dev_Years-2):Dev_Years
                    for z in 1:4
                        Benefit[i,j,(4*k+z-4)]=BS3.value[Dev_Years-2] # I assume this was quarterly?
                    end
                end
        end
    end
end
Benefit
typeof(Benefit)

# index 1 indicates false, 2 true
Eligible=[0, 1]
Work=[0,1]
Program=[0,1]




# budget(site,arm,NumChild,age0,quarter,eligible,program,work) (NS x NT x 3 x 17 x Q x 2 x 2 x 2)
budget1=zeros(4,3,3,Dev_Years,Dev_Years*4,2,2,2)
Foodstamps_receipt=zeros(4,3,3,Dev_Years,Dev_Years*4,2,2,2)
3*3*3*Dev_Years*Dev_Years*4*2*2*2

println("Checkpoint 3")
@time @inbounds @simd for a0 in 1:Dev_Years # outermost loop is for child initial age, 0-17 (shifted by 1-indexing)
    @inbounds @simd    for q in 1:(Dev_Years*4) # next loop is for quarter
    @inbounds @simd         for nk in 1:3 # num kids
    @inbounds @simd             for e in 1:2 # eligible or not
    @inbounds @simd                 for w in 1:2 # working or not
    @inbounds @simd                     for p in 1:2 # program, aka participating or not
    @inbounds @simd                         for site in 1:4 # loop over controls
                                                AFDC=Eligible[e]*max(Benefit[nk,site,q]-(1-0.33)*max(Earnings[site,q]*Work[w]-120,0),0)
                                                Foodstamps=max(SNAP[nk,site, q]-0.3*max(0.8*Earnings[site,q]*Work[w]+AFDC-134,0),0)
                                                budget1[site,1,nk,a0,q,e,p,w]=Program[p]*(AFDC+Foodstamps)+
                                                                            Earnings[site,q]*Work[w]
                                                Foodstamps_receipt[site,1,nk,a0,q,e,p,w]=Foodstamps

                                # CTJF treatment
                                if site==1 # next I fill the treatment arms in one by one
                                    Too_rich=0
                                    if Earnings[site,q]>Poverty[nk, site, q] # reminder that first index
                                        Too_rich=Work[w] # priced out of benefits AND food stamps if you pass a threshold
                                    end
                                    Foodstamps_C=max(SNAP[nk,site, q]-0.3*max(0.8*Earnings[site,q]*Work[w]+Benefit[nk,1,q]-134,0),0 )# I can probably do better than copy-pasting...

                                    Foodstamps_C=Eligible[e]*max(SNAP[nk,site, q]-0.3*max(0.8*Earnings[site,q]*Too_rich*Work[w]+(1-Too_rich)*Benefit[nk,1,q]-134,0),0 )+
                                                    (1-Eligible[e])*max(SNAP[nk,site, q]-0.3*max(0.8*Earnings[site,q]*Work[w]-134,0),0 ) # no benefit to subtract if ineligible

                                    CTJF=(Foodstamps_C+Benefit[nk,1,q])*(1-Too_rich)
                                    budget1[site,2,nk,a0,q,e,p,w]=Program[p]*(Eligible[e]*(CTJF)+(1-Eligible[e])*(Too_rich)*Foodstamps_C)+
                                                                            Earnings[site,q]*Work[w]
                                    budget1[site,3,nk,a0,q,e,p,w]=budget1[site,2,nk,a0,q,e,p,w]
                                    Foodstamps_receipt[site,2,nk,a0,q,e,p,w]=Program[p]*((1-Eligible[e])*(1-Too_rich)*Foodstamps_C+Eligible[e]*(Too_rich)*SNAP[nk,site, q])
                                    Foodstamps_receipt[site,3,nk,a0,q,e,p,w]=Foodstamps_receipt[site,2,nk,a0,q,e,p,w]
                                # FTP treatment
                                elseif site==2
                                Foodstamps_F=max(SNAP[nk,site, q]-0.3*max(0.8*Earnings[site,q]*Work[w]-134,0),0 )
                                FTP=max(Benefit[nk,site,q]-0.5*max(Earnings[site,q]*Work[w]-200,0),0)
                                budget1[site,2,nk,a0,q,e,p,w]=Program[p]*(FTP*Eligible[e]+Foodstamps_F)+
                                                                            Earnings[site,q]*Work[w]
                                budget1[site,3,nk,a0,q,e,p,w]=budget1[site,2,nk,a0,q,e,p,w]
                                Foodstamps_receipt[site,2,nk,a0,q,e,p,w]=Program[p]*Foodstamps_F
                                Foodstamps_receipt[site,3,nk,a0,q,e,p,w]=Foodstamps_receipt[site,2,nk,a0,q,e,p,w]
                            elseif site==3 || site==4
                            # MFIP treatment
                                Foodstamps_M=max(SNAP[nk,site, q]-0.3*max(0.8*Earnings[site,q]*Work[w]-134,0),0 )
                                D=Benefit[nk,site,q]+Foodstamps_M # notice the deviation from the formula in the document
                                MFIP=max( min(1.2*D-(1-0.38)*Earnings[site,q]*Work[w],D)  ,0)
                                budget1[site,2,nk,a0,q,e,p,w]=Program[p]*(MFIP*Eligible[e]+(1-Eligible[e])*Foodstamps_M)+
                                                                            Earnings[site,q]*Work[w]
                                budget1[site,3,nk,a0,q,e,p,w]=budget1[site,2,nk,a0,q,e,p,w]
                            end
                        end # end site loop

                    end
                end
            end
        end
    end
end
budget1[1,1,2,2,1,2,2]
#budget1=budget1.+0.0000001
#minimum(budget1)
#findmax(budget1[:,3,:,:,:,:,:,:])

Budget_Ageout=zeros(4,Dev_Years*4+1,2,2)
Foodstamps_receipt_ageout=zeros(4,Dev_Years*4+1,2,2)
@time @inbounds @simd for site in 1:4
    @inbounds @simd     for q in 1:(Dev_Years*4+1)
    @inbounds @simd         for w in 1:2
    @inbounds @simd             for p in 1:2
                        Budget_Ageout[site,q,p,w]=Earnings[site,q]*Work[w]+Program[p]*max(SNAP0[site,q]-0.3*max(0.8*Earnings[site,q]*Work[w]-134,0),0)

            end
        end
    end
end
Budget_Ageout
@time  Budget_Ageout[2,1,2,2]
#Budget_Ageout
#Budget_Ageout=Budget_Ageout.+0.000001


function AFDC(q, nk, earnings, eligible,participation,site; Pov_Guidelines=Poverty, Benefit=Benefit, SNAP=SNAP, SNAP0=SNAP0, ageout=0)

    Ben=0 # reminder that q is date
    FS=0
    if ageout==0
    FS=SNAP[nk,site, q]
    Ben=Benefit[nk,site,q]
    else
    FS=SNAP0[site,q]
    Ben=0
    end

    Welfare=participation*(eligible*max(Ben-(1-0.33)*max(earnings-120,0),0))*(1-ageout)
    FoodStamps=participation*(max(FS-0.3*max(0.8*earnings+Welfare-134,0),0))
    Budget=Welfare+FoodStamps+earnings

    return AFDC=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps)

end



function CTJF(q, nk, earnings, eligible,participation; Pov_Guidelines=Poverty, Benefit=Benefit, SNAP=SNAP,SNAP0=SNAP0, ageout=0)

    Ben=0
    FS=0
    if ageout==0
    FS=SNAP[nk,1, q]
    Ben=Benefit[nk,1,q]
    else
    FS=SNAP0[1,q]
    Ben=0
    end

    # First I create a 0-1 variable for the income cutoff
    TooRich=0
    if earnings>Pov_Guidelines[nk, 1,q]
        TooRich=1
    end


    FoodStamps=participation*(eligible*max(FS-0.3*max(0.8*earnings*TooRich+(1-TooRich)*Ben-134,0),0 )+
                    (1-eligible)*max(FS-0.3*max(0.8*earnings-134,0),0 ))

    Welfare=Ben*participation*(1-TooRich)*eligible

    Budget=earnings+Welfare+FoodStamps

    return CTJF=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps)

end



function FTP(q, nk, earnings, eligible,participation; Benefit=Benefit, SNAP=SNAP,SNAP0=SNAP0, ageout=0)
    Ben=0
    FS=0
    if ageout==0
    FS=SNAP[nk,2, q]
    Ben=Benefit[nk,2,q]
    else
    FS=SNAP0[2,q]
    Ben=0
    end

    FoodStamps=participation*max(FS-0.3*max(0.8*earnings-134,0),0 )
    Welfare=participation*eligible*max(Ben-0.5*max(earnings-200,0),0)
    Budget=earnings+Welfare+FoodStamps

    return FTP=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps)
end

function MFIP(q, nk, earnings, eligible,participation,s; Benefit=Benefit, SNAP=SNAP,SNAP0=SNAP0, ageout=0)
    Ben=0
    FS=0
    if ageout==0
    FS=SNAP[nk,s, q]
    Ben=Benefit[nk,s,q]
    else
    FS=SNAP0[s,q]
    Ben=0
    end


    FoodStamps=participation*max(FS-0.3*max(0.8*earnings-134,0),0 )
    D=participation*Ben+FoodStamps # notice the deviation from the formula in the document
    MFIP=max( min(1.2*D-(1-0.38)*earnings,D)  ,0)
    Welfare=(MFIP-FoodStamps)*eligible
    Budget=earnings+Welfare+FoodStamps

    return MFIP=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps)
end


function Budget_Function(s,tr,nk,age0,q,e,p,h, earnings; ageout=0)

    # initialize
    Budget=0
    Welfare=0
    FoodStamps=0
    numkids=nk


    eligible=1
    if e==1
        eligible=0
    end
    if ageout==1
        eligible=0
    end
    participation=1
    if p==1
        participation=0
    end

    #work=1
    if h==1
        earnings=0
    end

    if tr==1
            A=AFDC(q, numkids, earnings,eligible,participation, s; ageout=ageout)
            Budget=A[1]
            Welfare=A[2]
            FoodStamps=A[3]
    else
        if s==1
            C=CTJF(q, numkids, earnings,eligible,participation; ageout=ageout)
            Budget=C.Budget
            Welfare=C.Welfare
            FoodStamps=C.FoodStamps
        elseif s==2
            F=FTP(q, numkids, earnings,eligible,participation; ageout=ageout)
            Budget=F.Budget
            Welfare=F.Welfare
            FoodStamps=F.FoodStamps
        elseif s==3 || s==4
            M=MFIP(q, numkids, earnings,eligible,participation,s; ageout=ageout)
            Budget=M.Budget
            Welfare=M.Welfare
            FoodStamps=M.FoodStamps

        end
    end

    Both=Welfare+FoodStamps




    return Output=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps, Both=Both)

end






budget2=zeros(4,3,3,Dev_Years,Dev_Years*4,2,2,2)
Foodstamps_receipt2=zeros(4,3,3,Dev_Years,Dev_Years*4,2,2,2)
@time @inbounds @simd for a0 in 1:Dev_Years # outermost loop is for child initial age, 0-17 (shifted by 1-indexing)
    @inbounds @simd    for q in 1:(Dev_Years*4) # next loop is for quarter
    @inbounds @simd         for nk in 1:3 # num kids
    @inbounds @simd             for e in 1:2 # eligible or not
    @inbounds @simd                 for w in 1:2 # working or not
    @inbounds @simd                     for p in 1:2 # program, aka participating or not
    @inbounds @simd                         for site in 1:4 # loop over controls
                                                A=AFDC(q, nk, Earnings[site,q]*Work[w], Eligible[e],Program[p],site)
                                                budget2[site,1,nk,a0,q,e,p,w]=A.Budget
                                                Foodstamps_receipt2[site,1,nk,a0,q,e,p,w]=A.FoodStamps

                                # CTJF treatment
                                if site==1 # next I fill the treatment arms in one by one
                                    C=CTJF(q, nk, Earnings[site,q]*Work[w], Eligible[e],Program[p])
                                    budget2[site,2,nk,a0,q,e,p,w]=C.Budget
                                    budget2[site,3,nk,a0,q,e,p,w]=C.Budget
                                    Foodstamps_receipt2[site,2,nk,a0,q,e,p,w]=C.FoodStamps
                                    Foodstamps_receipt2[site,3,nk,a0,q,e,p,w]=C.FoodStamps
                                # FTP treatment
                                elseif site==2
                                    F=FTP(q, nk, Earnings[site,q]*Work[w],Eligible[e],Program[p])

                                budget2[site,2,nk,a0,q,e,p,w]=F.Budget
                                budget2[site,3,nk,a0,q,e,p,w]=F.Budget
                                Foodstamps_receipt2[site,2,nk,a0,q,e,p,w]=F.FoodStamps
                                Foodstamps_receipt2[site,3,nk,a0,q,e,p,w]=F.FoodStamps
                            elseif site==3 || site==4
                            # MFIP treatment
                                M=MFIP(q, nk, Earnings[site,q]*Work[w],Eligible[e],Program[p],site)

                                budget2[site,2,nk,a0,q,e,p,w]=M.Budget
                                budget2[site,3,nk,a0,q,e,p,w]=M.Budget
                            end
                        end # end site loop

                    end
                end
            end
        end
    end
end
budget2






Budget_Ageout2=zeros(4,Dev_Years*4+1,2,2)
Foodstamps_receipt_ageout=zeros(4,Dev_Years*4+1,2,2)
@time @inbounds @simd for site in 1:4
    @inbounds @simd     for q in 1:(Dev_Years*4+1)
    @inbounds @simd         for w in 1:2
    @inbounds @simd             for p in 1:2
                        Budget_Ageout2[site,q,p,w]=AFDC(q, 1, Earnings[site,q]*Work[w],0,Program[p], site; ageout=1).Budget

            end
        end
    end
end


# next I test if budget works--does it reproduce the old budget function?

@time M1=Budget_Function(1,1, 2, 2,   1,2,2,1,Earnings[1,1])
@time M1.Budget

@time Budget_Function(2,1, 1, 19,  1,2,2,2,Earnings[2,1]; ageout=1).Budget

@time Budget_Function(3,1, 1, 19,   1,2,2,2,Earnings[2,1]; ageout=1).Budget

budget2[:,:,:,:,:,2,:,:].-budget1[:,:,:,:,:,2,:,:]

findmax(budget2[2,:,:,:,:,:,:,:].-budget1[2,:,:,:,:,:,:,:])
findmin(budget2[2,:,:,:,:,:,:,:].-budget1[2,:,:,:,:,:,:,:])

findmax(budget2[3,:,:,:,:,:,:,:].-budget1[3,:,:,:,:,:,:,:])
findmin(budget2[3,:,:,:,:,:,:,:].-budget1[3,:,:,:,:,:,:,:])

findmax(budget2[1,:,:,:,:,:,:,:].-budget1[1,:,:,:,:,:,:,:])
findmin(budget2[1,:,:,:,:,:,:,:].-budget1[1,:,:,:,:,:,:,:])

Budget_Ageout2.-Budget_Ageout

# end testing




budget1 *= 3
budget2 *= 3
Budget_Ageout *= 3
Budget_Ageout2 *= 3
Earnings *= 3
Foodstamps_receipt *=3

writedlm("budget",budget2)
writedlm("budget_ageout",Budget_Ageout2)
writedlm("earnings",Earnings)
writedlm("foodstamps_receipt",Foodstamps_receipt)
TimeLimit_Ind=[false true true; false true true; false false false]

TimeLimits=[0 7 7; 0 8 8; 0 0 0]

Work_Reqs_Ind=[false true true; false true true; false false true]
