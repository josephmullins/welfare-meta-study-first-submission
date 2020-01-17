# budget matrix
cd("/Users/FilipB/github/welfare-meta-study/Code")
#cd("/Users/joseph/welfare-meta-study/Code")
using CSV
using GLM
using DataFrames
using DelimitedFiles



# Import CPI
CPI = CSV.read("../Data/CPIAUCSL.csv")
CPI[!,:year].=Dates.year.(CPI.DATE)
CPI[!,:Year].=Dates.year.(CPI.DATE)
Ninetyone=@where(CPI,:Year.==1991)


# Lump everything into annualized data
    # Data
Quarterly_Data=CSV.read("../Data/QuarterlyData_Extended.csv")
Annualized_Data=by(Quarterly_Data,
            [:Year,:Site,:Treatment_Num],
            [:LFP, :Earnings,:Welfare,:Receipt,:FoodStamps,:NEWWS,:Treatment_Num,:TotInc]=>
            x->(LFP=mean(x.LFP),  Earnings=mean(x.Earnings)*4, Welfare=mean(x.Welfare) ,
            Receipt=mean(x.Receipt)*4,FoodStamps=mean(x.FoodStamps),
            TotInc=mean(x.TotInc)*4,
            NEWWS=mean(x.NEWWS),  Treatment_Num=mean(x.Treatment_Num)  )      )
Annualized_Data.Treatment_Num

# but NEWSS was recorded annually
Annualized_Data[!,:Treatment].="Control"
@with Annualized_Data begin
    Annualized_Data.Earnings.=ifelse.((:NEWWS.==0.00), :Earnings, :Earnings./4)
    Annualized_Data.Receipt.=ifelse.((:NEWWS.==0.00), :Receipt, :Receipt./4)
    Annualized_Data.TotInc.=ifelse.((:NEWWS.==0.00), :TotInc, :TotInc./4)

    Annualized_Data.Treatment.=ifelse.((:Treatment_Num.!=0), "Treatment", "Control")

end


Annualized_Data[!,:Imputed_Earnings].=Annualized_Data.Earnings./(0.01.*Annualized_Data.LFP)
Annualized_Data=join(Annualized_Data,CPI, on=:Year)
Annualized_Data[!,:CPI].=Annualized_Data.CPIAUCSL./Ninetyone.CPIAUCSL[1]
CSV.write("../Data/Annualized_Data.csv",Annualized_Data)
    # Moments
#Q_moms = CSV.read("../Data/QuarterlyMoms.csv")

#Q_moms=by(Q_moms,
#                [:Year,:Site,:Treatment],
#                [:LFP, :Participation,:Receipt,:TotInc]=>
#                x->(LFP=mean(x.LFP),  Participation=mean(x.Participation), Receipt=mean(x.Receipt)*4 ,
#                TotInc=mean(x.TotInc)*4 )      )

A_moms=@select(Annualized_Data,:Year,:Site,:Treatment_Num,:LFP,:Earnings,:Welfare,:Receipt,:FoodStamps,:CPI,:TotInc)
A_moms=rename(A_moms,:Treatment_Num=>:Treatment, :Welfare=>:Participation)
A_moms=unique(A_moms)
CSV.write("../Data/Annualized_Moments.csv",A_moms)


# Guidelines
BenStd=CSV.read("../Data/WelfareRules/BenStd.csv")
PovGuideline=CSV.read("../Data/WelfareRules/PovGuideline.csv")
SNAPRules=CSV.read("../Data/WelfareRules/SNAPRules.csv")
CaliRules=CSV.read("../Data/WelfareRules/CaliforniaBenefitRules.csv")




println("Checkpoint 1")

#PovGuideline2=join(PovGuideline,CPI,on=:year, kind=:left)
#SNAPRules2=join(SNAPRules,CPI,on=:year, kind=:left)

# see later for deflating benefits

#=
I observe the following rules:

    - 1 indexes the control, 2 the treatment
    - 1 indexes Connecticut, 2 Florida, 3 Minneapolis
    - 1 indexes ineligible, 2 eligible

    There are 17*4=68 quarters

=#


Dev_Years=18

site=["CTJF" "FTP" "MFIP-LR" "MFIP-RA" "LA-GAIN" "LFA-Atlanta" "LFA-GR"]
Program_Lengths=[4,5,4,4,2,2,2] # note: needs to be the number of CALENDAR YEARS

Control = Annualized_Data[Annualized_Data[:, :Site] .== site[1], :]
First_Year=(minimum(Control.Year))


# I store the first year I observed a program. This will come in handy later.
First_Year=zeros(Int64,7)
for i in 1:7
    Control = Annualized_Data[Annualized_Data[:, :Site] .== site[i], :]
    First_Year[i]=(minimum(Control.Year))
end

Earnings=zeros(7,Dev_Years+1) # I add some extra years for the parents whose kids will age out

for i in 1:7
    Control = Annualized_Data[Annualized_Data[:, :Site] .== site[i], :]
    Control = Control[Control[:, :Treatment] .== "Control", :]
    Control.Annual_Earnings = Control.Earnings./(Control.LFP/100)
    Control.Period = 1:size(Control)[1]
    lm2 = lm(@formula(Annual_Earnings ~ Period), Control)
    observed=length(Control.Annual_Earnings)
    for j in 1:observed
        Earnings[i,j]=Control[j,:Annual_Earnings]
    end
    for j in (observed+1):length(Earnings[1,:])
        Earnings[i,j]=coef(lm2)[1]+j*coef(lm2)[2]
    end
end


#=

SNAP payouts next
These start monthly so I annualize them

=#


SNAP=zeros(3,7,Dev_Years) # first index is the number of kids, second is site, third is quarter
for i in 1:3
    SNAPRulesC = SNAPRules[SNAPRules[:, :NumChild] .== i, :]
    for j in 1:7
        SNAPRulesC2 = SNAPRulesC[SNAPRulesC[:, :year] .>= First_Year[j], :]
            for k in 1:Dev_Years
                SNAP[i,j,k]=SNAPRulesC2.MA[k]*12
            end
    end
end


# need to consider parents whose kids aged out
SnapRules_0=SNAPRules[SNAPRules[:, :NumChild] .== 0, :]
SNAP0=zeros(7,Dev_Years+1)
for i in 1:7
    SnapRules_02 = SnapRules_0[SnapRules_0[:, :year] .>= First_Year[i], :]
            for k in 1:19
                    SNAP0[i,k]=SnapRules_02.MA[k]*12
            end
end


#=

Poverty Cutoff--this was already annual

=#

Poverty=zeros(3,7,Dev_Years+1)# first index is the number of kids, second is site, third is quarter for compatibility
PovGuideline
for j in 1:7
    #PovGuideline2=PovGuideline[:,(i+1)]
    PovGuideline2=PovGuideline[PovGuideline[:, :year] .>= First_Year[j], :] # but we do site first
    for i in 1:3
        for k in 1:Dev_Years
            #for z in 1:4
            Poverty[i,j,k]=PovGuideline2[k,(i+2)] # the index of this spreadsheet is family size, so 2=1+1kid
            #end
        end
    end
end


#=

Benefits Guideline

=#
println("Checkpoint 2")
State=["Connecticut" "Florida" "Minnesota" "Minnesota" "California" "Georgia" "Michigan"]


BenStd
BS2=BenStd[BenStd[:, :NumChild] .== 1, :]
BS2=BS2[BS2[:, :state] .== State[5], :]
deletecols!(BS2, :NumChild)
BS3=melt(BS2,:state)
sort!(BS3, :variable)
BS3[!,:YearString]=string.(BS3.variable)
BS3[!,:Year]=parse.(Int64, BS3.YearString)
BS3=BS3[BS3[:, :Year] .>=  First_Year[5], :]
Maxyear=maximum(BS3.Year)
L=length(BS3.Year)
BS3.value[L]
BS3.value[1]
Price_index=CPI[CPI[:, :Year] .== BS3.Year[1], :].CPIAUCSL[1]

Benefit=zeros(3,7,Dev_Years)

for i in 1:3
    BS2=BenStd[BenStd[:, :NumChild] .== i, :]
    for j in 1:7
        BS3=BS2[BS2[:, :state] .== State[j], :]
        deletecols!(BS3,:NumChild)

        BS3=melt(BS3,:state)
        sort!(BS3, :variable)

        BS3.YearString=string.(BS3.variable)
        BS3.Year=parse.(Int64, BS3.YearString)

        BS3=BS3[BS3[:, :Year] .>= First_Year[j], :] # fix this for CT later on!!!

        Maxyear=maximum(BS3.Year)
        L=length(BS3.Year)

        for k in 1:Dev_Years # remember I don't have enough years for ctjf
            for z in 1:4
                if k<=L
                Price_index=CPI[CPI[:, :Year] .== BS3.Year[k], :].CPIAUCSL[1]
                Benefit[i,j,k]=BS3.value[k].*12
                #Benefit[i,j,k]=BS3.value[k].*12./Price_index
                else
                Price_index=CPI[CPI[:, :Year] .== BS3.Year[L], :].CPIAUCSL[1]
                Benefit[i,j,k]=BS3.value[L].*12
                end
            end
        end

    end
end


typeof(Benefit)

# index 1 indicates false, 2 true
Eligible=[0, 1]
Work=[0,1]
Program=[0,1]





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

    Error=Budget-earnings

    return AFDC=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps,Error=Error)

end





function AFDC_MN(q, nk, earnings, eligible,participation,site; Pov_Guidelines=Poverty, Benefit=Benefit, SNAP=SNAP, SNAP0=SNAP0, ageout=0)

    Ben=0 # reminder that q is date
    FS=0
    if ageout==0
    FS=SNAP[nk,site, q]
    Ben=Benefit[nk,site,q]
    else
    FS=SNAP0[site,q]
    Ben=0
    end

    Welfare=participation*(eligible*max((Ben-FS)-(1-0.33)*max(earnings-120,0),0))*(1-ageout)
    FoodStamps=participation*(max(FS-0.3*max(0.8*earnings+Welfare-134,0),0))
    Budget=Welfare+FoodStamps+earnings

    Error=Budget-earnings

    return AFDC=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps,Error=Error)

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


    #FoodStamps=participation*(eligible*max(FS-0.3*max((earnings)*TooRich+Ben-134,0),0 )+(1-eligible)*max(FS-0.3*max(0.8*earnings-134,0),0 ))

    FoodStamps=participation*(eligible*max(FS-0.3*max(TooRich*earnings+Ben-134,0),0 )+
                    (1-eligible)*max(FS-0.3*max(0.8*earnings-134,0),0 ))

    Welfare=Ben*participation*(1-TooRich)*eligible

    Budget=earnings+Welfare+FoodStamps

    return CTJF=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps, earn=earnings,pov=Pov_Guidelines[nk, 1,q])

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


    #FoodStamps=participation*max(FS-0.3*max(0.8*earnings-134,0),0 )
    #D=participation*Ben+FoodStamps # notice the deviation from the formula in the document
    #MFIP=max( min(1.2*D-(1-0.38)*earnings,D)  ,0)
    #Welfare=(MFIP-FoodStamps)*eligible
    #Budget=earnings+Welfare+FoodStamps

    Welfare=eligible*participation*
                max(min(1.2*Ben-(1-0.38)*earnings,Ben),0)
    Budget=earnings+Welfare
    FoodStamps=0

    return MFIP=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps)
end

#function LA_Gain

#end



# Connecticut budget function
Program_Lengths
CTJF_Budget=zeros(2,4,4,2,2) # arms, years, kids from 0 to 3, participation choice, work choice
CTJF_Budget_TL=zeros(2,4,4,2,2)


working=[0,1]
participating=[0,1]

for years in 1:Program_Lengths[1]
    for kids in 1:3
        for p in 1:2
            for w in 1:2
                        CTJF_Budget[1,years,kids,p, w]=0
                        CTJF_Budget[2,years,kids,p, w]=0
                        # =CTJF(years, kids, Earnings[1]*working[w], participating[p])
            end
        end
    end
end

1+1










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
            Budget=A.Budget
            Welfare=A.Welfare
            FoodStamps=A.FoodStamps
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



sites_considered=[1 2 3 4 6 7]


budget2=zeros(7,3,3,Dev_Years,Dev_Years,2,2,2) # site, arm, nk, a0, year/quarter, eligible, participating, working
Foodstamps_receipt2=zeros(7,3,3,Dev_Years,Dev_Years,2,2,2) # see above
@time @inbounds @simd for a0 in 1:Dev_Years # outermost loop is for child initial age, 0-17 (shifted by 1-indexing)
    @inbounds @simd    for q in 1:(Dev_Years) # next loop is for year (was quarter)
    @inbounds @simd         for nk in 1:3 # num kids
    @inbounds @simd             for e in 1:2 # eligible or not
    @inbounds @simd                 for w in 1:2 # working or not
    @inbounds @simd                     for p in 1:2 # program, aka participating or not
    @inbounds @simd                         for site in sites_considered # loop over controls
                                                if site==1 || site==2
                                                A=AFDC(q, nk, Earnings[site,q]*Work[w], Eligible[e],Program[p],site)
                                                budget2[site,1,nk,a0,q,e,p,w]=A.Budget
                                                Foodstamps_receipt2[site,1,nk,a0,q,e,p,w]=deepcopy(A.FoodStamps)
                                                else
                                                    A=AFDC_MN(q, nk, Earnings[site,q]*Work[w], Eligible[e],Program[p],site)
                                                    budget2[site,1,nk,a0,q,e,p,w]=A.Budget
                                                    Foodstamps_receipt2[site,1,nk,a0,q,e,p,w]=deepcopy(A.FoodStamps)
                                                end
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

                            else
                                A=AFDC(q, nk, Earnings[site,q]*Work[w], Eligible[e],Program[p],site)
                                budget2[site,1,nk,a0,q,e,p,w]=A.Budget
                                Foodstamps_receipt2[site,1,nk,a0,q,e,p,w]=deepcopy(A.FoodStamps)

                            end
                        end # end site loop

                    end
                end
            end
        end
    end
end
budget2






Budget_Ageout2=zeros(7,Dev_Years+1,2,2)

@time @inbounds @simd for site in sites_considered
    @inbounds @simd     for q in 1:(Dev_Years+1) # this is year now
    @inbounds @simd         for w in 1:2
    @inbounds @simd             for p in 1:2
                        Budget_Ageout2[site,q,p,w]=AFDC(q, 1, Earnings[site,q]*Work[w],0,Program[p], site; ageout=1).Budget

            end
        end
    end
end


# Everything is annualized now

#budget2 *= 3
#Budget_Ageout2 *= 3
#Earnings *=3
#Foodstamps_receipt2 *=3

minimum(Earnings)
minimum(budget2)
minimum(Budget_Ageout2)
minimum(Foodstamps_receipt2)


writedlm("budget",budget2)
writedlm("budget_ageout",Budget_Ageout2)
writedlm("earnings",Earnings)
writedlm("foodstamps_receipt",Foodstamps_receipt2)


TimeLimit_Ind=[false true true; false true true; false false false]


#TimeLimits=[0 7 7; 0 8 8; 0 0 0]
# OK, so we do lose some accuracy by going to years



Work_Reqs_Ind=[false true true; false true true; false false true]
