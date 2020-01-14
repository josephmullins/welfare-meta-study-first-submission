cd("/Users/FilipB/github/welfare-meta-study/Code")

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
A_moms=@select(Annualized_Data,:Year,:Site,:Treatment_Num,:LFP,:Earnings,:Welfare,:Receipt,:FoodStamps,:CPI,:TotInc)
A_moms=rename(A_moms,:Treatment_Num=>:Treatment, :Welfare=>:Participation)
A_moms=unique(A_moms)
CSV.write("../Data/Annualized_Moments.csv",A_moms)


# Guidelines
BenStd=CSV.read("../Data/WelfareRules/BenStd.csv")
PovGuideline=CSV.read("../Data/WelfareRules/PovGuideline.csv")
SNAPRules=CSV.read("../Data/WelfareRules/SNAPRules.csv")
CaliRules=CSV.read("../Data/WelfareRules/CaliforniaBenefitRules.csv")
    CaliRules.Year=CaliRules.Year.-1995 # everywhere else, "Year" will be indexed as 1,2....
    CaliRules


println("Checkpoint 1")


# see later for deflating benefits

#=
I observe the following rules:

    - 1 indexes the control, 2 the treatment
    - 1 indexes Connecticut, 2 Florida, 3/4 Minneapolis, 5 LA, 6 ATL, 7 GR, 8 Riverside
    - 1 indexes ineligible, 2 eligible

    There are 17*4=68 quarters

=#


Dev_Years=18

site=["CTJF" "FTP" "MFIP-LR" "MFIP-RA" "LA-GAIN" "Atlanta"  "GR" "Riverside"]
cdata=@where(Annualized_Data, :Treatment_Num.==0)
cdata[!,:Lengths].=1
Program_Lengths1=by(cdata, :Site, :Lengths=>sum)
Program_Lengths1=Program_Lengths1.Lengths_sum # note: needs to be the number of CALENDAR YEARS
Program_Lengths=[4, 5, 4, 4, 3, 5, 5, 5,]

# I store the first year I observed a program. This will come in handy later.
First_Year=zeros(Int64,8)
for i in 1:8
    Control = Annualized_Data[Annualized_Data[:, :Site] .== site[i], :]
    First_Year[i]=(minimum(Control.Year))
end

# I store labor market earnings for all the relevant years and locations
# I no longer extrapolate
Earnings=zeros(8,Dev_Years+1) # I add some extra years for the parents whose kids will age out
for i in 1:8
    Control = Annualized_Data[Annualized_Data[:, :Site] .== site[i], :]
    Control = Control[Control[:, :Treatment] .== "Control", :]
    Control.Annual_Earnings = Control.Earnings./(Control.LFP/100)
    Control.Period = 1:size(Control)[1]
    #lm2 = lm(@formula(Annual_Earnings ~ Period), Control)
    observed=length(Control.Annual_Earnings)
    for j in 1:observed
        Earnings[i,j]=Control[j,:Annual_Earnings]
    end
    for j in (observed+1):length(Earnings[1,:])
        #Earnings[i,j]=coef(lm2)[1]+j*coef(lm2)[2]
        Earnings[i,j]=Control[observed,:Annual_Earnings]
    end
end


#=
SNAP payouts next
These start monthly so I annualize them
=#


SNAP=zeros(4,8,Dev_Years) # first index is the number of kids, second is site, third is quarter
for i in 1:4
    SNAPRulesC = SNAPRules[SNAPRules[:, :NumChild] .== i-1, :]
    for j in 1:8
        SNAPRulesC2 = SNAPRulesC[SNAPRulesC[:, :year] .>= First_Year[j], :]
            for k in 1:Dev_Years
                SNAP[i,j,k]=SNAPRulesC2.MA[k]*12
            end
    end
end



#=
Poverty Cutoff--this was already annual
=#

Poverty=zeros(4,8,Dev_Years+1)# first index is the number of kids, second is site, third is quarter for compatibility
PovGuideline
for j in 1:8
    PovGuideline2=PovGuideline[PovGuideline[:, :year] .>= First_Year[j], :] # but we do site first
    for i in 1:4
        for k in 1:Dev_Years+1
            Poverty[i,j,k]=PovGuideline2[k,(i+1)] # the index of this spreadsheet is family size, so 2=1+1kid
        end
    end
end
Poverty




#=
Benefits Guideline
If no children, then benefits are zero
=#
println("Checkpoint 2")
State=["Connecticut" "Florida" "Minnesota" "Minnesota" "California" "Georgia" "Michigan" "California"]


BenStd


Benefit=zeros(4,8,Dev_Years)
for i in 2:4 # benefit with 0 kids is 0
    BS2=BenStd[BenStd[:, :NumChild] .== i-1, :]
    for j in 1:8
        BS3=BS2[BS2[:, :state] .== State[j], :]
        deletecols!(BS3,:NumChild)

        BS3=stack(BS3,2:ncol(BS3))
        #BS3=melt(BS3,:state)

        sort!(BS3, :variable)

        BS3.YearString=string.(BS3.variable)
        BS3.Year=parse.(Int64, BS3.YearString)

        BS3=BS3[BS3[:, :Year] .>= First_Year[j], :] # fix this for CT later on!!!

        Maxyear=maximum(BS3.Year)
        L=length(BS3.Year) # how many years do we observe a given site?

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

Benefit

#=

 Budget functions!!

=#



function AFDC(q, nk, earnings,participation,site; Pov_Guidelines=Poverty, Benefit=Benefit, SNAP=SNAP)


        FS=SNAP[nk,site, q]
        Ben=Benefit[nk,site,q]


    Welfare=participation*(max(Ben-(1-0.33)*max(earnings-120,0),0))
    FoodStamps=participation*(max(FS-0.3*max(0.8*earnings+Welfare-134,0),0))
    Budget=Welfare+FoodStamps+earnings

    Benefits=Welfare+FoodStamps


    return AFDC=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps, Benefits=Benefits)

end

function AFDC_MN(q, nk, earnings,participation,site; Pov_Guidelines=Poverty, Benefit=Benefit, SNAP=SNAP)

        FS=SNAP[nk,site, q]
        Ben=Benefit[nk,site,q]

    #=
    Note: the welfare payout calculation assumes that the benefit standard in MN
    sweeps in AFDC and food stamps, so to calculate AFDC for a family earning
    nothing I must subtract FS from Ben
    =#

    Welfare=participation*(max((Ben-FS)-(1-0.33)*max(earnings-120,0),0))
    FoodStamps=participation*(max(FS-0.3*max(0.8*earnings+Welfare-134,0),0))
    Budget=Welfare+FoodStamps+earnings

    Benefits=Welfare+FoodStamps

    return AFDC=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps, Benefits=Benefits)

end



function CTJF(q, nk, earnings,participation; Pov_Guidelines=Poverty, Benefit=Benefit, SNAP=SNAP, eligible=1)


        FS=SNAP[nk,1, q]
        Ben=Benefit[nk,1,q]

    # First I create a 0-1 variable for the income cutoff
    TooRich=0
    if earnings>Pov_Guidelines[nk, 1,q]
        TooRich=1
    end

    FoodStamps=participation*(max(FS-0.3*max(TooRich*earnings+Ben-134,0),0 )+
                    (1-eligible)*max(FS-0.3*max(0.8*earnings-134,0),0 ))

    Welfare=Ben*participation*(1-TooRich)*eligible

    Budget=earnings+Welfare+FoodStamps

    Benefits=Welfare+FoodStamps

    return CTJF=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps,Benefits=Benefits)

end



function FTP(q, nk, earnings,participation; Benefit=Benefit, SNAP=SNAP, eligible=1)

        FS=SNAP[nk,2, q]
        Ben=Benefit[nk,2,q]

    FoodStamps=participation*max(FS-0.3*max(0.8*earnings-134,0),0 )
    Welfare=participation*eligible*max(Ben-0.5*max(earnings-200,0),0)
    Budget=earnings+Welfare+FoodStamps
    Benefits=Welfare+FoodStamps

    return FTP=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps, Benefits=Benefits)
end



function MFIP(q, nk, earnings, participation; Benefit=Benefit, SNAP=SNAP)

        FS=SNAP[nk,3, q]
        Ben=Benefit[nk,3,q]


    Welfare=participation*
                max(min(1.2*Ben-(1-0.38)*earnings,Ben),0)
    Budget=earnings+Welfare
    FoodStamps=0
    Benefits=Welfare+FoodStamps


    return MFIP=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps, Benefits=Benefits)
end

c1=@where(CaliRules,:Year.==1,:NumKids.==1)
c1.MaxBenefit[1]


function LA_GAIN(year, nk, earnings, participation; Benefit=Benefit, SNAP=SNAP, CaliRules=CaliRules)

    NST=0
    MaxBen=0

    FS=SNAP[nk,5,year]

    if nk>1
        cal=@where(CaliRules,:Year.==year,:NumKids.==nk-1)
        NST=cal.NeedStandard[1]
        MaxBen=cal.MaxBenefit[1]
    end

    Welfare=participation*max(  min(MaxBen, NST-(1-0.33)*max(earnings-120,0)),0)
    FoodStamps=participation*(max(FS-0.3*max(0.8*earnings+Welfare-134,0),0))
    Budget=Welfare+FoodStamps+earnings
    Benefits=Welfare+FoodStamps

    return LA_Gain=(Budget=Budget, Welfare=Welfare, FoodStamps=FoodStamps, Benefits=Benefits)

end

working=[0,1]
participating=[0,1]



# Connecticut budget matrices
Program_Lengths
CTJF_Budget=zeros(2,Program_Lengths[1],4,2,2) # arms, years, kids from 0 to 3, participation choice, work choice
CTJF_Budget_Ineligible=zeros(2,Program_Lengths[1],4,2,2) # for when the time limit hits
CTJF_Benefits=zeros(2,Program_Lengths[1],4,2,2) # just welfare benefits and food stamps
CTJF_Benefits_Ineligible=zeros(2,Program_Lengths[1],4,2,2) # just welfare benefits and food stamps


AFDC(1, 2, Earnings[1,1]*working[2], participating[2],1)

for years in 1:Program_Lengths[1]
    for kids in 1:4 # kids=1 implies no kids
        for p in 1:2
            for w in 1:2
                        site=1
                            # Arm 1 is control, arm 2 is treatment
                            Arm1_E=AFDC(years, kids, Earnings[site,years]*working[w], participating[p],site)
                            Arm1_I=AFDC(years, kids, Earnings[site,years]*working[w], participating[p],site)
                            Arm2_E=CTJF(years, kids, Earnings[site,years]*working[w], participating[p])
                            Arm2_I=CTJF(years, kids, Earnings[site,years]*working[w], participating[p], eligible=0)

                        CTJF_Budget[1,years,kids,p, w]=Arm1_E.Budget
                        CTJF_Budget[2,years,kids,p, w]=Arm2_E.Budget

                        CTJF_Budget_Ineligible[1,years,kids,p, w]=Arm1_I.FoodStamps
                        CTJF_Budget_Ineligible[2,years,kids,p, w]=Arm2_E.Budget

                        CTJF_Benefits[1,years,kids,p, w]=Arm1_E.Benefits
                        CTJF_Benefits[2,years,kids,p, w]=Arm2_E.Benefits

                        CTJF_Benefits_Ineligible[1,years,kids,p, w]=Arm1_I.Benefits
                        CTJF_Benefits_Ineligible[2,years,kids,p, w]=Arm2_I.Benefits

            end
        end
    end
end

writedlm("Budgets/CTJF_MAIN",CTJF_Budget)
writedlm("Budgets/CTJF_INELIGIBLE",CTJF_Budget_Ineligible)
writedlm("Budgets/CTJF_BENEFITS",CTJF_Benefits)
writedlm("Budgets/CTJF_BENEFITS_INELIGIBLE",CTJF_Benefits_Ineligible)


# Florida budget matrices

FL_Budget=zeros(2,Program_Lengths[2],4,2,2) # arms, years, kids from 0 to 3, participation choice, work choice
FL_Budget_Ineligible=zeros(2,Program_Lengths[2],4,2,2) # for when the time limit hits
FL_Benefits=zeros(2,Program_Lengths[2],4,2,2) # just welfare benefits and food stamps
FL_Benefits_Ineligible=zeros(2,Program_Lengths[2],4,2,2) # for when the time limit hits




for years in 1:Program_Lengths[2]
    for kids in 1:4 # kids=1 implies no kids
        for p in 1:2
            for w in 1:2
                        site=2
                            # Arm 1 is control, arm 2 is treatment
                            Arm1_E=AFDC(years, kids, Earnings[site,years]*working[w], participating[p],site)
                            Arm1_I=AFDC(years, kids, Earnings[site,years]*working[w], participating[p],site)
                            Arm2_E=FTP(years, kids, Earnings[site,years]*working[w], participating[p])
                            Arm2_I=FTP(years, kids, Earnings[site,years]*working[w], participating[p], eligible=0)

                        FL_Budget[1,years,kids,p, w]=Arm1_E.Budget
                        FL_Budget[2,years,kids,p, w]=Arm2_E.Budget

                        FL_Budget_Ineligible[1,years,kids,p, w]=Arm1_I.FoodStamps
                        FL_Budget_Ineligible[2,years,kids,p, w]=Arm2_E.Budget

                        FL_Benefits[1,years,kids,p, w]=Arm1_E.Benefits
                        FL_Benefits[2,years,kids,p, w]=Arm2_E.Benefits

                        FL_Benefits_Ineligible[1,years,kids,p, w]=Arm1_I.Benefits
                        FL_Benefits_Ineligible[2,years,kids,p, w]=Arm2_I.Benefits


            end
        end
    end
end

writedlm("Budgets/FTP_MAIN",FL_Budget)
writedlm("Budgets/FTP_INELIGIBLE",FL_Budget_Ineligible)
writedlm("Budgets/FTP_BENEFITS",FL_Benefits)
writedlm("Budgets/FTP_BENEFITS_INELIGIBLE",FL_Benefits_Ineligible)


# MN budget matrices


MN_LR_Budget=zeros(2,Program_Lengths[3],4,2,2) # arms, years, kids from 0 to 3, participation choice, work choice
MN_RA_Budget=zeros(2,Program_Lengths[3],4,2,2) # arms, years, kids from 0 to 3, participation choice, work choice
MN_LR_Benefits=zeros(2,Program_Lengths[4],4,2,2) # just welfare benefits and food stamps
MN_RA_Benefits=zeros(2,Program_Lengths[4],4,2,2) # just welfare benefits and food stamps

Program_Lengths
for years in 1:Program_Lengths[3]
    for kids in 1:4 # kids=1 implies no kids
        for p in 1:2
            for w in 1:2
                            # Arm 1 is control, arm 2 is treatment
                        Arm1=AFDC_MN(years, kids, Earnings[3,years]*working[w], participating[p],3)
                        Arm2=MFIP(years, kids, Earnings[3,years]*working[w], participating[p])

                        MN_LR_Budget[1,years,kids,p, w]=Arm1.Budget
                        MN_LR_Budget[2,years,kids,p, w]=Arm2.Budget
                        MN_LR_Benefits[1,years,kids,p, w]=Arm1.Benefits
                        MN_LR_Benefits[2,years,kids,p, w]=Arm2.Benefits

                        # Site 3 is LR, 4 is RA
                        Arm1=AFDC_MN(years, kids, Earnings[4,years]*working[w], participating[p],4)
                        Arm2=MFIP(years, kids, Earnings[4,years]*working[w], participating[p])

                    MN_RA_Budget[1,years,kids,p, w]=Arm1.Budget
                    MN_RA_Budget[2,years,kids,p, w]=Arm2.Budget
                    MN_RA_Benefits[1,years,kids,p, w]=Arm1.Benefits
                    MN_RA_Benefits[2,years,kids,p, w]=Arm2.Benefits


            end
        end
    end
end

writedlm("Budgets/MFIP_LR_MAIN",MN_LR_Budget)
writedlm("Budgets/MFIP_LR_BENEFITS",MN_LR_Benefits)

writedlm("Budgets/MFIP_RA_MAIN",MN_RA_Budget)
writedlm("Budgets/MFIP_RA_BENEFITS",MN_RA_Benefits)


# LA budget matrices

LA_Budget=zeros(2,Program_Lengths[5],4,2,2) # arms, years, kids from 0 to 3, participation choice, work choice
LA_Benefits=zeros(2,Program_Lengths[5],4,2,2)
Program_Lengths[5]

for years in 1:3#Program_Lengths[5]
    for kids in 1:4 # kids=1 implies no kids
        for p in 1:2
            for w in 1:2
                        site=5

                            # Arm 1 is control, arm 2 is treatment
                            Arm1=AFDC(years, kids, Earnings[site,years]*working[w], participating[p],site)
                            Arm2=LA_GAIN(years, kids, Earnings[site,years]*working[w], participating[p])

                        LA_Budget[1,years,kids,p, w]=Arm1.Budget
                        LA_Budget[2,years,kids,p, w]=Arm2.Budget


                        LA_Benefits[1,years,kids,p, w]=Arm1.Benefits
                        LA_Benefits[2,years,kids,p, w]=Arm2.Benefits



            end
        end
    end
end

writedlm("Budgets/LAGAIN_MAIN",LA_Budget)
writedlm("Budgets/LAGAIN_BENEFITS",LA_Benefits)


# Atlanta budget matrices


ATL_Budget=zeros(1,Program_Lengths[6],4,2,2) # arms, years, kids from 0 to 3, participation choice, work choice
ATL_Benefits=zeros(1,Program_Lengths[6],4,2,2) # just welfare benefits and food stamps


for years in 1:Program_Lengths[6]
    for kids in 1:4 # kids=1 implies no kids
        for p in 1:2
            for w in 1:2
                        site=6
                            # Arm 1 is control, arm 2 is treatment
                            Arm1=AFDC(years, kids, Earnings[site,years]*working[w], participating[p],site)


                        ATL_Budget[1,years,kids,p, w]=Arm1.Budget
                        ATL_Benefits[1,years,kids,p, w]=Arm1.Benefits

            end
        end
    end
end
ATL_Budget
ATL_Benefits

writedlm("Budgets/NEWWS_A_MAIN",ATL_Budget)
writedlm("Budgets/NEWWS_A_BENEFITS",ATL_Benefits)


# GR budget matrices


GR_Budget=zeros(1,Program_Lengths[7],4,2,2) # arms, years, kids from 0 to 3, participation choice, work choice
GR_Benefits=zeros(1,Program_Lengths[7],4,2,2) # just welfare benefits and food stamps

for years in 1:Program_Lengths[7]
    for kids in 1:4 # kids=1 implies no kids
        for p in 1:2
            for w in 1:2
                        site=7
                            # Arm 1 is control, arm 2 is treatment
                            Arm1=AFDC(years, kids, Earnings[site,years]*working[w], participating[p],site)


                        GR_Budget[1,years,kids,p, w]=Arm1.Budget
                        GR_Benefits[1,years,kids,p, w]=Arm1.Benefits

            end
        end
    end
end

writedlm("Budgets/NEWWS_G_MAIN",GR_Budget)
writedlm("Budgets/NEWWS_G_BENEFITS",GR_Benefits)

# Riverside budget matrices

Riverside_Budget=zeros(1,Program_Lengths[8],4,2,2) # arms, years, kids from 0 to 3, participation choice, work choice
Riverside_Benefits=zeros(1,Program_Lengths[8],4,2,2) # just welfare benefits and food stamps

for years in 1:Program_Lengths[8]
    for kids in 1:4 # kids=1 implies no kids
        for p in 1:2
            for w in 1:2
                        site=8
                            # Arm 1 is control, arm 2 is treatment
                            Arm1=AFDC(years, kids, Earnings[site,years]*working[w], participating[p],site)


                        Riverside_Budget[1,years,kids,p, w]=Arm1.Budget
                        Riverside_Benefits[1,years,kids,p, w]=Arm1.Benefits

            end
        end
    end
end

writedlm("Budgets/NEWWS_R_MAIN",Riverside_Budget)
writedlm("Budgets/NEWWS_R_BENEFITS",Riverside_Benefits)
Riverside_Budget
