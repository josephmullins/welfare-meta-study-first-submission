# budget matrix
cd("/Users/FilipB/github/welfare-meta-study/Code")

using CSV
using GLM
using DataFrames
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


budget1=zeros(3,3,3,17,17*4,2,2,2)

site=["CTJF" "FTP" "MFIP-LR-F"]


# I store the first year I observed a program. This will come in handy later.
First_Year=zeros(Int64,3)
for i in 1:3
    Control = Quarterly_Data[Quarterly_Data[:, :Site] .== site[i], :]
    First_Year[i]=(minimum(Control.Year))
end

Earnings=zeros(3,68)

for i in 1:3
    Control = Quarterly_Data[Quarterly_Data[:, :Site] .== site[i], :]
    Control = Control[Control[:, :Treatment] .== "Control", :]
    Control[!,:Monthly_Earnings].=Control[!,:Earnings]./(3 .*Control[!,:LFP]./100)
    Control=hcat(Control, axes(Control, 1))
    Control[!, :Period]=(Control[!,:x1])
    lm2 = lm(@formula(Monthly_Earnings ~ Period), Control)
    observed=length(Control[!,:Monthly_Earnings])
    for j in 1:observed
        Earnings[i,j]=Control[j,:Monthly_Earnings]
    end
    for j in (observed+1):68
        Earnings[i,j]=coef(lm2)[1]+j*coef(lm2)[2]
    end
end

Earnings

#=

SNAP payouts next

Some of this will be a little ugly:
CTJF started 2 years after the others. Hence, I need to shift some years by 2.

=#


SNAP=zeros(3,3,68) # first index is the number of kids, second is site, third is quarter
for i in 1:3
    SNAPRulesC = SNAPRules[SNAPRules[:, :NumChild] .== i, :]
    for j in 1:3
        SNAPRulesC2 = SNAPRulesC[SNAPRulesC[:, :year] .>= First_Year[j], :]
            for k in 1:17
                for z in 1:4
                SNAP[i,j,(4*k+z-4)]=SNAPRulesC2.MA[k]
                end
            end
    end
end
SNAP




#=

Poverty Cutoff

=#

Poverty=zeros(3,3,68)# first index is the number of kids, second is site, third is quarter for compatibility
PovGuideline
for j in 1:3
    #PovGuideline2=PovGuideline[:,(i+1)]
    PovGuideline2=PovGuideline[PovGuideline[:, :year] .>= First_Year[j], :] # but we do site first
    for i in 1:3
        for k in 1:17
            for z in 1:4
            Poverty[i,j,(4*k+z-4)]=PovGuideline2[k,(i+1)]
            end
        end
    end
end


#=

Benefits Guideline

=#

State=["Connecticut" "Florida" "Minnesota"]

BenStd
BS2=BenStd[BenStd[:, :NumChild] .== 1, :]
BS2=BS2[BS2[:, :state] .== State[1], :]
select!(BS2, Not(:NumChild))
BS3=melt(BS2,:state)
sort!(BS3, :variable)
BS3[!,:YearString]=string.(BS3.variable)
BS3[!,:Year]=parse.(Int64, BS3.YearString)
BS3=BS3[BS3[:, :Year] .>= 1994, :]
BS3.value[1]

Benefit=zeros(3,3,68)
for i in 1:3
    BS2=BenStd[BenStd[:, :NumChild] .== i, :]
    for j in 1:3
        BS3=BS2[BS2[:, :state] .== State[j], :]
        select!(BS3, Not(:NumChild))
        BS3=melt(BS3,:state)
        sort!(BS3, :variable)
        BS3[!,:YearString]=string.(BS3.variable)
        BS3[!,:Year]=parse.(Int64, BS3.YearString)
        BS3=BS3[BS3[:, :Year] .>= 1994, :] # fix this for CT later on!!!
        for k in 1:17
            for z in 1:4
            Benefit[i,j,(4*k+z-4)]=BS3.value[k]
            end
        end
    end
end



# index 1 indicates false, 2 true
Eligible=[0, 1]
Work=[0,1]
Program=[0,1]

Disregard=zeros(3,3) # one fixed disregard for each site and arm
Disregard[:,:].=120 # first arm is AFDC
Disregard[2,2].=200 # FTP disregard higher
Disregard[2,3].=200
# CTJF, unbounded disregard irrelevant by 1.0 variable disregard



MTR= zeros(3,3) # this is the percent disregard, sort of like a marginal tax rate
MTR[:,1].=0.33 # 33 percent for AFDC, control arm
MTR[1,2]=1.0 # CTJF= disregard all
MTR[1,3]=1.0 # CTJF= disregard all
MTR[3,2]=0.38 # MFIP= 0.38
MTR[3,3]=0.38



Generosity=zeros(3,3,3,68) # nk, site, arm, quarter
Generosity[:,:,1,:].=Benefit#.+SNAP # for MFIP, FTP, CTJF control, it's the same as federal benefits
Generosity[:,:,2,:].=Benefit#.+SNAP
Generosity[:,:,3,:].=Benefit#.+SNAP


Generosity[:,1,2,:].=543 #.+SNAP[:,1,:] # for CT
Generosity[:,1,3,:].=543 #.+SNAP[:,1,:]

# budget(site,arm,NumChild,age0,quarter,program,work, eligible) (NS x NT x 3 x 17 x Q x 2 x 2 x 2)
budget1=zeros(3,3,3,17,17*4,2,2,2)
for a0 in 1:17 # outermost loop is for child initial age, 0-17 (shifted by 1-indexing)
    for q in 1:68 # next loop is for quarter
        for nk in 1:3 # num kids
            for e in 1:2 # eligible or not
                for w in 1:2 # working or not
                    for p in 1:2 # program, aka participating or not
                        for site in 1:3
                            for arm in 1:3
                            budget1[site,arm,nk,a0,q,p,w,e]=Earnings[site,q]*Work[w]+
                                                            Eligible[e]*Program[p]*
                                                            max( Generosity[nk,site,arm, q]+SNAP[nk, site, q]-
                                                               (Earnings[site,q]*Work[w]-Disregard[site,arm])*(1-MTR[site,arm])  ,0)


                            end
                        end
                    end
                end
            end
        end
    end
end
budget1[:,:,:,:,:,:,2,:]
