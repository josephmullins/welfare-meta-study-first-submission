using Revise
include("BaselineModel.jl")
using DelimitedFiles


Earnings = reshape(readdlm("earnings"),4,73)
budget = reshape(readdlm("budget"),4,3,3,18,18*4,2,2,2)
Budget_Ageout = reshape(readdlm("budget_ageout"),4,73,2,2)
Foodstamps_receipt=reshape(readdlm("foodstamps_receipt"),4,3,3,18,18*4,2,2,2)
#findmax(budget1[:,3,:,:,:,:,:,:])

TimeLimit_Ind=[false true true; false true true; false false false; false false false]

TimeLimits=[0 7 7; 0 8 8; 0 0 0; 0 0 0]

Work_Reqs_Ind=[false true true; false true true; false false true;false false true];

τ=ones(4,3).*0.5

@time Mod1=initialize_model()

Mod1.αc
UpdateSpecificParams!(Mod1; αc=0.3)
Mod1.αc
UpdateSpecificParams!(Mod1)
Mod1.αc

Mod2.work_prob[3,:,:,:,:,:,:]


Mod2.work_prob[2,:,:,:,:,:,:]


@time Mod2=initialize_model()



Mod2.αWR




#=


Below are a few sanity checks

=#



@time S1_c=Simulate(Mod2,10000,30,1,1)
@time S2_c=Simulate(Mod2,10000,30,2,1)
@time S3_c=Simulate(Mod2,10000,30,3,1)

S1_c.Participation
mean(S1_c.Participation)
mean(S2_c.Participation)
mean(S3_c.Participation)


mean(S1_c.Childcare)
mean(S2_c.Childcare)
mean(S3_c.Childcare)

@time S1_t=Simulate(Mod2,10000,30,1,3)
@time S2_t=Simulate(Mod2,10000,30,2,3)
@time S3_t=Simulate(Mod2,10000,30,3,3)
@time S3_t2=Simulate(Mod2,10000,30,3,2)

budget1[3,3,1,3,25,2,2,2]
Earnings[3,25]

mean(S1_t.Participation)
mean(S2_t.Participation)
mean(S3_t.Participation)


mean(S1_t.Childcare)
mean(S2_t.Childcare)
mean(S3_t.Childcare)

Mod2.welf_prob[1,1,:,:,:,:]

Mod2.welf_prob[1,2,:,:,:,:]


mean(S1_t.LFP)
mean(S1_t.Earned_Income)

findmax(budget1[:,3,:,:,:,:,:,:])

UpdateSpecificParams!(Mod2; pc=ones(length(Mod2.pc))*0.1,ϵ=ones(length(Mod2.ϵ))*0.1)

@time S1_t=Simulate(Mod2,100,30,1,3)

# Preparing Data


Quarterly_Data=CSV.read("../Data/QuarterlyData.csv")
Quarterly_Data[!,:Site2].=Quarterly_Data[!,:Site]
Quarterly_Data[!,:Treatment2].=Quarterly_Data[!,:Treatment]
for i in 1:length(Quarterly_Data[!,:Site2])

    if Quarterly_Data[i,:Site]=="MFIP-LR-F" && Quarterly_Data[i,:Treatment]=="Control"
        Quarterly_Data[i,:Site2]="MFIP_LR"
        Quarterly_Data[i,:Treatment2]="Control"

    elseif Quarterly_Data[i,:Site]=="MFIP-LR-I" && Quarterly_Data[i,:Treatment]=="Control"
            Quarterly_Data[i,:Site2]="MFIP_LR"
            Quarterly_Data[i,:Treatment2]="Control"

    elseif Quarterly_Data[i,:Site]=="MFIP-LR-F" && Quarterly_Data[i,:Treatment]=="Treatment"
        Quarterly_Data[i,:Site2]="MFIP_LR"
        Quarterly_Data[i,:Treatment2]="Treatment_WR"


    elseif Quarterly_Data[i,:Site]=="MFIP-LR-I" && Quarterly_Data[i,:Treatment]=="Treatment"
        Quarterly_Data[i,:Site2]="MFIP_LR"
        Quarterly_Data[i,:Treatment2]="Treatment"

        # recent joiners next

    elseif Quarterly_Data[i,:Site]=="MFIP-R-F" && Quarterly_Data[i,:Treatment]=="Control"
        Quarterly_Data[i,:Site2]="MFIP"
        Quarterly_Data[i,:Treatment2]="Control"

    elseif Quarterly_Data[i,:Site]=="MFIP-R-I" && Quarterly_Data[i,:Treatment]=="Control"
        Quarterly_Data[i,:Site2]="MFIP"
        Quarterly_Data[i,:Treatment2]="Control"

    elseif Quarterly_Data[i,:Site]=="MFIP-R-F" && Quarterly_Data[i,:Treatment]=="Treatment"
        Quarterly_Data[i,:Site2]="MFIP"
        Quarterly_Data[i,:Treatment2]="Treatment_WR"


    elseif Quarterly_Data[i,:Site]=="MFIP-R-I" && Quarterly_Data[i,:Treatment]=="Treatment"
        Quarterly_Data[i,:Site2]="MFIP"
        Quarterly_Data[i,:Treatment2]="Treatment"
    end
end

Q2=delete!(Quarterly_Data, :Site)
Q2=unique(Q2)


#=

Now I start coding up the moment generating function

Quarters:

CTJF=4*4
Measure=3 yr

FTP=4*3+2
Measure=4 year (will take last 2 quarters)


MFIP_LR_NWR=3+4*2+1
MFIP_LR_WR=3+4*2+1

MFIP_R_NWR=3+4*2+1
MFIP_R_WR=3+4*2+1

=#
@time S3_t2=Simulate(Mod2,10000,30,4,2)
S3_t2.AGE





Child=CSV.read("../Data/ChildOutcomes.csv")
Child=Child[1:26,:]

Child[!,:Measurement].="Hi"
for i in 1:length(Child[:,1])
    Child[i,:Measurement]=string(Child[i,:Site],"+",Child[i,:AgeMin],"+",Child[i,:AgeMax])
end

Meas=Child[(Child[:,:Treatment].=="C"),:]

Child[!,:Achievement2].=1.0
for i in 1:length(Child[:,1])


    if Child[i,:Treatment]=="C"
         Child[i,:Treatment]=="Control"
    elseif Child[i,:Treatment]=="T"
        Child[i,:Treatment]=="Treatment"
    elseif Child[i,:Treatment]=="I"
        Child[i,:Treatment]=="Treatment"
    elseif Child[i,:Treatment]=="F"
        Child[i,:Treatment]=="Treatment_WR"
    end

    if Child[i,:Site]=="MFIP-A" && Child[i,:Treatment]=="C"
        Child[i,:Site]=="MFIP"
        Child[i,:Treatment]="Control"
    end

    for j in 1:length(Meas.Measurement)
        if Child[i,:Measurement]==Meas[j,:Measurement]
            Child[i,:Achievement2]=Child[i,:Achievement]-Meas[j,:Achievement]
        end
    end
end




#=

The trick with the age limits:

MFIP has different measurements across sites and arms,
so I need different age vectors at each site-arm combination

MN arm 3 is the "full" version with work reqs

=#



Lengths=[4*4,4*3+2,3+4*2+1,3+4*2+1]
N_arms=[2,2,3,3]
Site_IDs=[1,2,3,4]

Site_Names=["CT", "FL", "MN", "MN_LR"]
Moment_Names=["Participation", "Receipt", "LFP"]
Arm_Names=["Treatment", "Control", "Control_WR"]

Age_Minima=[   [[5,13], [5,13]],
               [[5,13], [5,13]],
               [[5,5,9,13],  [5],     [5,5,9,13]],
               [[5,5,9,13],  [5],     [5,5,9,13]]     ]


Age_Maxima=[   [[12,17], [12,17]],
               [[12,17], [12,17]],
               [[8,12,12,18],  [12],     [8,12,12,18]],
               [[8,12,12,18],  [12],     [8,12,12,18]]     ]

Age_Minima[1,1,1]
Age_Minima[1][1][1]

Quarter_Meas_Min=[4*2+1,4*3+1,4*2+1,4*2+1]
Quarter_Meas_Max=[4*3,4*3+2,4*3,4*3] #note that florida is a bit short


Age_Minima_CT[1]



#=

This next loop stores the data analogue to "Quarterlyvalues"

It stores data by site, arm, and moment

=#

function Get_Simulated_Moments(M::Model, Lengths, N_arms, Site_IDs, Age_Minima, Age_Maxima, Quarter_Meas_Max, Quarter_Meas_Min, Data)
    x = zeros(0) # stores quarterly
    data_moms=zeros(0)
    Site_Names=["CTJF", "FTP", "MFIP", "MFIP_LR"]
    Site_samples=[1000,1000,1000,1000]
    Moment_Names=["Welfare", "Receipt", "LFP"]
    Arm_Names=["Control", "Treatment", "Treatment_WR"]



    xr=zeros(0)

    S_name=String[]
    M_name=String[]
    A_name=String[]
    wts=zeros(0)

    y=zeros(0) #stores outcomes
    @time for i in 1:4
        QD=Data[Data[:,:Site2].==Site_Names[i],:]
        for j in 1:N_arms[i]
            QD2=QD[QD[:,:Treatment2].==Arm_Names[j],:]
            Sim1=Simulate(M,10000,Lengths[i],i,j)
            # Next I append the quarterly moments, in order
            for k in 1:Lengths[i] # participation
                append!(x,mean(Sim1.Participation[k,:]))
                push!(S_name,Site_Names[i])
                push!(M_name, Moment_Names[1])
                push!(A_name, Arm_Names[j])
                append!(xr,QD2.Welfare[k])
                append!(wts,Site_samples[i])
            end
            for k in 1:Lengths[i]
                append!(x,mean(Sim1.Benefit_Receipt[k,:]))
                push!(S_name,Site_Names[i])
                push!(M_name, Moment_Names[2])
                push!(A_name, Arm_Names[j])
                append!(xr,QD2.Receipt[k])
                append!(wts,Site_samples[i])
            end
            for k in 1:Lengths[i]
                append!(x,mean(Sim1.LFP[k,:]))
                push!(S_name,Site_Names[i])
                push!(M_name, Moment_Names[3])
                push!(A_name, Arm_Names[j])
                append!(xr,QD2.LFP[k])
                append!(wts,Site_samples[i])

            end

            # Now the child outcome moments
            for l in 1:length(Age_Minima[i][j])
                S2=Sim1.AGE[Quarter_Meas_Min[i]:Quarter_Meas_Max[i],:]
                S3=Sim1.Skills[Quarter_Meas_Min[i]:Quarter_Meas_Max[i],:]
                It=(S2.>=Age_Minima[i][j][l]) .& (S2.<=Age_Maxima[i][j][l])

                append!(y,mean(S3[It]))
            end


        end
    end



    Distance=sum((x.-xr).^2 .*wts)

    return (x=x,y=y, S_name=S_name,M_name=M_name,A_name=A_name,xr=xr, Distance=Distance)
end


SM=Get_Simulated_Moments(Mod2,Lengths,N_arms,Site_IDs, Age_Minima, Age_Maxima, Quarter_Meas_Max, Quarter_Meas_Min,Q2)


SM.Distance










C2=[C1[:,1]; C1[:,2]; C1[:,3]]

function Get_Real_Moments()

    Quarterly_Data

end


Quarterly_Data

CT=Quarterly_Data[(Quarterly_Data[:,:Site].=="CTJF"),:]

CTC=CT[(CT[:,:Treatment].=="Control"),:]

Simulate(Mod2,10000,4*4,1,1)

unique(Quarterly_Data)


Q2=Q2[Q2[:,:Site2].!="LFA-Atlanta",:]

show(Quarterly_Data.Site2)
Site_Names=["CTJF", "FTP", "MFIP", "MFIP_LR"]
Moment_Names=["Welfare", "Receipt", "LFP"]
Arm_Names=["Control", "Treatment", "Treatment_WR"]
QD=Q2[Q2[:,:Site2].==Site_Names[1],:]

QD2=QD[QD[:,:Treatment2].==Arm_Names[1],:]

QD2.Welfare
