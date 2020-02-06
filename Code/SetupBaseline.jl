# this script sets up some of the baseline characteristics for estimation in the first stage
using CSV
using DataFrames
using DelimitedFiles

# ---- Set up Experiment Features
C = CSV.read("../Data/ChildCareMoms_estimated.csv")
#C.UsePaid = C.UsePaid.*C.Emp/100
num_sites = 8
site_list = [:CTJF,:FTP,:LAGAIN,:MFIPLR,:MFIPRA,:NEWWSA,:NEWWSG,:NEWWSR]
site_str = ["CTJF","FTP","LAGAIN","MFIP-LR","MFIP-RA","NEWWS-A","NEWWS-G","NEWWS-R"]
site_str2 = ["CTJF","FTP","LA-GAIN","MFIP-LR","MFIP-RA","Atlanta","GR","Riverside"]

n_arms = [2,2,2,3,3,2,2,2] #note; use only control arms for CTJF and FTP due to time limits
work_reqs = [[0,1,0,0,0,0,0,0] [1,0,1,1,1,1,1,1] zeros(8)] #<- FTP has work reqs in both
years = [4,4,3,4,4,5,5,5]
year_meas = [3,4,2,3,3,2,2,2]
yb = [1996,1994,1996,1994,1994,1991,1991,1991]
time_limits = zeros(8,3); time_limits[1,2] = 1; time_limits[2,2] = 1;
TLlength = zeros(Int64,8,3); TLlength[1,2] = 2; TLlength[2,2] = 2;

# below are the distributions across age of youngest 0-2,3-5,6+
πA = [[37.5 25.4 37.1];
    [43. 27.2 29.8];
    [52.3 25. 22.7];
    [35.7 35.8 28.5];
    [42.9 25.1 31.9];
    [0. 49.7 50.3];
    [29. 26.4 44.6];
    [0. 45.2 54.8]]

# below is the fraction of households with 2 or fewer children
πK = [67.7, 68.2, 73.1, 64.6, 82.4, 69.9, 82.2, 70.6]

# now construct the joint initial distribution
π0 = zeros(8,3,17)

for i=1:8
    πA0 = [πA[i,1]*ones(3)/300;πA[i,2]*ones(3)/300;πA[i,3]*ones(11)/1100]
    πk0 = [πK[i]*ones(2)/200;1-πK[i]/100]
    π0[i,:,:] = πk0*πA0'
end

# Prices:
Prices = zeros(num_sites,3)
for i=1:num_sites
    for a=0:n_arms[i]-1
        Prices[i,a+1] = C.Price[(C.Site.==site_str[i]) .& (C.Arm.==a)][1]
    end
end

site_features = (T = years,n_arms = n_arms,work_reqs = work_reqs,π0 = π0,prices = Prices, year_meas = year_meas, yb = yb, time_limits = time_limits,TLlength = TLlength, site_list = site_list)

# ------------ Set up budget function ------------------- #
# order of budget function is: A x T x NK x 2 x 2
ctjf = reshape(readdlm("Budgets/CTJF_MAIN"),2,4,4,2,2)/4
ctjf_fs = reshape(readdlm("Budgets/CTJF_BENEFITS_INELIGIBLE"),2,4,4,2,2)/4 #<- this isn't quite what we want
ctjf_ben = reshape(readdlm("Budgets/CTJF_BENEFITS"),2,4,4,2,2)/4 #<- this isn't quite what we want
ctjf_inel = reshape(readdlm("Budgets/CTJF_INELIGIBLE"),2,4,4,2,2)/4
ctjf_inel[2,:,:,:,1] .= 0
lagain = reshape(readdlm("Budgets/LAGAIN_MAIN"),2,3,4,2,2)/4
ftp = reshape(readdlm("Budgets/FTP_MAIN"),2,5,4,2,2)/4
ftp_inel = reshape(readdlm("Budgets/FTP_INELIGIBLE"),2,5,4,2,2)/4
ftp_inel[2,:,:,:,1] .= 0
# fix mfip budgets
mfiplr = reshape(readdlm("Budgets/MFIP_LR_MAIN"),2,4,4,2,2)/4
mfipra = reshape(readdlm("Budgets/MFIP_RA_MAIN"),2,4,4,2,2)/4
newwsa = reshape(readdlm("Budgets/NEWWS_A_MAIN"),1,5,4,2,2)/4
newwsg = reshape(readdlm("Budgets/NEWWS_G_MAIN"),1,5,4,2,2)/4
newwsr = reshape(readdlm("Budgets/NEWWS_R_MAIN"),1,5,4,2,2)/4
budget = (CTJF = ctjf, FTP = ftp, LAGAIN = lagain, MFIPLR = mfiplr, MFIPRA = mfipra, NEWWSA = newwsa, NEWWSG = newwsg, NEWWSR = newwsr,FTP_I = ftp_inel,CTJF_I = ctjf_inel)
num_sites = length(site_list)

# ------------ Set up moments and wghts ---------------- #
sample_size=[4803, 1405+1410,15683,3208,6009,4433,4554,8322]
D = CSV.read("../Data/Annualized_Moments.csv")

moms_collect = []
for i=1:num_sites
    moms = zeros(2*years[i]+1,n_arms[i])
    se = zeros(2*years[i]+1,n_arms[i])
    println(site_list[i])
    for a=1:n_arms[i]
        d = D[(D.Site.==site_str2[i]) .& (D.Treatment.==(a-1)),:]
        c = C[(C.Site.==site_str[i]) .& (C.Arm.==(a-1)),:]
        moms[:,a] = [d.Participation[1:years[i]]/100; d.LFP[1:years[i]]/100; c.UsePaid[1]]
        se[:,a] = sqrt.(2*(moms[:,a].*(1 .- moms[:,a]))/sample_size[i])
        se[end,a] = (moms[end,a]*(1-moms[end,a])/c.N[1])
    end
    append!(moms_collect,[(moms=moms,se=se)])
end
moments = (;zip(site_list,moms_collect)...)

#
#
# moms_collect = []
# for i=1:8
#     sname = site_list[i]
#     moms = []
#     T = site_features.T[i]
#     println(sname)
#     for a = 1:site_features.n_arms[i]
#         d = D[(D.Site.==site_str2[i]) .& (D.Treatment.==(a-1)),:]
#         d = d[1:T,:]
#         c = C[(C.Site.==site_str[i]) .& (C.Arm.==(a-1)),:]
#         append!(moms,[(Part = d.Participation/100,LFP = d.LFP/100,Care = c.UsePaid[1],Inc = d.TotInc/52)])
#     end
#     append!(moms_collect,[moms])
# end
# data_moments = (;zip(site_list,moms_collect)...)
