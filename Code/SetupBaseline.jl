# this script sets up some of the baseline characteristics for estimation in the first stage
using CSV
using DataFrames
using DelimitedFiles

# ---- Set up Experiment Features
C = CSV.read("../Data/ChildCareMoms_estimated.csv")

num_sites = 8
site_list = [:CTJF,:FTP,:LAGAIN,:MFIPLR,:MFIPRA,:NEWWSA,:NEWWSG,:NEWWSR]
site_str = ["CTJF","FTP","LAGAIN","MFIP-LR","MFIP-RA","NEWWS-A","NEWWS-G","NEWWS-R"]
n_arms = [1,1,2,3,3,2,2,2] #note; use only control arms for CTJF and FTP due to time limits
work_reqs = [zeros(8) [1,0,1,1,1,1,1,1] zeros(8)]
years = [4,5,3,4,4,5,5,5]
year_meas = [3,4,2,3,3,2,2,2]

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

site_features = (T = years,n_arms = n_arms,work_reqs = work_reqs,π0 = π0,prices = Prices, year_meas = year_meas)

# ------------ Set up budget function ------------------- #
# order of budget function is: A x T x NK x 2 x 2
ctjf = reshape(readdlm("Budgets/CTJF_MAIN"),2,4,4,2,2)/52
lagain = reshape(readdlm("Budgets/LAGAIN_MAIN"),2,3,4,2,2)/52
ftp = reshape(readdlm("Budgets/FTP_MAIN"),2,5,4,2,2)/52
mfiplr = reshape(readdlm("Budgets/MFIP_LR_MAIN"),2,4,4,2,2)/52
mfipra = reshape(readdlm("Budgets/MFIP_RA_MAIN"),2,4,4,2,2)/52
newwsa = reshape(readdlm("Budgets/NEWWS_A_MAIN"),1,5,4,2,2)/52
newwsg = reshape(readdlm("Budgets/NEWWS_G_MAIN"),1,5,4,2,2)/52
newwsr = reshape(readdlm("Budgets/NEWWS_R_MAIN"),1,5,4,2,2)/52
budget = (CTJF = ctjf, FTP = ftp, LAGAIN = lagain, MFIPLR = mfiplr, MFIPRA = mfipra, NEWWSA = newwsa, NEWWSG = newwsg, NEWWSR = newwsr)
num_sites = 8

# ------------ Set up moments and wghts ---------------- #
D = CSV.read("../Data/Annualized_Moments.csv")
# CTJF
ctjf = zeros(4*2+1,2)
w_ctjf = 100*ones(4*2+1,2)
for a=0:1
    d = D[(D.Site.=="CTJF") .& (D.Treatment.==a),:]
    c = C[(C.Site.=="CTJF") .& (C.Arm.==a),:]
    ctjf[:,a+1] = [d.Participation/100;d.LFP/100;c.UsePaid[1]]
end
# FTP
ftp = zeros(5*2+1,2)
w_ftp = 100*ones(5*2+1,2)
for a=0:1
    d = D[(D.Site.=="FTP") .& (D.Treatment.==a),:]
    c = C[(C.Site.=="FTP") .& (C.Arm.==a),:]
    ftp[:,a+1] = [d.Participation/100;d.LFP/100;c.UsePaid[1]]
end

# LAGAIN
LA = zeros(3*2+1,2)
w_LA = 100*ones(3*2+1,2)
for a=0:1
    d = D[(D.Site.=="LA-GAIN") .& (D.Treatment.==a),:]
    c = C[(C.Site.=="LAGAIN") .& (C.Arm.==a),:]
    LA[:,a+1] = [d.Participation/100;d.LFP/100;c.UsePaid[1]]
end

# MFIP-LR
mfiplr = zeros(4*2+1,3)
w_mfiplr = 100*ones(4*2+1,3)
for a=0:2
    d = D[(D.Site.=="MFIP-LR") .& (D.Treatment.==a),:]
    c = C[(C.Site.=="MFIP-LR") .& (C.Arm.==a),:]
    mfiplr[:,a+1] = [d.Participation/100;d.LFP/100;c.UsePaid[1]]
end

# MFIP-RA
mfipra = zeros(4*2+1,3)
w_mfipra = 100*ones(4*2+1,3)
for a=0:2
    d = D[(D.Site.=="MFIP-RA") .& (D.Treatment.==a),:]
    c = C[(C.Site.=="MFIP-RA") .& (C.Arm.==a),:]
    mfipra[:,a+1] = [d.Participation/100;d.LFP/100;c.UsePaid[1]]
end

# NEWWS-A
newwsa = zeros(5*2+1,2)
w_newwsa = 100*ones(5*2+1,2)
for a=0:1
    d = D[(D.Site.=="Atlanta") .& (D.Treatment.==a),:]
    c = C[(C.Site.=="NEWWS-A") .& (C.Arm.==a),:]
    newwsa[:,a+1] = [d.Participation/100;d.LFP/100;c.UsePaid[1]]
end

# NEWWS-G
newwsg = zeros(5*2+1,2)
w_newwsg = 100*ones(5*2+1,2)
for a=0:1
    d = D[(D.Site.=="GR") .& (D.Treatment.==a),:]
    c = C[(C.Site.=="NEWWS-G") .& (C.Arm.==a),:]
    newwsg[:,a+1] = [d.Participation/100;d.LFP/100;c.UsePaid[1]]
end

# NEWWS-R
newwsr = zeros(5*2+1,2)
w_newwsr = 100*ones(5*2+1,2)
for a=0:1
    d = D[(D.Site.=="Riverside") .& (D.Treatment.==a),:]
    c = C[(C.Site.=="NEWWS-R") .& (C.Arm.==a),:]
    newwsr[:,a+1] = [d.Participation/100;d.LFP/100;c.UsePaid[1]]
end

moments = (CTJF = ctjf, FTP = ftp, LAGAIN = LA, MFIPLR = mfiplr, MFIPRA = mfipra, NEWWSA = newwsa, NEWWSG = newwsg, NEWWSR = newwsr)
# later one we can let the relative sample sizes inform the weights if we want
wghts = (CTJF = w_ctjf, FTP = w_ftp, LAGAIN = w_LA, MFIPLR = w_mfiplr, MFIPRA = w_mfipra, NEWWSA = w_newwsa, NEWWSG = w_newwsg, NEWWSR = w_newwsr)
