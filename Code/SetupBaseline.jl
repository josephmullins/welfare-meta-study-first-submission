# this script sets up some of the baseline characteristics for estimation in the first stage
using CSV
using DataFrames

# ---- Set up Experiment Features
C = CSV.read("../Data/ChildCareMoms_estimated.csv")

num_sites = 8
site_list = [:CTJF,:FTP,:LAGAIN,:MFIPLR,:MFIPRA,:NEWWSA,:NEWWSG,:NEWWSR]
site_str = ["CTJF","FTP","LAGAIN","MFIP-LR","MFIP-RA","NEWWS-A","NEWWS-G","NEWWS-R"]
n_arms = [1,1,2,3,3,2,2,2] #note; use only control arms for CTJF and FTP due to time limits
work_reqs = [zeros(8) [0,0,1,1,1,1,1,1] zeros(8)]
years = [3,4,2,3,3,5,5,5]
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
    πA0 = [πA[i,1]*ones(3)/3;πA[i,2]*ones(3)/3;πA[i,3]*ones(11)/11]
    πk0 = [πK[i]*ones(2)/3;1-πK[i]]
    π0[i,:,:] = πK*πA'
end

# Prices:
Prices = zeros(num_sites,3)
for i=1:num_sites
    for a=0:n_arms[i]-1
        Prices[i,a+1] = C.Price[C.Site.==site_str[i] .& C.Arm.==a]
    end
end

site_features = (T = years,n_arms = n_arms,work_reqs = work_reqs,π0 = π0,prices = Prices)

# --- Set up budget function

# --- Set up moments and wghts