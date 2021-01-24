#!python
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
import seaborn as sns
import argparse
from matplotlib.pyplot import figure
from random import seed
from random import randint
import sys
from numpy import genfromtxt
pd.set_option('display.max_rows', 10000) #some default pandas settings will chop data
pd.set_option('display.max_columns', 10000)
pd.set_option('display.width', 10000)
np.set_printoptions(threshold=sys.maxsize) #if printed out, matrix won't be truncated


### ====================================== ARGPARSE ==================================== ###
parser = argparse.ArgumentParser(description='Compare ct data with RNAfold')
parser.add_argument('--verbose','-v', action='store_true',
                    default=False,
                    help='be verbose')

parser.add_argument('--graphics','-g', action='store_true',
                    default=False,
                    help='be verbose')

parser.add_argument('--printout','-p', action='store_true',
                    default=False,
                    help='print values')

parser.add_argument('--test', '-t', action='store_true',
                    default=False,
                    help='run with test csv')

parser.add_argument('--label','-l', action='store',
                    default="mod_unblock_comp_m_n",
                    help='parse label to graph')

parser.add_argument('--output','-o', action='store_true',
                    default=False,
                    help='recalculate boot/csvs')

args=parser.parse_args()

### ==================================== NAMES & IMPORT ================================ ###

n_input = "n"
m_input = "m"
mu_input = "mu"

### mod means
n_mean_m = "comp_" + m_input + "_n_fin_n_mean.csv"
c_mean_m = "comp_" + m_input + "_c_fin_n_mean.csv"
p_mean_m = "comp_" + m_input + "_p_fin_n_mean.csv"
cp_mean_m = "comp_" + m_input + "_cp_fin_n_mean.csv"

### mod bootstrapping
n_boot_m = "comp_" + m_input + "_n_fin_n_boot.csv"
c_boot_m = "comp_" + m_input + "_c_fin_n_boot.csv"
p_boot_m = "comp_" + m_input + "_p_fin_n_boot.csv"
cp_boot_m = "comp_" + m_input + "_cp_fin_n_boot.csv"

### mod unblock means
n_mean_mu = "comp_" + mu_input + "_n_fin_n_mean.csv"
c_mean_mu = "comp_" + mu_input + "_c_fin_n_mean.csv"
p_mean_mu = "comp_" + mu_input + "_p_fin_n_mean.csv"
cp_mean_mu = "comp_" + mu_input + "_cp_fin_n_mean.csv"

### mod unblock bootstrapping
n_boot_mu = "comp_" + mu_input + "_n_fin_n_boot.csv"
c_boot_mu = "comp_" + mu_input + "_c_fin_n_boot.csv"
p_boot_mu = "comp_" + mu_input + "_p_fin_n_boot.csv"
cp_boot_mu = "comp_" + mu_input + "_cp_fin_n_boot.csv"

### nat means
n_mean_n = "comp_" + n_input + "_n_fin_n_mean.csv"
c_mean_n = "comp_" + n_input + "_c_fin_n_mean.csv"
p_mean_n = "comp_" + n_input + "_p_fin_n_mean.csv"
cp_mean_n = "comp_" + n_input + "_cp_fin_n_mean.csv"

### nat bootstrapping
n_boot_n = "comp_" + n_input + "_n_fin_n_boot.csv"
c_boot_n = "comp_" + n_input + "_c_fin_n_boot.csv"
p_boot_n = "comp_" + n_input + "_p_fin_n_boot.csv"
cp_boot_n = "comp_" + n_input + "_cp_fin_n_boot.csv"


### Import of data
data_n_mean_m = genfromtxt(n_mean_m, delimiter = ",")
data_c_mean_m = genfromtxt(c_mean_m, delimiter = ",")
data_p_mean_m = genfromtxt(p_mean_m, delimiter = ",")
data_cp_mean_m = genfromtxt(cp_mean_m, delimiter = ",")

data_n_boot_m = genfromtxt(n_boot_m, delimiter = ",")
data_c_boot_m = genfromtxt(c_boot_m, delimiter = ",")
data_p_boot_m = genfromtxt(p_boot_m, delimiter = ",")
data_cp_boot_m = genfromtxt(cp_boot_m, delimiter = ",")


data_n_mean_n = genfromtxt(n_mean_n, delimiter = ",")
data_c_mean_n = genfromtxt(c_mean_n, delimiter = ",")
data_p_mean_n = genfromtxt(p_mean_n, delimiter = ",")
data_cp_mean_n = genfromtxt(cp_mean_n, delimiter = ",")

data_n_boot_n = genfromtxt(n_boot_n, delimiter = ",")
data_c_boot_n = genfromtxt(c_boot_n, delimiter = ",")
data_p_boot_n = genfromtxt(p_boot_n, delimiter = ",")
data_cp_boot_n = genfromtxt(cp_boot_n, delimiter = ",")


data_n_mean_mu = genfromtxt(n_mean_mu, delimiter = ",")
data_c_mean_mu = genfromtxt(c_mean_mu, delimiter = ",")
data_p_mean_mu = genfromtxt(p_mean_mu, delimiter = ",")
data_cp_mean_mu = genfromtxt(cp_mean_mu, delimiter = ",")

data_n_boot_mu = genfromtxt(n_boot_mu, delimiter = ",")
data_c_boot_mu = genfromtxt(c_boot_mu, delimiter = ",")
data_p_boot_mu = genfromtxt(p_boot_mu, delimiter = ",")
data_cp_boot_mu = genfromtxt(cp_boot_mu, delimiter = ",")



### ==================================== VARIABLES ================================== ###

mean_m_mcc_fold_mfe = data_n_mean_m[16]
mean_m_mcc_drtr_mfe = data_n_mean_m[17]
mean_m_mcc_drtr_prob = data_n_mean_m[18]
mean_m_mcc_fold_bpp = data_n_mean_m[19]
mean_m_mcc_drtr_bpp = data_n_mean_m[20]

mean_m_mcc_fold_mfe = mean_m_mcc_fold_mfe[1]
mean_m_mcc_drtr_mfe = mean_m_mcc_drtr_mfe[1]
mean_m_mcc_drtr_prob = mean_m_mcc_drtr_prob[1]
mean_m_mcc_fold_bpp = mean_m_mcc_fold_bpp[1]
mean_m_mcc_drtr_bpp = mean_m_mcc_drtr_bpp[1]


mean_mu_mcc_fold_mfe = data_n_mean_mu[16]
mean_mu_mcc_drtr_mfe = data_n_mean_mu[17]
mean_mu_mcc_drtr_prob = data_n_mean_mu[18]
mean_mu_mcc_fold_bpp = data_n_mean_mu[19]
mean_mu_mcc_drtr_bpp = data_n_mean_mu[20]

mean_mu_mcc_fold_mfe = mean_mu_mcc_fold_mfe[1]
mean_mu_mcc_drtr_mfe = mean_mu_mcc_drtr_mfe[1]
mean_mu_mcc_drtr_prob = mean_mu_mcc_drtr_prob[1]
mean_mu_mcc_fold_bpp = mean_mu_mcc_fold_bpp[1]
mean_mu_mcc_drtr_bpp = mean_mu_mcc_drtr_bpp[1]


mean_n_mcc_fold_mfe = data_n_mean_n[16]
mean_n_mcc_drtr_mfe = data_n_mean_n[17]
mean_n_mcc_drtr_prob = data_n_mean_n[18]
mean_n_mcc_fold_bpp = data_n_mean_n[19]
mean_n_mcc_drtr_bpp = data_n_mean_n[20]

mean_n_mcc_fold_mfe = mean_n_mcc_fold_mfe[1]
mean_n_mcc_drtr_mfe = mean_n_mcc_drtr_mfe[1]
mean_n_mcc_drtr_prob = mean_n_mcc_drtr_prob[1]
mean_n_mcc_fold_bpp = mean_n_mcc_fold_bpp[1]
mean_n_mcc_drtr_bpp = mean_n_mcc_drtr_bpp[1]



boot_m_mcc_fold_mfe = data_n_boot_m[16]
boot_m_mcc_drtr_mfe = data_n_boot_m[17]
boot_m_mcc_drtr_prob = data_n_boot_m[18]
boot_m_mcc_fold_bpp = data_n_boot_m[19]
boot_m_mcc_drtr_bpp = data_n_boot_m[20]

boot_m_mcc_fold_mfe = boot_m_mcc_fold_mfe[1]
boot_m_mcc_drtr_mfe = boot_m_mcc_drtr_mfe[1]
boot_m_mcc_drtr_prob = boot_m_mcc_drtr_prob[1]
boot_m_mcc_fold_bpp = boot_m_mcc_fold_bpp[1]
boot_m_mcc_drtr_bpp = boot_m_mcc_drtr_bpp[1]


boot_mu_mcc_fold_mfe = data_n_boot_mu[16]
boot_mu_mcc_drtr_mfe = data_n_boot_mu[17]
boot_mu_mcc_drtr_prob = data_n_boot_mu[18]
boot_mu_mcc_fold_bpp = data_n_boot_mu[19]
boot_mu_mcc_drtr_bpp = data_n_boot_mu[20]

boot_mu_mcc_fold_mfe = boot_mu_mcc_fold_mfe[1]
boot_mu_mcc_drtr_mfe = boot_mu_mcc_drtr_mfe[1]
boot_mu_mcc_drtr_prob = boot_mu_mcc_drtr_prob[1]
boot_mu_mcc_fold_bpp = boot_mu_mcc_fold_bpp[1]
boot_mu_mcc_drtr_bpp = boot_mu_mcc_drtr_bpp[1]


boot_n_mcc_fold_mfe = data_n_boot_n[16]
boot_n_mcc_drtr_mfe = data_n_boot_n[17]
boot_n_mcc_drtr_prob = data_n_boot_n[18]
boot_n_mcc_fold_bpp = data_n_boot_n[19]
boot_n_mcc_drtr_bpp = data_n_boot_n[20]

boot_n_mcc_fold_mfe = boot_n_mcc_fold_mfe[1]
boot_n_mcc_drtr_mfe = boot_n_mcc_drtr_mfe[1]
boot_n_mcc_drtr_prob = boot_n_mcc_drtr_prob[1]
boot_n_mcc_fold_bpp = boot_n_mcc_fold_bpp[1]
boot_n_mcc_drtr_bpp = boot_n_mcc_drtr_bpp[1]



### ======================== DIFFS ============================ ###
diff_mean_nm_fold_mfe = mean_m_mcc_fold_mfe - mean_n_mcc_fold_mfe
diff_mean_nm_fold_bpp = mean_m_mcc_fold_bpp - mean_n_mcc_fold_bpp

diff_mean_nmu_fold_mfe = mean_mu_mcc_fold_mfe - mean_n_mcc_fold_mfe
diff_mean_nmu_fold_bpp = mean_mu_mcc_fold_bpp - mean_n_mcc_fold_bpp
diff_mean_nmu_drtr_mfe = mean_mu_mcc_drtr_mfe - mean_n_mcc_drtr_mfe
diff_mean_nmu_drtr_prob = mean_mu_mcc_drtr_prob - mean_n_mcc_drtr_prob
diff_mean_nmu_drtr_bpp = mean_mu_mcc_drtr_bpp - mean_n_mcc_drtr_bpp

diff_mean_mmu_fold_mfe = mean_mu_mcc_fold_mfe - mean_m_mcc_fold_mfe
diff_mean_mmu_fold_bpp = mean_mu_mcc_fold_bpp - mean_m_mcc_fold_bpp
diff_mean_mmu_drtr_mfe = mean_mu_mcc_drtr_mfe - mean_m_mcc_drtr_mfe
diff_mean_mmu_drtr_prob = mean_mu_mcc_drtr_prob - mean_m_mcc_drtr_prob
diff_mean_mmu_drtr_bpp = mean_mu_mcc_drtr_bpp - mean_m_mcc_drtr_bpp



diff_means_nm = [
    diff_mean_nm_fold_mfe,
    diff_mean_nm_fold_bpp
    ]


diff_means_nmu = [
                diff_mean_nmu_fold_mfe,
                diff_mean_nmu_fold_bpp,
                ]

diff_means_mmu = [
                # diff_mean_mmu_fold_mfe,
                # diff_mean_mmu_fold_bpp,
                diff_mean_mmu_drtr_mfe,
                diff_mean_mmu_drtr_prob,
                diff_mean_mmu_drtr_bpp,

]

### ======================= plots ========================= ###


indexer = [
        # "fold_mfe",
        # "fold_bpp",
        "drtr_mfe",
        "drtr_prob",
        "drtr_bpp"
        ]

indexer_unblock = [
        "fold_mfe",
        "fold_bpp"
        ]


zero_line_x = [
-1,
1,
2,
3,
5
]

zero_line_y = [
0,
0,
0,
0,
0
]

line_y_1 = [
-0.005,
-0.005,
-0.005,
-0.005,
-0.005
]

line_y_2 = [
0.005,
0.005,
0.005,
0.005,
0.005,
]




fig = plt.figure()
ax = fig.add_subplot()

plt.scatter(indexer_unblock, diff_means_nmu, color = 'blue', label='mu - n')
plt.scatter(indexer, diff_means_mmu, color = 'red', label='mu - m')
plt.scatter(indexer_unblock, diff_means_nm, color = 'green', label='m - n')

# print(diff_means_nmu)
# print(diff_means_nm)

plt.plot(zero_line_x, zero_line_y, linewidth = 0.1, color="black")
plt.plot(zero_line_x, line_y_1, linewidth = 0.1, color="black")
plt.plot(zero_line_x, line_y_2, linewidth = 0.1, color="black")

# ax.fill_between(indexer, diff_boots_lq_c_mcc, diff_boots_hq_c_mcc, color = 'green', alpha=0.1)
plt.axis([-1, 5, -0.01, 0.01])
ax.set_xlabel('Algorithm', fontsize = 12)
ax.set_ylabel('MCC difference', fontsize = 12)
# ax.set_title('Effect of removing non canonical bps and pseudoknots \n on performance based MCC', fontsize=14)
plt.gcf().subplots_adjust(bottom=0.12, left=0.15, right=0.94, top=0.90)
plt.legend()
plt.title("Effect of modified bases on prediction performance \n MCC difference")
# print(diff_means_mn)
# print(diff_means_mmn)
# print(diff_means_mun)
# print(diff_means_nmu)

if args.graphics==True:
    plt.show()
