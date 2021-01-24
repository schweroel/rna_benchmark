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

### outputsettings

pd.set_option('display.max_rows', 10000) #some default pandas settings will chop data
pd.set_option('display.max_columns', 10000)
pd.set_option('display.width', 10000)
np.set_printoptions(threshold=sys.maxsize) #if printed out, matrix won't be truncated


### ====================================== ARGPARSE ==================================== ###
parser = argparse.ArgumentParser(description='Compare ct data with RNAfold')
# parser.add_argument('--verbose','-v', action='store_true',
#                     default=False,
#                     help='be verbose')

parser.add_argument('--graphics','-g', action='store_true',
                    default=False,
                    help='make grap')

parser.add_argument('--printout','-p', action='store_true',
                    default=False,
                    help='print values')

parser.add_argument('--test', '-t', action='store_true',
                    default=False,
                    help='run with test csv')

parser.add_argument('--label','-l', action='store',
                    default="mod_unblock_comp_m_n",
                    help='parse label to graph')

parser.add_argument('--input','-i', action='store',
                    default="n",
                    help='mu, m or n')

parser.add_argument('--output','-o', action='store_true',
                    default=False,
                    help='recalculate boot/csvs')

args=parser.parse_args()

### ==================================== NAMES & IMPORT ================================ ###
input = args.input

input_name = "comp_" + args.input + "_n_fin_s_boot.csv"
basename = input_name [:-10]

label = args.input

i_n_mean = "comp_" + args.input + "_n_fin_s_mean.csv"
i_c_mean = "comp_" + args.input + "_c_fin_s_mean.csv"
i_p_mean = "comp_" + args.input + "_p_fin_s_mean.csv"
i_cp_mean = "comp_" + args.input + "_cp_fin_s_mean.csv"

i_n_boot = "comp_" + args.input + "_n_fin_s_boot.csv"
i_c_boot = "comp_" + args.input + "_c_fin_s_boot.csv"
i_p_boot = "comp_" + args.input + "_p_fin_s_boot.csv"
i_cp_boot = "comp_" + args.input + "_cp_fin_s_boot.csv"

data_i_n_mean = genfromtxt(i_n_mean, delimiter = ",")
data_i_c_mean = genfromtxt(i_c_mean, delimiter = ",")
data_i_p_mean = genfromtxt(i_p_mean, delimiter = ",")
data_i_cp_mean = genfromtxt(i_cp_mean, delimiter = ",")

data_i_n_boot = genfromtxt(i_n_boot, delimiter = ",")
data_i_c_boot = genfromtxt(i_c_boot, delimiter = ",")
data_i_p_boot = genfromtxt(i_p_boot, delimiter = ",")
data_i_cp_boot = genfromtxt(i_cp_boot, delimiter = ",")

# print(data_i_n_mean)

data_200 = data_i_n_mean[1,1:]
data_500 = data_i_n_mean[2,1:]
data_1000 = data_i_n_mean[3,1:]
data_1000_plus = data_i_n_mean[4,1:]

data_200_cp = data_i_cp_mean[1,1:]
data_500_cp = data_i_cp_mean[2,1:]
data_1000_cp = data_i_cp_mean[3,1:]
data_1000_plus_cp = data_i_cp_mean[4,1:]

### ================================================================= SET DATA VARIABLES ============================================================ ###

ppv_200_fold_mfe = data_200[0]
ppv_200_fold_bpp = data_200[1]
ppv_200_drtr_mfe = data_200[2]
ppv_200_drtr_prob = data_200[3]
ppv_200_drtr_bpp= data_200[4]
tpr_200_fold_mfe = data_200[5]
tpr_200_fold_bpp = data_200[6]
tpr_200_drtr_mfe = data_200[7]
tpr_200_drtr_prob = data_200[8]
tpr_200_drtr_bpp= data_200[9]
f1_200_fold_mfe = data_200[10]
f1_200_fold_bpp = data_200[11]
f1_200_drtr_mfe = data_200[12]
f1_200_drtr_prob = data_200[13]
f1_200_drtr_bpp= data_200[14]
mcc_200_fold_mfe = data_200[15]
mcc_200_fold_bpp = data_200[16]
mcc_200_drtr_mfe = data_200[17]
mcc_200_drtr_prob = data_200[18]
mcc_200_drtr_bpp= data_200[19]

mcc_200_fold_mfe_cp = data_200_cp[15]
mcc_200_fold_bpp_cp = data_200_cp[16]
mcc_200_drtr_mfe_cp = data_200_cp[17]
mcc_200_drtr_prob_cp = data_200_cp[18]
mcc_200_drtr_bpp_cp = data_200_cp[19]


ppv_500_fold_mfe = data_500[0]
ppv_500_fold_bpp = data_500[1]
ppv_500_drtr_mfe = data_500[2]
ppv_500_drtr_prob = data_500[3]
ppv_500_drtr_bpp= data_500[4]
tpr_500_fold_mfe = data_500[5]
tpr_500_fold_bpp = data_500[6]
tpr_500_drtr_mfe = data_500[7]
tpr_500_drtr_prob = data_500[8]
tpr_500_drtr_bpp= data_500[9]
f1_500_fold_mfe = data_500[10]
f1_500_fold_bpp = data_500[11]
f1_500_drtr_mfe = data_500[12]
f1_500_drtr_prob = data_500[13]
f1_500_drtr_bpp= data_500[14]
mcc_500_fold_mfe = data_500[15]
mcc_500_fold_bpp = data_500[16]
mcc_500_drtr_mfe = data_500[17]
mcc_500_drtr_prob = data_500[18]
mcc_500_drtr_bpp= data_500[19]

mcc_500_fold_mfe_cp = data_500_cp[15]
mcc_500_fold_bpp_cp = data_500_cp[16]
mcc_500_drtr_mfe_cp = data_500_cp[17]
mcc_500_drtr_prob_cp = data_500_cp[18]
mcc_500_drtr_bpp_cp = data_500_cp[19]


ppv_1000_fold_mfe = data_1000[0]
ppv_1000_fold_bpp = data_1000[1]
ppv_1000_drtr_mfe = data_1000[2]
ppv_1000_drtr_prob = data_1000[3]
ppv_1000_drtr_bpp= data_1000[4]
tpr_1000_fold_mfe = data_1000[5]
tpr_1000_fold_bpp = data_1000[6]
tpr_1000_drtr_mfe = data_1000[7]
tpr_1000_drtr_prob = data_1000[8]
tpr_1000_drtr_bpp= data_1000[9]
f1_1000_fold_mfe = data_1000[10]
f1_1000_fold_bpp = data_1000[11]
f1_1000_drtr_mfe = data_1000[12]
f1_1000_drtr_prob = data_1000[13]
f1_1000_drtr_bpp= data_1000[14]
mcc_1000_fold_mfe = data_1000[15]
mcc_1000_fold_bpp = data_1000[16]
mcc_1000_drtr_mfe = data_1000[17]
mcc_1000_drtr_prob = data_1000[18]
mcc_1000_drtr_bpp= data_1000[19]

mcc_1000_fold_mfe_cp = data_1000_cp[15]
mcc_1000_fold_bpp_cp = data_1000_cp[16]
mcc_1000_drtr_mfe_cp = data_1000_cp[17]
mcc_1000_drtr_prob_cp = data_1000_cp[18]
mcc_1000_drtr_bpp_cp = data_1000_cp[19]


ppv_1000_plus_fold_mfe = data_1000_plus[0]
ppv_1000_plus_fold_bpp = data_1000_plus[1]
ppv_1000_plus_drtr_mfe = data_1000_plus[2]
ppv_1000_plus_drtr_prob = data_1000_plus[3]
ppv_1000_plus_drtr_bpp= data_1000_plus[4]
tpr_1000_plus_fold_mfe = data_1000_plus[5]
tpr_1000_plus_fold_bpp = data_1000_plus[6]
tpr_1000_plus_drtr_mfe = data_1000_plus[7]
tpr_1000_plus_drtr_prob = data_1000_plus[8]
tpr_1000_plus_drtr_bpp= data_1000_plus[9]
f1_1000_plus_fold_mfe = data_1000_plus[10]
f1_1000_plus_fold_bpp = data_1000_plus[11]
f1_1000_plus_drtr_mfe = data_1000_plus[12]
f1_1000_plus_drtr_prob = data_1000_plus[13]
f1_1000_plus_drtr_bpp= data_1000_plus[14]
mcc_1000_plus_fold_mfe = data_1000_plus[15]
mcc_1000_plus_fold_bpp = data_1000_plus[16]
mcc_1000_plus_drtr_mfe = data_1000_plus[17]
mcc_1000_plus_drtr_prob = data_1000_plus[18]
mcc_1000_plus_drtr_bpp= data_1000_plus[19]

mcc_1000_plus_fold_mfe_cp = data_1000_plus_cp[15]
mcc_1000_plus_fold_bpp_cp = data_1000_plus_cp[16]
mcc_1000_plus_drtr_mfe_cp = data_1000_plus_cp[17]
mcc_1000_plus_drtr_prob_cp = data_1000_plus_cp[18]
mcc_1000_plus_drtr_bpp_cp = data_1000_plus_cp[19]

###========================================================== MAKING LISTS OF VARIABLES FOR PLOTTING =========================================== ###
fold_mfe_mcc_bars = [
mcc_200_fold_mfe,
mcc_500_fold_mfe,
mcc_1000_fold_mfe,
mcc_1000_plus_fold_mfe
]

fold_bpp_mcc_bars = [
mcc_200_fold_bpp,
mcc_500_fold_bpp,
mcc_1000_fold_bpp,
mcc_1000_plus_fold_bpp
]

drtr_mfe_mcc_bars = [
mcc_200_drtr_mfe,
mcc_500_drtr_mfe,
mcc_1000_drtr_mfe,
mcc_1000_plus_drtr_mfe
]

drtr_prob_mcc_bars = [
mcc_200_drtr_prob,
mcc_500_drtr_prob,
mcc_1000_drtr_prob,
mcc_1000_plus_drtr_prob
]

drtr_bpp_mcc_bars = [
mcc_200_drtr_bpp,
mcc_500_drtr_bpp,
mcc_1000_drtr_bpp,
mcc_1000_plus_drtr_bpp
]



fold_mfe_mcc_bars_cp = [
mcc_200_fold_mfe_cp,
mcc_500_fold_mfe_cp,
mcc_1000_fold_mfe_cp,
mcc_1000_plus_fold_mfe_cp
]

fold_bpp_mcc_bars_cp = [
mcc_200_fold_bpp_cp,
mcc_500_fold_bpp_cp,
mcc_1000_fold_bpp_cp,
mcc_1000_plus_fold_bpp_cp
]

drtr_mfe_mcc_bars_cp = [
mcc_200_drtr_mfe_cp,
mcc_500_drtr_mfe_cp,
mcc_1000_drtr_mfe_cp,
mcc_1000_plus_drtr_mfe_cp
]

drtr_prob_mcc_bars_cp = [
mcc_200_drtr_prob_cp,
mcc_500_drtr_prob_cp,
mcc_1000_drtr_prob_cp,
mcc_1000_plus_drtr_prob_cp
]

drtr_bpp_mcc_bars_cp = [
mcc_200_drtr_bpp_cp,
mcc_500_drtr_bpp_cp,
mcc_1000_drtr_bpp_cp,
mcc_1000_plus_drtr_bpp_cp
]

if args.graphics==True:


    barwidth = 0.18

    r_1 = np.arange(len(fold_mfe_mcc_bars))
    r_2 = [x + barwidth for x in r_1]
    r_3 = [x + barwidth for x in r_2]
    r_4 = [x + barwidth for x in r_3]
    r_5 = [x + barwidth for x in r_4]

    plt.bar(r_1, fold_mfe_mcc_bars_cp, width=barwidth, color='#a8a8a8', label='canonical bps and pks removed')
    plt.bar(r_2, drtr_mfe_mcc_bars_cp, width=barwidth, color='#a8a8a8')
    plt.bar(r_3, drtr_prob_mcc_bars_cp, width=barwidth, color='#a8a8a8')
    plt.bar(r_4, fold_bpp_mcc_bars_cp, width=barwidth, color='#a8a8a8')
    plt.bar(r_5, drtr_bpp_mcc_bars_cp, width=barwidth, color='#a8a8a8')

    plt.bar(r_1, fold_mfe_mcc_bars, width=barwidth, color='#1d93ad', label='RNAfold: MFE')
    plt.bar(r_2, drtr_mfe_mcc_bars, width=barwidth, color='#BE3151', label='DrTransformer: MFE')
    plt.bar(r_3, drtr_prob_mcc_bars, width=barwidth, color='#1dad78', label='DrTransformer: max prop')
    plt.bar(r_4, fold_bpp_mcc_bars, width=barwidth, color='#e3ed21', label='RNAfold: BPP')
    plt.bar(r_5, drtr_bpp_mcc_bars, width=barwidth, color='#0E3151', label='DrTransformer: BPP')

    plt.gcf().subplots_adjust(bottom=0.32, left=0.11, right=0.91, top=0.91)

    plt.xticks([r + barwidth for r in range (len(fold_mfe_mcc_bars))],
                                                   ["200 (n=1428)",
                                                   "500 (n=1411)",
                                                   "1000 (n=104)",
                                                   "1000+ (n=6)"],

                                                   rotation=55,
                                                   ha="right",
                                                   rotation_mode="anchor",
                                                   fontsize=12)

    plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    plt.title("MCC grouped by sequence-length (" + label + ")" )
    plt.ylabel("MCC", fontsize=12)
    plt.xlabel("Number of nucleotides", fontsize=12)
    plt.legend()
    plt.show()

if args.printout==True:
    print("fold_mfe")
    print(fold_mfe_mcc_bars)
    print("")
    print("fold_bpp")
    print(fold_bpp_mcc_bars)
    print("")
    print("drtr_mfe")
    print(drtr_mfe_mcc_bars)
    print("")
    print("drtr_bpp")
    print(drtr_bpp_mcc_bars)
    print("")
    print("drtr_prob")
    print(drtr_prob_mcc_bars)
    print("")
