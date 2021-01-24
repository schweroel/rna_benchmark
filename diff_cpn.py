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

parser.add_argument('--input','-i', action='store',
                    default="n",
                    help='mu, m or n')

parser.add_argument('--output','-o', action='store_true',
                    default=False,
                    help='recalculate boot/csvs')

args=parser.parse_args()

### ==================================== NAMES & IMPORT ================================ ###
input = args.input

input_name = "comp_" + args.input + "_n_mfe_prob_n_boot.csv"
basename = input_name [:-10]

label = args.input

i_n_mean = "comp_" + args.input + "_n_mfe_prob_n_mean.csv"
i_c_mean = "comp_" + args.input + "_c_mfe_prob_n_mean.csv"
i_p_mean = "comp_" + args.input + "_p_mfe_prob_n_mean.csv"
i_cp_mean = "comp_" + args.input + "_cp_mfe_prob_n_mean.csv"

i_n_boot = "comp_" + args.input + "_n_mfe_prob_n_boot.csv"
i_c_boot = "comp_" + args.input + "_c_mfe_prob_n_boot.csv"
i_p_boot = "comp_" + args.input + "_p_mfe_prob_n_boot.csv"
i_cp_boot = "comp_" + args.input + "_cp_mfe_prob_n_boot.csv"

data_i_n_mean = genfromtxt(i_n_mean, delimiter = ",")
data_i_c_mean = genfromtxt(i_c_mean, delimiter = ",")
data_i_p_mean = genfromtxt(i_p_mean, delimiter = ",")
data_i_cp_mean = genfromtxt(i_cp_mean, delimiter = ",")

data_i_n_boot = genfromtxt(i_n_boot, delimiter = ",")
data_i_c_boot = genfromtxt(i_c_boot, delimiter = ",")
data_i_p_boot = genfromtxt(i_p_boot, delimiter = ",")
data_i_cp_boot = genfromtxt(i_cp_boot, delimiter = ",")

mean_ppv_n_fold_mfe = data_i_n_mean[1]
mean_ppv_n_drtr_mfe = data_i_n_mean[2]
mean_ppv_n_drtr_prob = data_i_n_mean[3]
mean_ppv_n_fold_bpp = data_i_n_mean[4]
mean_ppv_n_drtr_bpp = data_i_n_mean[5]

mean_tpr_n_fold_mfe = data_i_n_mean[6]
mean_tpr_n_drtr_mfe = data_i_n_mean[7]
mean_tpr_n_drtr_prob = data_i_n_mean[8]
mean_tpr_n_fold_bpp = data_i_n_mean[9]
mean_tpr_n_drtr_bpp = data_i_n_mean[10]

mean_f_1_n_fold_mfe = data_i_n_mean[11]
mean_f_1_n_drtr_mfe = data_i_n_mean[12]
mean_f_1_n_drtr_prob = data_i_n_mean[13]
mean_f_1_n_fold_bpp = data_i_n_mean[14]
mean_f_1_n_drtr_bpp = data_i_n_mean[15]

mean_mcc_n_fold_mfe = data_i_n_mean[16]
mean_mcc_n_drtr_mfe = data_i_n_mean[17]
mean_mcc_n_drtr_prob = data_i_n_mean[18]
mean_mcc_n_fold_bpp = data_i_n_mean[19]
mean_mcc_n_drtr_bpp = data_i_n_mean[20]

mean_ppv_c_fold_mfe = data_i_c_mean[1]
mean_ppv_c_drtr_mfe = data_i_c_mean[2]
mean_ppv_c_drtr_prob = data_i_c_mean[3]
mean_ppv_c_fold_bpp = data_i_c_mean[4]
mean_ppv_c_drtr_bpp = data_i_c_mean[5]

mean_tpr_c_fold_mfe = data_i_c_mean[6]
mean_tpr_c_drtr_mfe = data_i_c_mean[7]
mean_tpr_c_drtr_prob = data_i_c_mean[8]
mean_tpr_c_fold_bpp = data_i_c_mean[9]
mean_tpr_c_drtr_bpp = data_i_c_mean[10]

mean_f_1_c_fold_mfe = data_i_c_mean[11]
mean_f_1_c_drtr_mfe = data_i_c_mean[12]
mean_f_1_c_drtr_prob = data_i_c_mean[13]
mean_f_1_c_fold_bpp = data_i_c_mean[14]
mean_f_1_c_drtr_bpp = data_i_c_mean[15]

mean_mcc_c_fold_mfe = data_i_c_mean[16]
mean_mcc_c_drtr_mfe = data_i_c_mean[17]
mean_mcc_c_drtr_prob = data_i_c_mean[18]
mean_mcc_c_fold_bpp = data_i_c_mean[19]
mean_mcc_c_drtr_bpp = data_i_c_mean[20]

mean_ppv_p_fold_mfe = data_i_p_mean[1]
mean_ppv_p_drtr_mfe = data_i_p_mean[2]
mean_ppv_p_drtr_prob = data_i_p_mean[3]
mean_ppv_p_fold_bpp = data_i_p_mean[4]
mean_ppv_p_drtr_bpp = data_i_p_mean[5]

mean_tpr_p_fold_mfe = data_i_p_mean[6]
mean_tpr_p_drtr_mfe = data_i_p_mean[7]
mean_tpr_p_drtr_prob = data_i_p_mean[8]
mean_tpr_p_fold_bpp = data_i_p_mean[9]
mean_tpr_p_drtr_bpp = data_i_p_mean[10]

mean_f_1_p_fold_mfe = data_i_p_mean[11]
mean_f_1_p_drtr_mfe = data_i_p_mean[12]
mean_f_1_p_drtr_prob = data_i_p_mean[13]
mean_f_1_p_fold_bpp = data_i_p_mean[14]
mean_f_1_p_drtr_bpp = data_i_p_mean[15]

mean_mcc_p_fold_mfe = data_i_p_mean[16]
mean_mcc_p_drtr_mfe = data_i_p_mean[17]
mean_mcc_p_drtr_prob = data_i_p_mean[18]
mean_mcc_p_fold_bpp = data_i_p_mean[19]
mean_mcc_p_drtr_bpp = data_i_p_mean[20]

mean_ppv_cp_fold_mfe = data_i_cp_mean[1]
mean_ppv_cp_drtr_mfe = data_i_cp_mean[2]
mean_ppv_cp_drtr_prob = data_i_cp_mean[3]
mean_ppv_cp_fold_bpp = data_i_cp_mean[4]
mean_ppv_cp_drtr_bpp = data_i_cp_mean[5]

mean_tpr_cp_fold_mfe = data_i_cp_mean[6]
mean_tpr_cp_drtr_mfe = data_i_cp_mean[7]
mean_tpr_cp_drtr_prob = data_i_cp_mean[8]
mean_tpr_cp_fold_bpp = data_i_cp_mean[9]
mean_tpr_cp_drtr_bpp = data_i_cp_mean[10]

mean_f_1_cp_fold_mfe = data_i_cp_mean[11]
mean_f_1_cp_drtr_mfe = data_i_cp_mean[12]
mean_f_1_cp_drtr_prob = data_i_cp_mean[13]
mean_f_1_cp_fold_bpp = data_i_cp_mean[14]
mean_f_1_cp_drtr_bpp = data_i_cp_mean[15]

mean_mcc_cp_fold_mfe = data_i_cp_mean[16]
mean_mcc_cp_drtr_mfe = data_i_cp_mean[17]
mean_mcc_cp_drtr_prob = data_i_cp_mean[18]
mean_mcc_cp_fold_bpp = data_i_cp_mean[19]
mean_mcc_cp_drtr_bpp = data_i_cp_mean[20]

### boot n

boot_n_ppv_fold_mfe = data_i_n_boot[1]
boot_n_ppv_fold_mfe = np.sort(boot_n_ppv_fold_mfe)
lq_n_ppv_fold_mfe = boot_n_ppv_fold_mfe[249]
hq_n_ppv_fold_mfe = boot_n_ppv_fold_mfe[749]

boot_n_ppv_drtr_mfe = data_i_n_boot[2]
boot_n_ppv_drtr_mfe = np.sort(boot_n_ppv_drtr_mfe)
lq_n_ppv_drtr_mfe = boot_n_ppv_drtr_mfe[249]
hq_n_ppv_drtr_mfe = boot_n_ppv_drtr_mfe[749]

boot_n_ppv_drtr_prob = data_i_n_boot[3]
boot_n_ppv_drtr_prob = np.sort(boot_n_ppv_drtr_prob)
lq_n_ppv_drtr_prob = boot_n_ppv_drtr_prob[249]
hq_n_ppv_drtr_prob = boot_n_ppv_drtr_prob[749]

boot_n_ppv_fold_bpp = data_i_n_boot[4]
boot_n_ppv_fold_bpp = np.sort(boot_n_ppv_fold_bpp)
lq_n_ppv_fold_bpp = boot_n_ppv_fold_bpp[249]
hq_n_ppv_fold_bpp = boot_n_ppv_fold_bpp[749]

boot_n_ppv_drtr_bpp = data_i_n_boot[5]
boot_n_ppv_drtr_bpp = np.sort(boot_n_ppv_drtr_bpp)
lq_n_ppv_drtr_bpp = boot_n_ppv_drtr_bpp[249]
hq_n_ppv_drtr_bpp = boot_n_ppv_drtr_bpp[749]


boot_n_tpr_fold_mfe = data_i_n_boot[6]
boot_n_tpr_fold_mfe = np.sort(boot_n_tpr_fold_mfe)
lq_n_tpr_fold_mfe = boot_n_tpr_fold_mfe[249]
hq_n_tpr_fold_mfe = boot_n_tpr_fold_mfe[749]

boot_n_tpr_drtr_mfe = data_i_n_boot[7]
boot_n_tpr_drtr_mfe = np.sort(boot_n_tpr_drtr_mfe)
lq_n_tpr_drtr_mfe = boot_n_tpr_drtr_mfe[249]
hq_n_tpr_drtr_mfe = boot_n_tpr_drtr_mfe[749]

boot_n_tpr_drtr_prob = data_i_n_boot[8]
boot_n_tpr_drtr_prob = np.sort(boot_n_tpr_drtr_prob)
lq_n_tpr_drtr_prob = boot_n_tpr_drtr_prob[249]
hq_n_tpr_drtr_prob = boot_n_tpr_drtr_prob[749]

boot_n_tpr_fold_bpp = data_i_n_boot[9]
boot_n_tpr_fold_bpp = np.sort(boot_n_tpr_fold_bpp)
lq_n_tpr_fold_bpp = boot_n_tpr_fold_bpp[249]
hq_n_tpr_fold_bpp = boot_n_tpr_fold_bpp[749]

boot_n_tpr_drtr_bpp = data_i_n_boot[10]
boot_n_tpr_drtr_bpp = np.sort(boot_n_tpr_drtr_bpp)
lq_n_tpr_drtr_bpp = boot_n_tpr_drtr_bpp[249]
hq_n_tpr_drtr_bpp = boot_n_tpr_drtr_bpp[749]


boot_n_f_1_fold_mfe = data_i_n_boot[11]
boot_n_f_1_fold_mfe = np.sort(boot_n_f_1_fold_mfe)
lq_n_f_1_fold_mfe = boot_n_f_1_fold_mfe[249]
hq_n_f_1_fold_mfe = boot_n_f_1_fold_mfe[749]

boot_n_f_1_drtr_mfe = data_i_n_boot[12]
boot_n_f_1_drtr_mfe = np.sort(boot_n_f_1_drtr_mfe)
lq_n_f_1_drtr_mfe = boot_n_f_1_drtr_mfe[249]
hq_n_f_1_drtr_mfe = boot_n_f_1_drtr_mfe[749]

boot_n_f_1_drtr_prob = data_i_n_boot[13]
boot_n_f_1_drtr_prob = np.sort(boot_n_f_1_drtr_prob)
lq_n_f_1_drtr_prob = boot_n_f_1_drtr_prob[249]
hq_n_f_1_drtr_prob = boot_n_f_1_drtr_prob[749]

boot_n_f_1_fold_bpp = data_i_n_boot[14]
boot_n_f_1_fold_bpp = np.sort(boot_n_f_1_fold_bpp)
lq_n_f_1_fold_bpp = boot_n_f_1_fold_bpp[249]
hq_n_f_1_fold_bpp = boot_n_f_1_fold_bpp[749]

boot_n_f_1_drtr_bpp = data_i_n_boot[15]
boot_n_f_1_drtr_bpp = np.sort(boot_n_f_1_drtr_bpp)
lq_n_f_1_drtr_bpp = boot_n_f_1_drtr_bpp[249]
hq_n_f_1_drtr_bpp = boot_n_f_1_drtr_bpp[749]


boot_n_mcc_fold_mfe = data_i_n_boot[16]
boot_n_mcc_fold_mfe = np.sort(boot_n_mcc_fold_mfe)
lq_n_mcc_fold_mfe = boot_n_mcc_fold_mfe[249]
hq_n_mcc_fold_mfe = boot_n_mcc_fold_mfe[749]

boot_n_mcc_drtr_mfe = data_i_n_boot[17]
boot_n_mcc_drtr_mfe = np.sort(boot_n_mcc_drtr_mfe)
lq_n_mcc_drtr_mfe = boot_n_mcc_drtr_mfe[249]
hq_n_mcc_drtr_mfe = boot_n_mcc_drtr_mfe[749]

boot_n_mcc_drtr_prob = data_i_n_boot[18]
boot_n_mcc_drtr_prob = np.sort(boot_n_mcc_drtr_prob)
lq_n_mcc_drtr_prob = boot_n_mcc_drtr_prob[249]
hq_n_mcc_drtr_prob = boot_n_mcc_drtr_prob[749]

boot_n_mcc_fold_bpp = data_i_n_boot[19]
boot_n_mcc_fold_bpp = np.sort(boot_n_mcc_fold_bpp)
lq_n_mcc_fold_bpp = boot_n_mcc_fold_bpp[249]
hq_n_mcc_fold_bpp = boot_n_mcc_fold_bpp[749]

boot_n_mcc_drtr_bpp = data_i_n_boot[20]
boot_n_mcc_drtr_bpp = np.sort(boot_n_mcc_drtr_bpp)
lq_n_mcc_drtr_bpp = boot_n_mcc_drtr_bpp[249]
hq_n_mcc_drtr_bpp = boot_n_mcc_drtr_bpp[749]

### boot c

boot_c_ppv_fold_mfe = data_i_c_boot[1]
boot_c_ppv_fold_mfe = np.sort(boot_c_ppv_fold_mfe)
lq_c_ppv_fold_mfe = boot_c_ppv_fold_mfe[249]
hq_c_ppv_fold_mfe = boot_c_ppv_fold_mfe[749]

boot_c_ppv_drtr_mfe = data_i_c_boot[2]
boot_c_ppv_drtr_mfe = np.sort(boot_c_ppv_drtr_mfe)
lq_c_ppv_drtr_mfe = boot_c_ppv_drtr_mfe[249]
hq_c_ppv_drtr_mfe = boot_c_ppv_drtr_mfe[749]

boot_c_ppv_drtr_prob = data_i_c_boot[3]
boot_c_ppv_drtr_prob = np.sort(boot_c_ppv_drtr_prob)
lq_c_ppv_drtr_prob = boot_c_ppv_drtr_prob[249]
hq_c_ppv_drtr_prob = boot_c_ppv_drtr_prob[749]

boot_c_ppv_fold_bpp = data_i_c_boot[4]
boot_c_ppv_fold_bpp = np.sort(boot_c_ppv_fold_bpp)
lq_c_ppv_fold_bpp = boot_c_ppv_fold_bpp[249]
hq_c_ppv_fold_bpp = boot_c_ppv_fold_bpp[749]

boot_c_ppv_drtr_bpp = data_i_c_boot[5]
boot_c_ppv_drtr_bpp = np.sort(boot_c_ppv_drtr_bpp)
lq_c_ppv_drtr_bpp = boot_c_ppv_drtr_bpp[249]
hq_c_ppv_drtr_bpp = boot_c_ppv_drtr_bpp[749]


boot_c_tpr_fold_mfe = data_i_c_boot[6]
boot_c_tpr_fold_mfe = np.sort(boot_c_tpr_fold_mfe)
lq_c_tpr_fold_mfe = boot_c_tpr_fold_mfe[249]
hq_c_tpr_fold_mfe = boot_c_tpr_fold_mfe[749]

boot_c_tpr_drtr_mfe = data_i_c_boot[7]
boot_c_tpr_drtr_mfe = np.sort(boot_c_tpr_drtr_mfe)
lq_c_tpr_drtr_mfe = boot_c_tpr_drtr_mfe[249]
hq_c_tpr_drtr_mfe = boot_c_tpr_drtr_mfe[749]

boot_c_tpr_drtr_prob = data_i_c_boot[8]
boot_c_tpr_drtr_prob = np.sort(boot_c_tpr_drtr_prob)
lq_c_tpr_drtr_prob = boot_c_tpr_drtr_prob[249]
hq_c_tpr_drtr_prob = boot_c_tpr_drtr_prob[749]

boot_c_tpr_fold_bpp = data_i_c_boot[9]
boot_c_tpr_fold_bpp = np.sort(boot_c_tpr_fold_bpp)
lq_c_tpr_fold_bpp = boot_c_tpr_fold_bpp[249]
hq_c_tpr_fold_bpp = boot_c_tpr_fold_bpp[749]

boot_c_tpr_drtr_bpp = data_i_c_boot[10]
boot_c_tpr_drtr_bpp = np.sort(boot_c_tpr_drtr_bpp)
lq_c_tpr_drtr_bpp = boot_c_tpr_drtr_bpp[249]
hq_c_tpr_drtr_bpp = boot_c_tpr_drtr_bpp[749]


boot_c_f_1_fold_mfe = data_i_c_boot[11]
boot_c_f_1_fold_mfe = np.sort(boot_c_f_1_fold_mfe)
lq_c_f_1_fold_mfe = boot_c_f_1_fold_mfe[249]
hq_c_f_1_fold_mfe = boot_c_f_1_fold_mfe[749]

boot_c_f_1_drtr_mfe = data_i_c_boot[12]
boot_c_f_1_drtr_mfe = np.sort(boot_c_f_1_drtr_mfe)
lq_c_f_1_drtr_mfe = boot_c_f_1_drtr_mfe[249]
hq_c_f_1_drtr_mfe = boot_c_f_1_drtr_mfe[749]

boot_c_f_1_drtr_prob = data_i_c_boot[13]
boot_c_f_1_drtr_prob = np.sort(boot_c_f_1_drtr_prob)
lq_c_f_1_drtr_prob = boot_c_f_1_drtr_prob[249]
hq_c_f_1_drtr_prob = boot_c_f_1_drtr_prob[749]

boot_c_f_1_fold_bpp = data_i_c_boot[14]
boot_c_f_1_fold_bpp = np.sort(boot_c_f_1_fold_bpp)
lq_c_f_1_fold_bpp = boot_c_f_1_fold_bpp[249]
hq_c_f_1_fold_bpp = boot_c_f_1_fold_bpp[749]

boot_c_f_1_drtr_bpp = data_i_c_boot[15]
boot_c_f_1_drtr_bpp = np.sort(boot_c_f_1_drtr_bpp)
lq_c_f_1_drtr_bpp = boot_c_f_1_drtr_bpp[249]
hq_c_f_1_drtr_bpp = boot_c_f_1_drtr_bpp[749]


boot_c_mcc_fold_mfe = data_i_c_boot[16]
boot_c_mcc_fold_mfe = np.sort(boot_c_mcc_fold_mfe)
lq_c_mcc_fold_mfe = boot_c_mcc_fold_mfe[249]
hq_c_mcc_fold_mfe = boot_c_mcc_fold_mfe[749]

boot_c_mcc_drtr_mfe = data_i_c_boot[17]
boot_c_mcc_drtr_mfe = np.sort(boot_c_mcc_drtr_mfe)
lq_c_mcc_drtr_mfe = boot_c_mcc_drtr_mfe[249]
hq_c_mcc_drtr_mfe = boot_c_mcc_drtr_mfe[749]

boot_c_mcc_drtr_prob = data_i_c_boot[18]
boot_c_mcc_drtr_prob = np.sort(boot_c_mcc_drtr_prob)
lq_c_mcc_drtr_prob = boot_c_mcc_drtr_prob[249]
hq_c_mcc_drtr_prob = boot_c_mcc_drtr_prob[749]

boot_c_mcc_fold_bpp = data_i_c_boot[19]
boot_c_mcc_fold_bpp = np.sort(boot_c_mcc_fold_bpp)
lq_c_mcc_fold_bpp = boot_c_mcc_fold_bpp[249]
hq_c_mcc_fold_bpp = boot_c_mcc_fold_bpp[749]

boot_c_mcc_drtr_bpp = data_i_c_boot[20]
boot_c_mcc_drtr_bpp = np.sort(boot_c_mcc_drtr_bpp)
lq_c_mcc_drtr_bpp = boot_c_mcc_drtr_bpp[249]
hq_c_mcc_drtr_bpp = boot_c_mcc_drtr_bpp[749]


### boot p

boot_p_ppv_fold_mfe = data_i_p_boot[1]
boot_p_ppv_fold_mfe = np.sort(boot_p_ppv_fold_mfe)
lq_p_ppv_fold_mfe = boot_p_ppv_fold_mfe[249]
hq_p_ppv_fold_mfe = boot_p_ppv_fold_mfe[749]

boot_p_ppv_drtr_mfe = data_i_p_boot[2]
boot_p_ppv_drtr_mfe = np.sort(boot_p_ppv_drtr_mfe)
lq_p_ppv_drtr_mfe = boot_p_ppv_drtr_mfe[249]
hq_p_ppv_drtr_mfe = boot_p_ppv_drtr_mfe[749]

boot_p_ppv_drtr_prob = data_i_p_boot[3]
boot_p_ppv_drtr_prob = np.sort(boot_p_ppv_drtr_prob)
lq_p_ppv_drtr_prob = boot_p_ppv_drtr_prob[249]
hq_p_ppv_drtr_prob = boot_p_ppv_drtr_prob[749]

boot_p_ppv_fold_bpp = data_i_p_boot[4]
boot_p_ppv_fold_bpp = np.sort(boot_p_ppv_fold_bpp)
lq_p_ppv_fold_bpp = boot_p_ppv_fold_bpp[249]
hq_p_ppv_fold_bpp = boot_p_ppv_fold_bpp[749]

boot_p_ppv_drtr_bpp = data_i_p_boot[5]
boot_p_ppv_drtr_bpp = np.sort(boot_p_ppv_drtr_bpp)
lq_p_ppv_drtr_bpp = boot_p_ppv_drtr_bpp[249]
hq_p_ppv_drtr_bpp = boot_p_ppv_drtr_bpp[749]


boot_p_tpr_fold_mfe = data_i_p_boot[6]
boot_p_tpr_fold_mfe = np.sort(boot_p_tpr_fold_mfe)
lq_p_tpr_fold_mfe = boot_p_tpr_fold_mfe[249]
hq_p_tpr_fold_mfe = boot_p_tpr_fold_mfe[749]

boot_p_tpr_drtr_mfe = data_i_p_boot[7]
boot_p_tpr_drtr_mfe = np.sort(boot_p_tpr_drtr_mfe)
lq_p_tpr_drtr_mfe = boot_p_tpr_drtr_mfe[249]
hq_p_tpr_drtr_mfe = boot_p_tpr_drtr_mfe[749]

boot_p_tpr_drtr_prob = data_i_p_boot[8]
boot_p_tpr_drtr_prob = np.sort(boot_p_tpr_drtr_prob)
lq_p_tpr_drtr_prob = boot_p_tpr_drtr_prob[249]
hq_p_tpr_drtr_prob = boot_p_tpr_drtr_prob[749]

boot_p_tpr_fold_bpp = data_i_p_boot[9]
boot_p_tpr_fold_bpp = np.sort(boot_p_tpr_fold_bpp)
lq_p_tpr_fold_bpp = boot_p_tpr_fold_bpp[249]
hq_p_tpr_fold_bpp = boot_p_tpr_fold_bpp[749]

boot_p_tpr_drtr_bpp = data_i_p_boot[10]
boot_p_tpr_drtr_bpp = np.sort(boot_p_tpr_drtr_bpp)
lq_p_tpr_drtr_bpp = boot_p_tpr_drtr_bpp[249]
hq_p_tpr_drtr_bpp = boot_p_tpr_drtr_bpp[749]


boot_p_f_1_fold_mfe = data_i_p_boot[11]
boot_p_f_1_fold_mfe = np.sort(boot_p_f_1_fold_mfe)
lq_p_f_1_fold_mfe = boot_p_f_1_fold_mfe[249]
hq_p_f_1_fold_mfe = boot_p_f_1_fold_mfe[749]

boot_p_f_1_drtr_mfe = data_i_p_boot[12]
boot_p_f_1_drtr_mfe = np.sort(boot_p_f_1_drtr_mfe)
lq_p_f_1_drtr_mfe = boot_p_f_1_drtr_mfe[249]
hq_p_f_1_drtr_mfe = boot_p_f_1_drtr_mfe[749]

boot_p_f_1_drtr_prob = data_i_p_boot[13]
boot_p_f_1_drtr_prob = np.sort(boot_p_f_1_drtr_prob)
lq_p_f_1_drtr_prob = boot_p_f_1_drtr_prob[249]
hq_p_f_1_drtr_prob = boot_p_f_1_drtr_prob[749]

boot_p_f_1_fold_bpp = data_i_p_boot[14]
boot_p_f_1_fold_bpp = np.sort(boot_p_f_1_fold_bpp)
lq_p_f_1_fold_bpp = boot_p_f_1_fold_bpp[249]
hq_p_f_1_fold_bpp = boot_p_f_1_fold_bpp[749]

boot_p_f_1_drtr_bpp = data_i_p_boot[15]
boot_p_f_1_drtr_bpp = np.sort(boot_p_f_1_drtr_bpp)
lq_p_f_1_drtr_bpp = boot_p_f_1_drtr_bpp[249]
hq_p_f_1_drtr_bpp = boot_p_f_1_drtr_bpp[749]


boot_p_mcc_fold_mfe = data_i_p_boot[16]
boot_p_mcc_fold_mfe = np.sort(boot_p_mcc_fold_mfe)
lq_p_mcc_fold_mfe = boot_p_mcc_fold_mfe[249]
hq_p_mcc_fold_mfe = boot_p_mcc_fold_mfe[749]

boot_p_mcc_drtr_mfe = data_i_p_boot[17]
boot_p_mcc_drtr_mfe = np.sort(boot_p_mcc_drtr_mfe)
lq_p_mcc_drtr_mfe = boot_p_mcc_drtr_mfe[249]
hq_p_mcc_drtr_mfe = boot_p_mcc_drtr_mfe[749]

boot_p_mcc_drtr_prob = data_i_p_boot[18]
boot_p_mcc_drtr_prob = np.sort(boot_p_mcc_drtr_prob)
lq_p_mcc_drtr_prob = boot_p_mcc_drtr_prob[249]
hq_p_mcc_drtr_prob = boot_p_mcc_drtr_prob[749]

boot_p_mcc_fold_bpp = data_i_p_boot[19]
boot_p_mcc_fold_bpp = np.sort(boot_p_mcc_fold_bpp)
lq_p_mcc_fold_bpp = boot_p_mcc_fold_bpp[249]
hq_p_mcc_fold_bpp = boot_p_mcc_fold_bpp[749]

boot_p_mcc_drtr_bpp = data_i_p_boot[20]
boot_p_mcc_drtr_bpp = np.sort(boot_p_mcc_drtr_bpp)
lq_p_mcc_drtr_bpp = boot_p_mcc_drtr_bpp[249]
hq_p_mcc_drtr_bpp = boot_p_mcc_drtr_bpp[749]


### boot cp


boot_cp_ppv_fold_mfe = data_i_cp_boot[1]
boot_cp_ppv_fold_mfe = np.sort(boot_cp_ppv_fold_mfe)
lq_cp_ppv_fold_mfe = boot_cp_ppv_fold_mfe[249]
hq_cp_ppv_fold_mfe = boot_cp_ppv_fold_mfe[749]

boot_cp_ppv_drtr_mfe = data_i_cp_boot[2]
boot_cp_ppv_drtr_mfe = np.sort(boot_cp_ppv_drtr_mfe)
lq_cp_ppv_drtr_mfe = boot_cp_ppv_drtr_mfe[249]
hq_cp_ppv_drtr_mfe = boot_cp_ppv_drtr_mfe[749]

boot_cp_ppv_drtr_prob = data_i_cp_boot[3]
boot_cp_ppv_drtr_prob = np.sort(boot_cp_ppv_drtr_prob)
lq_cp_ppv_drtr_prob = boot_cp_ppv_drtr_prob[249]
hq_cp_ppv_drtr_prob = boot_cp_ppv_drtr_prob[749]

boot_cp_ppv_fold_bpp = data_i_cp_boot[4]
boot_cp_ppv_fold_bpp = np.sort(boot_cp_ppv_fold_bpp)
lq_cp_ppv_fold_bpp = boot_cp_ppv_fold_bpp[249]
hq_cp_ppv_fold_bpp = boot_cp_ppv_fold_bpp[749]

boot_cp_ppv_drtr_bpp = data_i_cp_boot[5]
boot_cp_ppv_drtr_bpp = np.sort(boot_cp_ppv_drtr_bpp)
lq_cp_ppv_drtr_bpp = boot_cp_ppv_drtr_bpp[249]
hq_cp_ppv_drtr_bpp = boot_cp_ppv_drtr_bpp[749]


boot_cp_tpr_fold_mfe = data_i_cp_boot[6]
boot_cp_tpr_fold_mfe = np.sort(boot_cp_tpr_fold_mfe)
lq_cp_tpr_fold_mfe = boot_cp_tpr_fold_mfe[249]
hq_cp_tpr_fold_mfe = boot_cp_tpr_fold_mfe[749]

boot_cp_tpr_drtr_mfe = data_i_cp_boot[7]
boot_cp_tpr_drtr_mfe = np.sort(boot_cp_tpr_drtr_mfe)
lq_cp_tpr_drtr_mfe = boot_cp_tpr_drtr_mfe[249]
hq_cp_tpr_drtr_mfe = boot_cp_tpr_drtr_mfe[749]

boot_cp_tpr_drtr_prob = data_i_cp_boot[8]
boot_cp_tpr_drtr_prob = np.sort(boot_cp_tpr_drtr_prob)
lq_cp_tpr_drtr_prob = boot_cp_tpr_drtr_prob[249]
hq_cp_tpr_drtr_prob = boot_cp_tpr_drtr_prob[749]

boot_cp_tpr_fold_bpp = data_i_cp_boot[9]
boot_cp_tpr_fold_bpp = np.sort(boot_cp_tpr_fold_bpp)
lq_cp_tpr_fold_bpp = boot_cp_tpr_fold_bpp[249]
hq_cp_tpr_fold_bpp = boot_cp_tpr_fold_bpp[749]

boot_cp_tpr_drtr_bpp = data_i_cp_boot[10]
boot_cp_tpr_drtr_bpp = np.sort(boot_cp_tpr_drtr_bpp)
lq_cp_tpr_drtr_bpp = boot_cp_tpr_drtr_bpp[249]
hq_cp_tpr_drtr_bpp = boot_cp_tpr_drtr_bpp[749]


boot_cp_f_1_fold_mfe = data_i_cp_boot[11]
boot_cp_f_1_fold_mfe = np.sort(boot_cp_f_1_fold_mfe)
lq_cp_f_1_fold_mfe = boot_cp_f_1_fold_mfe[249]
hq_cp_f_1_fold_mfe = boot_cp_f_1_fold_mfe[749]

boot_cp_f_1_drtr_mfe = data_i_cp_boot[12]
boot_cp_f_1_drtr_mfe = np.sort(boot_cp_f_1_drtr_mfe)
lq_cp_f_1_drtr_mfe = boot_cp_f_1_drtr_mfe[249]
hq_cp_f_1_drtr_mfe = boot_cp_f_1_drtr_mfe[749]

boot_cp_f_1_drtr_prob = data_i_cp_boot[13]
boot_cp_f_1_drtr_prob = np.sort(boot_cp_f_1_drtr_prob)
lq_cp_f_1_drtr_prob = boot_cp_f_1_drtr_prob[249]
hq_cp_f_1_drtr_prob = boot_cp_f_1_drtr_prob[749]

boot_cp_f_1_fold_bpp = data_i_cp_boot[14]
boot_cp_f_1_fold_bpp = np.sort(boot_cp_f_1_fold_bpp)
lq_cp_f_1_fold_bpp = boot_cp_f_1_fold_bpp[249]
hq_cp_f_1_fold_bpp = boot_cp_f_1_fold_bpp[749]

boot_cp_f_1_drtr_bpp = data_i_cp_boot[15]
boot_cp_f_1_drtr_bpp = np.sort(boot_cp_f_1_drtr_bpp)
lq_cp_f_1_drtr_bpp = boot_cp_f_1_drtr_bpp[249]
hq_cp_f_1_drtr_bpp = boot_cp_f_1_drtr_bpp[749]


boot_cp_mcc_fold_mfe = data_i_cp_boot[16]
boot_cp_mcc_fold_mfe = np.sort(boot_cp_mcc_fold_mfe)
lq_cp_mcc_fold_mfe = boot_cp_mcc_fold_mfe[249]
hq_cp_mcc_fold_mfe = boot_cp_mcc_fold_mfe[749]

boot_cp_mcc_drtr_mfe = data_i_cp_boot[17]
boot_cp_mcc_drtr_mfe = np.sort(boot_cp_mcc_drtr_mfe)
lq_cp_mcc_drtr_mfe = boot_cp_mcc_drtr_mfe[249]
hq_cp_mcc_drtr_mfe = boot_cp_mcc_drtr_mfe[749]

boot_cp_mcc_drtr_prob = data_i_cp_boot[18]
boot_cp_mcc_drtr_prob = np.sort(boot_cp_mcc_drtr_prob)
lq_cp_mcc_drtr_prob = boot_cp_mcc_drtr_prob[249]
hq_cp_mcc_drtr_prob = boot_cp_mcc_drtr_prob[749]

boot_cp_mcc_fold_bpp = data_i_cp_boot[19]
boot_cp_mcc_fold_bpp = np.sort(boot_cp_mcc_fold_bpp)
lq_cp_mcc_fold_bpp = boot_cp_mcc_fold_bpp[249]
hq_cp_mcc_fold_bpp = boot_cp_mcc_fold_bpp[749]

boot_cp_mcc_drtr_bpp = data_i_cp_boot[20]
boot_cp_mcc_drtr_bpp = np.sort(boot_cp_mcc_drtr_bpp)
lq_cp_mcc_drtr_bpp = boot_cp_mcc_drtr_bpp[249]
hq_cp_mcc_drtr_bpp = boot_cp_mcc_drtr_bpp[749]


### ========== MEDS
med_ppv_n_fold_mfe = boot_n_ppv_fold_mfe[499]
med_ppv_n_fold_bpp = boot_n_ppv_fold_bpp[499]
med_ppv_n_drtr_mfe = boot_n_ppv_drtr_mfe[499]
med_ppv_n_drtr_prob = boot_n_ppv_drtr_prob[499]
med_ppv_n_drtr_bpp = boot_n_ppv_drtr_bpp[499]


med_tpr_n_fold_mfe = boot_n_tpr_fold_mfe[499]
med_tpr_n_fold_bpp = boot_n_tpr_fold_bpp[499]
med_tpr_n_drtr_mfe = boot_n_tpr_drtr_mfe[499]
med_tpr_n_drtr_prob = boot_n_tpr_drtr_prob[499]
med_tpr_n_drtr_bpp = boot_n_tpr_drtr_bpp[499]


med_f_1_n_fold_mfe = boot_n_f_1_fold_mfe[499]
med_f_1_n_fold_bpp = boot_n_f_1_fold_bpp[499]
med_f_1_n_drtr_mfe = boot_n_f_1_drtr_mfe[499]
med_f_1_n_drtr_prob = boot_n_f_1_drtr_prob[499]
med_f_1_n_drtr_bpp = boot_n_f_1_drtr_bpp[499]


med_mcc_n_fold_mfe = boot_n_mcc_fold_mfe[499]
med_mcc_n_fold_bpp = boot_n_mcc_fold_bpp[499]
med_mcc_n_drtr_mfe = boot_n_mcc_drtr_mfe[499]
med_mcc_n_drtr_prob = boot_n_mcc_drtr_prob[499]
med_mcc_n_drtr_bpp = boot_n_mcc_drtr_bpp[499]


###
med_ppv_c_fold_mfe = boot_c_ppv_fold_mfe[499]
med_ppv_c_fold_bpp = boot_c_ppv_fold_bpp[499]
med_ppv_c_drtr_mfe = boot_c_ppv_drtr_mfe[499]
med_ppv_c_drtr_prob = boot_c_ppv_drtr_prob[499]
med_ppv_c_drtr_bpp = boot_c_ppv_drtr_bpp[499]


med_tpr_c_fold_mfe = boot_c_tpr_fold_mfe[499]
med_tpr_c_fold_bpp = boot_c_tpr_fold_bpp[499]
med_tpr_c_drtr_mfe = boot_c_tpr_drtr_mfe[499]
med_tpr_c_drtr_prob = boot_c_tpr_drtr_prob[499]
med_tpr_c_drtr_bpp = boot_c_tpr_drtr_bpp[499]


med_f_1_c_fold_mfe = boot_c_f_1_fold_mfe[499]
med_f_1_c_fold_bpp = boot_c_f_1_fold_bpp[499]
med_f_1_c_drtr_mfe = boot_c_f_1_drtr_mfe[499]
med_f_1_c_drtr_prob = boot_c_f_1_drtr_prob[499]
med_f_1_c_drtr_bpp = boot_c_f_1_drtr_bpp[499]


med_mcc_c_fold_mfe = boot_c_mcc_fold_mfe[499]
med_mcc_c_fold_bpp = boot_c_mcc_fold_bpp[499]
med_mcc_c_drtr_mfe = boot_c_mcc_drtr_mfe[499]
med_mcc_c_drtr_prob = boot_c_mcc_drtr_prob[499]
med_mcc_c_drtr_bpp = boot_c_mcc_drtr_bpp[499]


#
med_ppv_p_fold_mfe = boot_p_ppv_fold_mfe[499]
med_ppv_p_fold_bpp = boot_p_ppv_fold_bpp[499]
med_ppv_p_drtr_mfe = boot_p_ppv_drtr_mfe[499]
med_ppv_p_drtr_prob = boot_p_ppv_drtr_prob[499]
med_ppv_p_drtr_bpp = boot_p_ppv_drtr_bpp[499]


med_tpr_p_fold_mfe = boot_p_tpr_fold_mfe[499]
med_tpr_p_fold_bpp = boot_p_tpr_fold_bpp[499]
med_tpr_p_drtr_mfe = boot_p_tpr_drtr_mfe[499]
med_tpr_p_drtr_prob = boot_p_tpr_drtr_prob[499]
med_tpr_p_drtr_bpp = boot_p_tpr_drtr_bpp[499]


med_f_1_p_fold_mfe = boot_p_f_1_fold_mfe[499]
med_f_1_p_fold_bpp = boot_p_f_1_fold_bpp[499]
med_f_1_p_drtr_mfe = boot_p_f_1_drtr_mfe[499]
med_f_1_p_drtr_prob = boot_p_f_1_drtr_prob[499]
med_f_1_p_drtr_bpp = boot_p_f_1_drtr_bpp[499]


med_mcc_p_fold_mfe = boot_p_mcc_fold_mfe[499]
med_mcc_p_fold_bpp = boot_p_mcc_fold_bpp[499]
med_mcc_p_drtr_mfe = boot_p_mcc_drtr_mfe[499]
med_mcc_p_drtr_prob = boot_p_mcc_drtr_prob[499]
med_mcc_p_drtr_bpp = boot_p_mcc_drtr_bpp[499]


#
med_ppv_cp_fold_mfe = boot_cp_ppv_fold_mfe[499]
med_ppv_cp_fold_bpp = boot_cp_ppv_fold_bpp[499]
med_ppv_cp_drtr_mfe = boot_cp_ppv_drtr_mfe[499]
med_ppv_cp_drtr_prob = boot_cp_ppv_drtr_prob[499]
med_ppv_cp_drtr_bpp = boot_cp_ppv_drtr_bpp[499]


med_tpr_cp_fold_mfe = boot_cp_tpr_fold_mfe[499]
med_tpr_cp_fold_bpp = boot_cp_tpr_fold_bpp[499]
med_tpr_cp_drtr_mfe = boot_cp_tpr_drtr_mfe[499]
med_tpr_cp_drtr_prob = boot_cp_tpr_drtr_prob[499]
med_tpr_cp_drtr_bpp = boot_cp_tpr_drtr_bpp[499]


med_f_1_cp_fold_mfe = boot_cp_f_1_fold_mfe[499]
med_f_1_cp_fold_bpp = boot_cp_f_1_fold_bpp[499]
med_f_1_cp_drtr_mfe = boot_cp_f_1_drtr_mfe[499]
med_f_1_cp_drtr_prob = boot_cp_f_1_drtr_prob[499]
med_f_1_cp_drtr_bpp = boot_cp_f_1_drtr_bpp[499]


med_mcc_cp_fold_mfe = boot_cp_mcc_fold_mfe[499]
med_mcc_cp_fold_bpp = boot_cp_mcc_fold_bpp[499]
med_mcc_cp_drtr_mfe = boot_cp_mcc_drtr_mfe[499]
med_mcc_cp_drtr_prob = boot_cp_mcc_drtr_prob[499]
med_mcc_cp_drtr_bpp = boot_cp_mcc_drtr_bpp[499]




lq_wh_n_ppv_fold_mfe = boot_n_ppv_fold_mfe[24]
hq_wh_n_ppv_fold_mfe = boot_n_ppv_fold_mfe[974]

lq_wh_n_ppv_drtr_mfe = boot_n_ppv_drtr_mfe[24]
hq_wh_n_ppv_drtr_mfe = boot_n_ppv_drtr_mfe[974]

lq_wh_n_ppv_drtr_prob = boot_n_ppv_drtr_prob[24]
hq_wh_n_ppv_drtr_prob = boot_n_ppv_drtr_prob[974]

lq_wh_n_ppv_fold_bpp = boot_n_ppv_fold_bpp[24]
hq_wh_n_ppv_fold_bpp = boot_n_ppv_fold_bpp[974]

lq_wh_n_ppv_drtr_bpp = boot_n_ppv_drtr_bpp[24]
hq_wh_n_ppv_drtr_bpp = boot_n_ppv_drtr_bpp[974]


lq_wh_n_tpr_fold_mfe = boot_n_tpr_fold_mfe[24]
hq_wh_n_tpr_fold_mfe = boot_n_tpr_fold_mfe[974]

lq_wh_n_tpr_drtr_mfe = boot_n_tpr_drtr_mfe[24]
hq_wh_n_tpr_drtr_mfe = boot_n_tpr_drtr_mfe[974]

lq_wh_n_tpr_drtr_prob = boot_n_tpr_drtr_prob[24]
hq_wh_n_tpr_drtr_prob = boot_n_tpr_drtr_prob[974]

lq_wh_n_tpr_fold_bpp = boot_n_tpr_fold_bpp[24]
hq_wh_n_tpr_fold_bpp = boot_n_tpr_fold_bpp[974]

lq_wh_n_tpr_drtr_bpp = boot_n_tpr_drtr_bpp[24]
hq_wh_n_tpr_drtr_bpp = boot_n_tpr_drtr_bpp[974]


lq_wh_n_f_1_fold_mfe = boot_n_f_1_fold_mfe[24]
hq_wh_n_f_1_fold_mfe = boot_n_f_1_fold_mfe[974]

lq_wh_n_f_1_drtr_mfe = boot_n_f_1_drtr_mfe[24]
hq_wh_n_f_1_drtr_mfe = boot_n_f_1_drtr_mfe[974]

lq_wh_n_f_1_drtr_prob = boot_n_f_1_drtr_prob[24]
hq_wh_n_f_1_drtr_prob = boot_n_f_1_drtr_prob[974]

lq_wh_n_f_1_fold_bpp = boot_n_f_1_fold_bpp[24]
hq_wh_n_f_1_fold_bpp = boot_n_f_1_fold_bpp[974]

lq_wh_n_f_1_drtr_bpp = boot_n_f_1_drtr_bpp[24]
hq_wh_n_f_1_drtr_bpp = boot_n_f_1_drtr_bpp[974]


lq_wh_n_mcc_fold_mfe = boot_n_mcc_fold_mfe[24]
hq_wh_n_mcc_fold_mfe = boot_n_mcc_fold_mfe[974]

lq_wh_n_mcc_drtr_mfe = boot_n_mcc_drtr_mfe[24]
hq_wh_n_mcc_drtr_mfe = boot_n_mcc_drtr_mfe[974]

lq_wh_n_mcc_drtr_prob = boot_n_mcc_drtr_prob[24]
hq_wh_n_mcc_drtr_prob = boot_n_mcc_drtr_prob[974]

lq_wh_n_mcc_fold_bpp = boot_n_mcc_fold_bpp[24]
hq_wh_n_mcc_fold_bpp = boot_n_mcc_fold_bpp[974]

lq_wh_n_mcc_drtr_bpp = boot_n_mcc_drtr_bpp[24]
hq_wh_n_mcc_drtr_bpp = boot_n_mcc_drtr_bpp[974]

### boot c

lq_wh_c_ppv_fold_mfe = boot_c_ppv_fold_mfe[24]
hq_wh_c_ppv_fold_mfe = boot_c_ppv_fold_mfe[974]

lq_wh_c_ppv_drtr_mfe = boot_c_ppv_drtr_mfe[24]
hq_wh_c_ppv_drtr_mfe = boot_c_ppv_drtr_mfe[974]

lq_wh_c_ppv_drtr_prob = boot_c_ppv_drtr_prob[24]
hq_wh_c_ppv_drtr_prob = boot_c_ppv_drtr_prob[974]

lq_wh_c_ppv_fold_bpp = boot_c_ppv_fold_bpp[24]
hq_wh_c_ppv_fold_bpp = boot_c_ppv_fold_bpp[974]

lq_wh_c_ppv_drtr_bpp = boot_c_ppv_drtr_bpp[24]
hq_wh_c_ppv_drtr_bpp = boot_c_ppv_drtr_bpp[974]


lq_wh_c_tpr_fold_mfe = boot_c_tpr_fold_mfe[24]
hq_wh_c_tpr_fold_mfe = boot_c_tpr_fold_mfe[974]

lq_wh_c_tpr_drtr_mfe = boot_c_tpr_drtr_mfe[24]
hq_wh_c_tpr_drtr_mfe = boot_c_tpr_drtr_mfe[974]

lq_wh_c_tpr_drtr_prob = boot_c_tpr_drtr_prob[24]
hq_wh_c_tpr_drtr_prob = boot_c_tpr_drtr_prob[974]

lq_wh_c_tpr_fold_bpp = boot_c_tpr_fold_bpp[24]
hq_wh_c_tpr_fold_bpp = boot_c_tpr_fold_bpp[974]

lq_wh_c_tpr_drtr_bpp = boot_c_tpr_drtr_bpp[24]
hq_wh_c_tpr_drtr_bpp = boot_c_tpr_drtr_bpp[974]


lq_wh_c_f_1_fold_mfe = boot_c_f_1_fold_mfe[24]
hq_wh_c_f_1_fold_mfe = boot_c_f_1_fold_mfe[974]

lq_wh_c_f_1_drtr_mfe = boot_c_f_1_drtr_mfe[24]
hq_wh_c_f_1_drtr_mfe = boot_c_f_1_drtr_mfe[974]

lq_wh_c_f_1_drtr_prob = boot_c_f_1_drtr_prob[24]
hq_wh_c_f_1_drtr_prob = boot_c_f_1_drtr_prob[974]

lq_wh_c_f_1_fold_bpp = boot_c_f_1_fold_bpp[24]
hq_wh_c_f_1_fold_bpp = boot_c_f_1_fold_bpp[974]

lq_wh_c_f_1_drtr_bpp = boot_c_f_1_drtr_bpp[24]
hq_wh_c_f_1_drtr_bpp = boot_c_f_1_drtr_bpp[974]


lq_wh_c_mcc_fold_mfe = boot_c_mcc_fold_mfe[24]
hq_wh_c_mcc_fold_mfe = boot_c_mcc_fold_mfe[974]

lq_wh_c_mcc_drtr_mfe = boot_c_mcc_drtr_mfe[24]
hq_wh_c_mcc_drtr_mfe = boot_c_mcc_drtr_mfe[974]

lq_wh_c_mcc_drtr_prob = boot_c_mcc_drtr_prob[24]
hq_wh_c_mcc_drtr_prob = boot_c_mcc_drtr_prob[974]

lq_wh_c_mcc_fold_bpp = boot_c_mcc_fold_bpp[24]
hq_wh_c_mcc_fold_bpp = boot_c_mcc_fold_bpp[974]

lq_wh_c_mcc_drtr_bpp = boot_c_mcc_drtr_bpp[24]
hq_wh_c_mcc_drtr_bpp = boot_c_mcc_drtr_bpp[974]


### boot p

lq_wh_p_ppv_fold_mfe = boot_p_ppv_fold_mfe[24]
hq_wh_p_ppv_fold_mfe = boot_p_ppv_fold_mfe[974]

lq_wh_p_ppv_drtr_mfe = boot_p_ppv_drtr_mfe[24]
hq_wh_p_ppv_drtr_mfe = boot_p_ppv_drtr_mfe[974]

lq_wh_p_ppv_drtr_prob = boot_p_ppv_drtr_prob[24]
hq_wh_p_ppv_drtr_prob = boot_p_ppv_drtr_prob[974]

lq_wh_p_ppv_fold_bpp = boot_p_ppv_fold_bpp[24]
hq_wh_p_ppv_fold_bpp = boot_p_ppv_fold_bpp[974]

lq_wh_p_ppv_drtr_bpp = boot_p_ppv_drtr_bpp[24]
hq_wh_p_ppv_drtr_bpp = boot_p_ppv_drtr_bpp[974]


lq_wh_p_tpr_fold_mfe = boot_p_tpr_fold_mfe[24]
hq_wh_p_tpr_fold_mfe = boot_p_tpr_fold_mfe[974]

lq_wh_p_tpr_drtr_mfe = boot_p_tpr_drtr_mfe[24]
hq_wh_p_tpr_drtr_mfe = boot_p_tpr_drtr_mfe[974]

lq_wh_p_tpr_drtr_prob = boot_p_tpr_drtr_prob[24]
hq_wh_p_tpr_drtr_prob = boot_p_tpr_drtr_prob[974]

lq_wh_p_tpr_fold_bpp = boot_p_tpr_fold_bpp[24]
hq_wh_p_tpr_fold_bpp = boot_p_tpr_fold_bpp[974]

lq_wh_p_tpr_drtr_bpp = boot_p_tpr_drtr_bpp[24]
hq_wh_p_tpr_drtr_bpp = boot_p_tpr_drtr_bpp[974]


lq_wh_p_f_1_fold_mfe = boot_p_f_1_fold_mfe[24]
hq_wh_p_f_1_fold_mfe = boot_p_f_1_fold_mfe[974]

lq_wh_p_f_1_drtr_mfe = boot_p_f_1_drtr_mfe[24]
hq_wh_p_f_1_drtr_mfe = boot_p_f_1_drtr_mfe[974]

lq_wh_p_f_1_drtr_prob = boot_p_f_1_drtr_prob[24]
hq_wh_p_f_1_drtr_prob = boot_p_f_1_drtr_prob[974]

lq_wh_p_f_1_fold_bpp = boot_p_f_1_fold_bpp[24]
hq_wh_p_f_1_fold_bpp = boot_p_f_1_fold_bpp[974]

lq_wh_p_f_1_drtr_bpp = boot_p_f_1_drtr_bpp[24]
hq_wh_p_f_1_drtr_bpp = boot_p_f_1_drtr_bpp[974]


lq_wh_p_mcc_fold_mfe = boot_p_mcc_fold_mfe[24]
hq_wh_p_mcc_fold_mfe = boot_p_mcc_fold_mfe[974]

lq_wh_p_mcc_drtr_mfe = boot_p_mcc_drtr_mfe[24]
hq_wh_p_mcc_drtr_mfe = boot_p_mcc_drtr_mfe[974]

lq_wh_p_mcc_drtr_prob = boot_p_mcc_drtr_prob[24]
hq_wh_p_mcc_drtr_prob = boot_p_mcc_drtr_prob[974]

lq_wh_p_mcc_fold_bpp = boot_p_mcc_fold_bpp[24]
hq_wh_p_mcc_fold_bpp = boot_p_mcc_fold_bpp[974]

lq_wh_p_mcc_drtr_bpp = boot_p_mcc_drtr_bpp[24]
hq_wh_p_mcc_drtr_bpp = boot_p_mcc_drtr_bpp[974]


### boot cp


lq_wh_cp_ppv_fold_mfe = boot_cp_ppv_fold_mfe[24]
hq_wh_cp_ppv_fold_mfe = boot_cp_ppv_fold_mfe[974]

lq_wh_cp_ppv_drtr_mfe = boot_cp_ppv_drtr_mfe[24]
hq_wh_cp_ppv_drtr_mfe = boot_cp_ppv_drtr_mfe[974]

lq_wh_cp_ppv_drtr_prob = boot_cp_ppv_drtr_prob[24]
hq_wh_cp_ppv_drtr_prob = boot_cp_ppv_drtr_prob[974]

lq_wh_cp_ppv_fold_bpp = boot_cp_ppv_fold_bpp[24]
hq_wh_cp_ppv_fold_bpp = boot_cp_ppv_fold_bpp[974]

lq_wh_cp_ppv_drtr_bpp = boot_cp_ppv_drtr_bpp[24]
hq_wh_cp_ppv_drtr_bpp = boot_cp_ppv_drtr_bpp[974]


lq_wh_cp_tpr_fold_mfe = boot_cp_tpr_fold_mfe[24]
hq_wh_cp_tpr_fold_mfe = boot_cp_tpr_fold_mfe[974]

lq_wh_cp_tpr_drtr_mfe = boot_cp_tpr_drtr_mfe[24]
hq_wh_cp_tpr_drtr_mfe = boot_cp_tpr_drtr_mfe[974]

lq_wh_cp_tpr_drtr_prob = boot_cp_tpr_drtr_prob[24]
hq_wh_cp_tpr_drtr_prob = boot_cp_tpr_drtr_prob[974]

lq_wh_cp_tpr_fold_bpp = boot_cp_tpr_fold_bpp[24]
hq_wh_cp_tpr_fold_bpp = boot_cp_tpr_fold_bpp[974]

lq_wh_cp_tpr_drtr_bpp = boot_cp_tpr_drtr_bpp[24]
hq_wh_cp_tpr_drtr_bpp = boot_cp_tpr_drtr_bpp[974]


lq_wh_cp_f_1_fold_mfe = boot_cp_f_1_fold_mfe[24]
hq_wh_cp_f_1_fold_mfe = boot_cp_f_1_fold_mfe[974]

lq_wh_cp_f_1_drtr_mfe = boot_cp_f_1_drtr_mfe[24]
hq_wh_cp_f_1_drtr_mfe = boot_cp_f_1_drtr_mfe[974]

lq_wh_cp_f_1_drtr_prob = boot_cp_f_1_drtr_prob[24]
hq_wh_cp_f_1_drtr_prob = boot_cp_f_1_drtr_prob[974]

lq_wh_cp_f_1_fold_bpp = boot_cp_f_1_fold_bpp[24]
hq_wh_cp_f_1_fold_bpp = boot_cp_f_1_fold_bpp[974]

lq_wh_cp_f_1_drtr_bpp = boot_cp_f_1_drtr_bpp[24]
hq_wh_cp_f_1_drtr_bpp = boot_cp_f_1_drtr_bpp[974]


lq_wh_cp_mcc_fold_mfe = boot_cp_mcc_fold_mfe[24]
hq_wh_cp_mcc_fold_mfe = boot_cp_mcc_fold_mfe[974]

lq_wh_cp_mcc_drtr_mfe = boot_cp_mcc_drtr_mfe[24]
hq_wh_cp_mcc_drtr_mfe = boot_cp_mcc_drtr_mfe[974]

lq_wh_cp_mcc_drtr_prob = boot_cp_mcc_drtr_prob[24]
hq_wh_cp_mcc_drtr_prob = boot_cp_mcc_drtr_prob[974]

lq_wh_cp_mcc_fold_bpp = boot_cp_mcc_fold_bpp[24]
hq_wh_cp_mcc_fold_bpp = boot_cp_mcc_fold_bpp[974]

lq_wh_cp_mcc_drtr_bpp = boot_cp_mcc_drtr_bpp[24]
hq_wh_cp_mcc_drtr_bpp = boot_cp_mcc_drtr_bpp[974]


### =========================== DIFFS ================================= ###


means_n = [
        mean_mcc_n_fold_mfe[1],
        mean_mcc_n_drtr_mfe[1],
        mean_mcc_n_drtr_prob[1],
        mean_mcc_n_fold_bpp[1],
        mean_mcc_n_drtr_bpp[1],
        ]

means_c = [
        mean_mcc_c_fold_mfe[1],
        mean_mcc_c_drtr_mfe[1],
        mean_mcc_c_drtr_prob[1],
        mean_mcc_c_fold_bpp[1],
        mean_mcc_c_drtr_bpp[1],
        ]

means_p = [
        mean_mcc_p_fold_mfe[1],
        mean_mcc_p_drtr_mfe[1],
        mean_mcc_p_drtr_prob[1],
        mean_mcc_p_fold_bpp[1],
        mean_mcc_p_drtr_bpp[1],
        ]

means_cp = [
        mean_mcc_cp_fold_mfe[1],
        mean_mcc_cp_drtr_mfe[1],
        mean_mcc_cp_drtr_prob[1],
        mean_mcc_cp_fold_bpp[1],
        mean_mcc_cp_drtr_bpp[1],
        ]


diff_mean_c_mcc_fold_mfe = means_c[0] -  means_n[0]

diff_mean_c_mcc_drtr_mfe = means_c[1] -  means_n[1]

diff_mean_c_mcc_drtr_prob = means_c[2] -  means_n[2]

diff_mean_c_mcc_fold_bpp = means_c[3] -  means_n[3]

diff_mean_c_mcc_drtr_bpp = means_c[4] -  means_n[4]

diff_means_c_mcc = [
            diff_mean_c_mcc_drtr_bpp,
            diff_mean_c_mcc_drtr_prob,
            diff_mean_c_mcc_fold_bpp,
            diff_mean_c_mcc_drtr_mfe,
            diff_mean_c_mcc_fold_mfe,
            ]

diff_mean_p_mcc_fold_mfe = means_p[0] -  means_n[0]

diff_mean_p_mcc_drtr_mfe = means_p[1] -  means_n[1]

diff_mean_p_mcc_drtr_prob = means_p[2] -  means_n[2]

diff_mean_p_mcc_fold_bpp = means_p[3] -  means_n[3]

diff_mean_p_mcc_drtr_bpp = means_p[4] -  means_n[4]

diff_means_p_mcc = [
            diff_mean_p_mcc_drtr_bpp,
            diff_mean_p_mcc_drtr_prob,
            diff_mean_p_mcc_fold_bpp,
            diff_mean_p_mcc_drtr_mfe,
            diff_mean_p_mcc_fold_mfe,
            ]


diff_mean_cp_mcc_fold_mfe = means_cp[0] - means_n[0]

diff_mean_cp_mcc_drtr_mfe = means_cp[1] - means_n[1]

diff_mean_cp_mcc_drtr_prob = means_cp[2] - means_n[2]

diff_mean_cp_mcc_fold_bpp = means_cp[3] - means_n[3]

diff_mean_cp_mcc_drtr_bpp = means_cp[4] - means_n[4]

diff_means_cp_mcc = [
            diff_mean_cp_mcc_drtr_bpp,
            diff_mean_cp_mcc_drtr_prob,
            diff_mean_cp_mcc_fold_bpp,
            diff_mean_cp_mcc_drtr_mfe,
            diff_mean_cp_mcc_fold_mfe,
            ]

diff_med_c_mcc_fold_mfe = med_mcc_c_fold_mfe - med_mcc_n_fold_mfe
diff_med_c_mcc_fold_bpp = med_mcc_c_fold_bpp - med_mcc_n_fold_bpp
diff_med_c_mcc_drtr_mfe = med_mcc_c_drtr_mfe - med_mcc_n_drtr_mfe
diff_med_c_mcc_drtr_prob = med_mcc_c_drtr_prob - med_mcc_n_drtr_prob
diff_med_c_mcc_drtr_bpp = med_mcc_c_drtr_bpp - med_mcc_n_drtr_bpp


diff_med_cp_mcc_fold_mfe = med_mcc_cp_fold_mfe - med_mcc_n_fold_mfe
diff_med_cp_mcc_fold_bpp = med_mcc_cp_fold_bpp - med_mcc_n_fold_bpp
diff_med_cp_mcc_drtr_mfe = med_mcc_cp_drtr_mfe - med_mcc_n_drtr_mfe
diff_med_cp_mcc_drtr_prob = med_mcc_cp_drtr_prob - med_mcc_n_drtr_prob
diff_med_cp_mcc_drtr_bpp = med_mcc_cp_drtr_bpp - med_mcc_n_drtr_bpp


diff_med_p_mcc_fold_mfe = med_mcc_p_fold_mfe - med_mcc_n_fold_mfe
diff_med_p_mcc_fold_bpp = med_mcc_p_fold_bpp - med_mcc_n_fold_bpp
diff_med_p_mcc_drtr_mfe = med_mcc_p_drtr_mfe - med_mcc_n_drtr_mfe
diff_med_p_mcc_drtr_prob = med_mcc_p_drtr_prob - med_mcc_n_drtr_prob
diff_med_p_mcc_drtr_bpp = med_mcc_p_drtr_bpp - med_mcc_n_drtr_bpp



diff_boot_lq_c_mcc_fold_mfe = lq_c_mcc_fold_mfe - means_n[0]
diff_boot_hq_c_mcc_fold_mfe = hq_c_mcc_fold_mfe - means_n[0]

diff_boot_lq_c_mcc_drtr_mfe = lq_c_mcc_drtr_mfe - means_n[1]
diff_boot_hq_c_mcc_drtr_mfe = hq_c_mcc_drtr_mfe - means_n[1]

diff_boot_lq_c_mcc_drtr_prob = lq_c_mcc_drtr_prob - means_n[2]
diff_boot_hq_c_mcc_drtr_prob = hq_c_mcc_drtr_prob - means_n[2]

diff_boot_lq_c_mcc_fold_bpp = lq_c_mcc_fold_bpp - means_n[3]
diff_boot_hq_c_mcc_fold_bpp = hq_c_mcc_fold_bpp - means_n[3]

diff_boot_lq_c_mcc_drtr_bpp = lq_c_mcc_drtr_bpp - means_n[4]
diff_boot_hq_c_mcc_drtr_bpp = hq_c_mcc_drtr_bpp - means_n[4]


diff_boot_lq_wh_c_mcc_fold_mfe = lq_wh_c_mcc_fold_mfe - means_n[0]
diff_boot_hq_wh_c_mcc_fold_mfe = hq_wh_c_mcc_fold_mfe - means_n[0]

diff_boot_lq_wh_c_mcc_drtr_mfe = lq_wh_c_mcc_drtr_mfe - means_n[1]
diff_boot_hq_wh_c_mcc_drtr_mfe = hq_wh_c_mcc_drtr_mfe - means_n[1]

diff_boot_lq_wh_c_mcc_drtr_prob = lq_wh_c_mcc_drtr_prob - means_n[2]
diff_boot_hq_wh_c_mcc_drtr_prob = hq_wh_c_mcc_drtr_prob - means_n[2]

diff_boot_lq_wh_c_mcc_fold_bpp = lq_wh_c_mcc_fold_bpp - means_n[3]
diff_boot_hq_wh_c_mcc_fold_bpp = hq_wh_c_mcc_fold_bpp - means_n[3]

diff_boot_lq_wh_c_mcc_drtr_bpp = lq_wh_c_mcc_drtr_bpp - means_n[4]
diff_boot_hq_wh_c_mcc_drtr_bpp = hq_wh_c_mcc_drtr_bpp - means_n[4]


diff_boots_lq_c_mcc = [
            diff_boot_lq_c_mcc_drtr_bpp,
            diff_boot_lq_c_mcc_drtr_prob,
            diff_boot_lq_c_mcc_fold_bpp,
            diff_boot_lq_c_mcc_drtr_mfe,
            diff_boot_lq_c_mcc_fold_mfe
            ]

diff_boots_hq_c_mcc = [
            diff_boot_hq_c_mcc_drtr_bpp,
            diff_boot_hq_c_mcc_drtr_prob,
            diff_boot_hq_c_mcc_fold_bpp,
            diff_boot_hq_c_mcc_drtr_mfe,
            diff_boot_hq_c_mcc_fold_mfe
            ]


diff_boot_lq_p_mcc_fold_mfe = lq_p_mcc_fold_mfe - means_n[0]
diff_boot_hq_p_mcc_fold_mfe = hq_p_mcc_fold_mfe - means_n[0]

diff_boot_lq_p_mcc_drtr_mfe = lq_p_mcc_drtr_mfe - means_n[1]
diff_boot_hq_p_mcc_drtr_mfe = hq_p_mcc_drtr_mfe - means_n[1]

diff_boot_lq_p_mcc_drtr_prob = lq_p_mcc_drtr_prob - means_n[2]
diff_boot_hq_p_mcc_drtr_prob = hq_p_mcc_drtr_prob - means_n[2]

diff_boot_lq_p_mcc_fold_bpp = lq_p_mcc_fold_bpp - means_n[3]
diff_boot_hq_p_mcc_fold_bpp = hq_p_mcc_fold_bpp - means_n[3]

diff_boot_lq_p_mcc_drtr_bpp = lq_p_mcc_drtr_bpp - means_n[4]
diff_boot_hq_p_mcc_drtr_bpp = hq_p_mcc_drtr_bpp - means_n[4]


diff_boot_lq_wh_p_mcc_fold_mfe = lq_wh_p_mcc_fold_mfe - means_n[0]
diff_boot_hq_wh_p_mcc_fold_mfe = hq_wh_p_mcc_fold_mfe - means_n[0]

diff_boot_lq_wh_p_mcc_drtr_mfe = lq_wh_p_mcc_drtr_mfe - means_n[1]
diff_boot_hq_wh_p_mcc_drtr_mfe = hq_wh_p_mcc_drtr_mfe - means_n[1]

diff_boot_lq_wh_p_mcc_drtr_prob = lq_wh_p_mcc_drtr_prob - means_n[2]
diff_boot_hq_wh_p_mcc_drtr_prob = hq_wh_p_mcc_drtr_prob - means_n[2]

diff_boot_lq_wh_p_mcc_fold_bpp = lq_wh_p_mcc_fold_bpp - means_n[3]
diff_boot_hq_wh_p_mcc_fold_bpp = hq_wh_p_mcc_fold_bpp - means_n[3]

diff_boot_lq_wh_p_mcc_drtr_bpp = lq_wh_p_mcc_drtr_bpp - means_n[4]
diff_boot_hq_wh_p_mcc_drtr_bpp = hq_wh_p_mcc_drtr_bpp - means_n[4]


diff_boots_lq_p_mcc = [
            diff_boot_lq_p_mcc_drtr_bpp,
            diff_boot_lq_p_mcc_drtr_prob,
            diff_boot_lq_p_mcc_fold_bpp,
            diff_boot_lq_p_mcc_drtr_mfe,
            diff_boot_lq_p_mcc_fold_mfe
            ]

diff_boots_hq_p_mcc = [
            diff_boot_hq_p_mcc_drtr_bpp,
            diff_boot_hq_p_mcc_drtr_prob,
            diff_boot_hq_p_mcc_fold_bpp,
            diff_boot_hq_p_mcc_drtr_mfe,
            diff_boot_hq_p_mcc_fold_mfe
            ]


diff_boot_lq_cp_mcc_fold_mfe = lq_cp_mcc_fold_mfe - means_n[0]
diff_boot_hq_cp_mcc_fold_mfe = hq_cp_mcc_fold_mfe - means_n[0]

diff_boot_lq_cp_mcc_drtr_mfe = lq_cp_mcc_drtr_mfe - means_n[1]
diff_boot_hq_cp_mcc_drtr_mfe = hq_cp_mcc_drtr_mfe - means_n[1]

diff_boot_lq_cp_mcc_drtr_prob = lq_cp_mcc_drtr_prob - means_n[2]
diff_boot_hq_cp_mcc_drtr_prob = hq_cp_mcc_drtr_prob - means_n[2]

diff_boot_lq_cp_mcc_fold_bpp = lq_cp_mcc_fold_bpp - means_n[3]
diff_boot_hq_cp_mcc_fold_bpp = hq_cp_mcc_fold_bpp - means_n[3]

diff_boot_lq_cp_mcc_drtr_bpp = lq_cp_mcc_drtr_bpp - means_n[4]
diff_boot_hq_cp_mcc_drtr_bpp = hq_cp_mcc_drtr_bpp - means_n[4]


diff_boot_lq_wh_cp_mcc_fold_mfe = lq_wh_cp_mcc_fold_mfe - means_n[0]
diff_boot_hq_wh_cp_mcc_fold_mfe = hq_wh_cp_mcc_fold_mfe - means_n[0]

diff_boot_lq_wh_cp_mcc_drtr_mfe = lq_wh_cp_mcc_drtr_mfe - means_n[1]
diff_boot_hq_wh_cp_mcc_drtr_mfe = hq_wh_cp_mcc_drtr_mfe - means_n[1]

diff_boot_lq_wh_cp_mcc_drtr_prob = lq_wh_cp_mcc_drtr_prob - means_n[2]
diff_boot_hq_wh_cp_mcc_drtr_prob = hq_wh_cp_mcc_drtr_prob - means_n[2]

diff_boot_lq_wh_cp_mcc_fold_bpp = lq_wh_cp_mcc_fold_bpp - means_n[3]
diff_boot_hq_wh_cp_mcc_fold_bpp = hq_wh_cp_mcc_fold_bpp - means_n[3]

diff_boot_lq_wh_cp_mcc_drtr_bpp = lq_wh_cp_mcc_drtr_bpp - means_n[4]
diff_boot_hq_wh_cp_mcc_drtr_bpp = hq_wh_cp_mcc_drtr_bpp - means_n[4]


diff_boots_lq_cp_mcc = [
            diff_boot_lq_cp_mcc_drtr_bpp,
            diff_boot_lq_cp_mcc_drtr_prob,
            diff_boot_lq_cp_mcc_fold_bpp,
            diff_boot_lq_cp_mcc_drtr_mfe,
            diff_boot_lq_cp_mcc_fold_mfe
            ]

diff_boots_hq_cp_mcc = [
            diff_boot_hq_cp_mcc_drtr_bpp,
            diff_boot_hq_cp_mcc_drtr_prob,
            diff_boot_hq_cp_mcc_fold_bpp,
            diff_boot_hq_cp_mcc_drtr_mfe,
            diff_boot_hq_cp_mcc_fold_mfe
            ]

indexer = []
for i in range(0, len(means_n)):
    indexer.append(i)
# print(indexer)
# print(diff_mean_c_mcc_fold_mfe)

indexer = [
        "drtr_bpp",
        "drtr_prob",
        "fold_bpp",
        "drtr_mfe",
        "fold_mfe"
        ]

if args.graphics==True:
    item_1 = {}
    # item_1["label"] ='RNAfold:MFE/cp'
    # item_1["mean"] = diff_mean_cp_mcc_fold_mfe
    item_1["med"] = diff_med_cp_mcc_fold_mfe
    item_1["q1"] = diff_boot_lq_cp_mcc_fold_mfe
    item_1["q3"] = diff_boot_hq_cp_mcc_fold_mfe
    item_1["whislo"] = diff_boot_lq_wh_cp_mcc_fold_mfe # required
    item_1["whishi"] = diff_boot_hq_wh_cp_mcc_fold_mfe # required
    item_1["fliers"] = []
    stats_1 = [item_1]

    item_2 = {}
    # item_2["label"] ='RNAfold:BPP/cp'
    item_2["med"] = diff_med_cp_mcc_fold_bpp
    item_2["q1"] = diff_boot_lq_cp_mcc_fold_bpp
    item_2["q3"] = diff_boot_hq_cp_mcc_fold_bpp
    item_2["whislo"] = diff_boot_lq_wh_cp_mcc_fold_bpp # required
    item_2["whishi"] = diff_boot_hq_wh_cp_mcc_fold_bpp # required
    item_2["fliers"] = []
    stats_2 = [item_2]

    item_3 = {}
    # item_3["label"] ='DrTr:MFE/cp'
    item_3["med"] = diff_med_cp_mcc_drtr_mfe
    item_3["q1"] = diff_boot_lq_cp_mcc_drtr_mfe
    item_3["q3"] = diff_boot_hq_cp_mcc_drtr_mfe
    item_3["whislo"] = diff_boot_lq_wh_cp_mcc_drtr_mfe # required
    item_3["whishi"] = diff_boot_hq_wh_cp_mcc_drtr_mfe # required
    item_3["fliers"] = []
    stats_3 = [item_3]

    item_4 = {}
    # item_4["label"] ='DrTr:PROB/cp'
    item_4["med"] = diff_med_cp_mcc_drtr_prob
    item_4["q1"] = diff_boot_lq_cp_mcc_drtr_prob
    item_4["q3"] = diff_boot_hq_cp_mcc_drtr_prob
    item_4["whislo"] = diff_boot_lq_wh_cp_mcc_drtr_prob # required
    item_4["whishi"] = diff_boot_hq_wh_cp_mcc_drtr_prob # required
    item_4["fliers"] = []
    stats_4 = [item_4]

    item_5 = {}
    # item_5["label"] ='DrTr:BPP/cp'
    item_5["med"] = diff_med_cp_mcc_drtr_bpp
    item_5["q1"] = diff_boot_lq_cp_mcc_drtr_bpp
    item_5["q3"] = diff_boot_hq_cp_mcc_drtr_bpp
    item_5["whislo"] = diff_boot_lq_wh_cp_mcc_drtr_bpp # required
    item_5["whishi"] = diff_boot_hq_wh_cp_mcc_drtr_bpp # required
    item_5["fliers"] = []
    stats_5 = [item_5]

    item_6 = {}
    # item_6["label"] ='RNAfold:MFE/c'
    item_6["med"] = diff_med_c_mcc_fold_mfe
    item_6["q1"] = diff_boot_lq_c_mcc_fold_mfe
    item_6["q3"] = diff_boot_hq_c_mcc_fold_mfe
    item_6["whislo"] = diff_boot_lq_wh_c_mcc_fold_mfe # required
    item_6["whishi"] = diff_boot_hq_wh_c_mcc_fold_mfe # required
    item_6["fliers"] = []
    stats_6 = [item_6]

    item_7 = {}
    # item_7["label"] ='RNAfold:BPP/c'
    item_7["med"] = diff_med_c_mcc_fold_bpp
    item_7["q1"] = diff_boot_lq_c_mcc_fold_bpp
    item_7["q3"] = diff_boot_hq_c_mcc_fold_bpp
    item_7["whislo"] = diff_boot_lq_wh_c_mcc_fold_bpp # required
    item_7["whishi"] = diff_boot_hq_wh_c_mcc_fold_bpp # required
    item_7["fliers"] = []
    stats_7 = [item_7]

    item_8 = {}
    # item_8["label"] ='DrTr:MFE/c'
    item_8["med"] = diff_med_c_mcc_drtr_mfe
    item_8["q1"] = diff_boot_lq_c_mcc_drtr_mfe
    item_8["q3"] = diff_boot_hq_c_mcc_drtr_mfe
    item_8["whislo"] = diff_boot_lq_wh_c_mcc_drtr_mfe # required
    item_8["whishi"] = diff_boot_hq_wh_c_mcc_drtr_mfe # required
    item_8["fliers"] = []
    stats_8 = [item_8]

    item_9 = {}
    # item_9["label"] ='DrTr:PROB/c'
    item_9["med"] = diff_med_c_mcc_drtr_prob
    item_9["q1"] = diff_boot_lq_c_mcc_drtr_prob
    item_9["q3"] = diff_boot_hq_c_mcc_drtr_prob
    item_9["whislo"] = diff_boot_lq_wh_c_mcc_drtr_prob # required
    item_9["whishi"] = diff_boot_hq_wh_c_mcc_drtr_prob # required
    item_9["fliers"] = []
    stats_9 = [item_9]

    item_10 = {}
    # item_10["label"] ='DrTr:BPP/c'
    item_10["med"] = diff_med_p_mcc_drtr_bpp
    item_10["q1"] = diff_boot_lq_p_mcc_drtr_bpp
    item_10["q3"] = diff_boot_hq_p_mcc_drtr_bpp
    item_10["whislo"] = diff_boot_lq_wh_p_mcc_drtr_bpp # required
    item_10["whishi"] = diff_boot_hq_wh_p_mcc_drtr_bpp # required
    item_10["fliers"] = []
    stats_10 = [item_10]

    item_11 = {}
    # item_11["label"] ='RNAfold:MFE/p'
    item_11["med"] = diff_med_p_mcc_fold_mfe
    item_11["q1"] = diff_boot_lq_p_mcc_fold_mfe
    item_11["q3"] = diff_boot_hq_p_mcc_fold_mfe
    item_11["whislo"] = diff_boot_lq_wh_p_mcc_fold_mfe # required
    item_11["whishi"] = diff_boot_hq_wh_p_mcc_fold_mfe # required
    item_11["fliers"] = []
    stats_11 = [item_11]

    item_12 = {}
    # item_12["label"] ='RNAfold:BPP/p'
    item_12["med"] = diff_med_p_mcc_fold_bpp
    item_12["q1"] = diff_boot_lq_p_mcc_fold_bpp
    item_12["q3"] = diff_boot_hq_p_mcc_fold_bpp
    item_12["whislo"] = diff_boot_lq_wh_p_mcc_fold_bpp # required
    item_12["whishi"] = diff_boot_hq_wh_p_mcc_fold_bpp # required
    item_12["fliers"] = []
    stats_12 = [item_12]

    item_13 = {}
    # item_13["label"] ='DrTr:MFE/p'
    item_13["med"] = diff_med_p_mcc_drtr_mfe
    item_13["q1"] = diff_boot_lq_p_mcc_drtr_mfe
    item_13["q3"] = diff_boot_hq_p_mcc_drtr_mfe
    item_13["whislo"] = diff_boot_lq_wh_p_mcc_drtr_mfe # required
    item_13["whishi"] = diff_boot_hq_wh_p_mcc_drtr_mfe # required
    item_13["fliers"] = []
    stats_13 = [item_13]

    item_14 = {}
    # item_14["label"] ='DrTr:PROB/p'
    item_14["med"] = diff_med_p_mcc_drtr_prob
    item_14["q1"] = diff_boot_lq_p_mcc_drtr_prob
    item_14["q3"] = diff_boot_hq_p_mcc_drtr_prob
    item_14["whislo"] = diff_boot_lq_wh_p_mcc_drtr_prob # required
    item_14["whishi"] = diff_boot_hq_wh_p_mcc_drtr_prob # required
    item_14["fliers"] = []
    stats_14 = [item_14]

    item_15 = {}
    # item_15["label"] ='DrTr:BPP/p'
    item_15["med"] = diff_med_p_mcc_drtr_bpp
    item_15["q1"] = diff_boot_lq_p_mcc_drtr_bpp
    item_15["q3"] = diff_boot_hq_p_mcc_drtr_bpp
    item_15["whislo"] = diff_boot_lq_wh_p_mcc_drtr_bpp # required
    item_15["whishi"] = diff_boot_hq_wh_p_mcc_drtr_bpp # required
    item_15["fliers"] = []
    stats_15 = [item_15]
    fig, ax = plt.subplots()
    # ax.xticks("2")

    # ax.plot(names)
    # ax.plot(indexer, diff_boots_lq_cp_mcc, color = 'red', linewidth=0.5, alpha=0.5)
    # ax.plot(indexer, diff_boots_hq_cp_mcc, color = 'red', linewidth=0.5, alpha=0.5)
    # plt.axis([0, 4, 0, 0.04])
    # plt.yticks(["x"])
    # ax.set_xlabel('Algorithm', fontsize = 12)
    # ax.set_ylabel('MCC difference', fontsize = 12)
    # ax.set_title('Effect of removing non canonical bps and pseudoknots \n on performance based MCC', fontsize=14)
    # plt.legend()

    ax.bxp(stats_1, positions=[12], widths=0.50)
    ax.bxp(stats_2, positions=[13], widths=0.50)
    ax.bxp(stats_3, positions=[14], widths=0.50)
    ax.bxp(stats_4, positions=[15], widths=0.50)
    ax.bxp(stats_5, positions=[16], widths=0.50)
    ax.bxp(stats_6, positions=[0], widths=0.50)
    ax.bxp(stats_7, positions=[1], widths=0.50)
    ax.bxp(stats_8, positions=[2], widths=0.50)
    ax.bxp(stats_9, positions=[3], widths=0.50)
    ax.bxp(stats_10, positions=[4], widths=0.50)
    ax.bxp(stats_11, positions=[6], widths=0.50)
    ax.bxp(stats_12, positions=[7], widths=0.50)
    ax.bxp(stats_13, positions=[8], widths=0.50)
    ax.bxp(stats_14, positions=[9], widths=0.50)
    ax.bxp(stats_15, positions=[10], widths=0.50)
    # ax.set_xticks("2","2")
    values=[0,1,2,3,4,6,7,8,9,10,12,13,14,15,16]
    values_y=[
              diff_mean_cp_mcc_fold_mfe,
              diff_mean_cp_mcc_fold_bpp,
              diff_mean_cp_mcc_drtr_mfe,
              diff_mean_cp_mcc_drtr_prob,
              diff_mean_cp_mcc_drtr_bpp,
              diff_mean_c_mcc_fold_mfe,
              diff_mean_c_mcc_fold_bpp,
              diff_mean_c_mcc_drtr_mfe,
              diff_mean_c_mcc_drtr_prob,
              diff_mean_c_mcc_drtr_bpp,
              diff_mean_p_mcc_fold_mfe,
              diff_mean_p_mcc_fold_bpp,
              diff_mean_p_mcc_drtr_mfe,
              diff_mean_p_mcc_drtr_prob,
              diff_mean_p_mcc_drtr_bpp
              ]
    names=["fold_mfe:cp","fold_bpp:cp","DrTr_mfe:cp","DrTr_prob:cp","DrTr_bpp:cp",
    "RNAfold_mfe:c","RNAfold_bpp:c","DrTr_mfe:c","DrTr_prob:c","DrTr_bpp:c",
    "RNAfold_mfe:p","RNAfold_bpp:p","DrTr_mfe:p","DrTr_prob:p","DrTr_bpp:p"]
    # plt.scatter(values,values_y, marker="_", color="red")
    # plt.xticks(names)
    plt.gca().set_xticklabels(names, rotation=55, ha="right", rotation_mode="anchor")
    plt.gcf().subplots_adjust(bottom=0.28, left=0.15, right=0.94, top=0.88)
    if label=="m":
        x = "modified"
    elif label=="mu":
        x = "RTunblock"
    elif label=="n":
        x = "RNAf:n; Drtr:m"
    label = "Effect of removing non canonical base pairs \nand pseudoknots from reference structure"
    plt.title(label)
    plt.xlabel("Algorithm")
    plt.ylabel("deltaMCC")
    plt.show()
    # ax.scatter(x[4],y[4], color = 'red', label='drtr: bpp')
    #
    # ax.plot(x[4],lq_tpr_drtr_bpp, color = 'red', linestyle='dashed', marker = '_' , markersize=25)
    # ax.plot(x[4],hq_tpr_drtr_bpp, color = 'red', linestyle='dashed', marker = '_' , markersize=25)
    # ax.plot([x[4],x[4]], [lq_tpr_drtr_bpp,hq_tpr_drtr_bpp], color = 'red', linewidth=1)
    #
    # ax.plot(lq_ppv_drtr_bpp, y[4], color = 'red', linestyle='dashed', marker = '|' , markersize=25)
    # ax.plot(hq_ppv_drtr_bpp, y[4], color = 'red', linestyle='dashed', marker = '|' , markersize=25)
    # ax.plot([lq_ppv_drtr_bpp,hq_ppv_drtr_bpp],[y[4],y[4]], color = 'red', linewidth=1)
    #
    # ax.set_ylabel('TPR / Sensitivity')
    # ax.set_xlabel('PPV / Precision')
    # ax.set_title('PPV vs TPR over whole set')
    # plt.legend()
    # plt.show()
    # fig = plt.figure()
    # ax = fig.add_subplot()
    # ax.plot(indexer, diff_means_c_mcc, color = 'green', label='canonical bps removed')
    # ax.fill_between(indexer, diff_boots_lq_c_mcc, diff_boots_hq_c_mcc, color = 'green', alpha=0.1)
    # # ax.plot(indexer, diff_boots_lq_c_mcc, color = 'green', linewidth=0.5, alpha=0.5)
    # # ax.plot(indexer, diff_boots_hq_c_mcc, color = 'green', linewidth=0.5, alpha=0.5)
    #
    # ax.plot(indexer, diff_means_p_mcc, color = 'blue', label='pseudoknots removed')
    # ax.fill_between(indexer, diff_boots_lq_p_mcc, diff_boots_hq_p_mcc, color = 'blue', alpha=0.1)
    # # ax.plot( color = 'blue', linewidth=0.5, alpha=0.5)
    #
    # ax.plot(indexer, diff_means_cp_mcc, color = 'red', label='canonical bps and pseudoknots removed')
    # ax.fill_between(indexer, diff_boots_lq_cp_mcc, diff_boots_hq_cp_mcc, color = 'red', alpha=0.1)
    # #
    # # print(lq_tpr_fold_bpp)
    # # print(hq_tpr_fold_bpp)
    # # print(mean_tpr_fold_bpp[1])
    #
    # item_2 = {}
    # item_2["label"] ='drtr_mfe_tpr'
    # item_2["med"] = m_drtr_tpr_mfe[1]
    # item_2["q1"] = m_drtr_tpr_mfe[1]
    # item_2["q3"] = m_drtr_tpr_mfe[1]
    # item_2["whislo"] = lq_drtr_tpr_mfe # required
    # item_2["whishi"] = hq_drtr_tpr_mfe # required
    # item_2["fliers"] = []
    # stats_2 = [item_2]
    #
    #
    # # fig2, axes2 = plt.subplots(1, 1)
    # axes.bxp(stats, showbox=False)
    # # axes2.bxp(stats2, showbox=False)
    # axes.bxp(stats_2, positions=[1])
    # axes.bxp(stats_3, positions=[2])
    # axes.bxp(stats_4, positions=[3])
    # axes.bxp(stats_5, positions=[4])
    # axes.bxp(stats_6, positions=[5])
    # axes.bxp(stats_7, positions=[6])
    # axes.bxp(stats_8, positions=[7])
    # axes.bxp(stats_9, positions=[8])
    # axes.bxp(stats_10, positions=[9])
    # axes.bxp(stats_11, positions=[10])
    # axes.bxp(stats_12, positions=[11])
    # axes.bxp(stats_13, positions=[12])
    # axes.bxp(stats_14, positions=[13])
    # axes.bxp(stats_15, positions=[14])
    # axes.bxp(stats_16, positions=[15])
    # axes.bxp(stats_17, positions=[16])
    # axes.bxp(stats_18, positions=[17])
    # axes.bxp(stats_19, positions=[18])
    # axes.bxp(stats_20, positions=[19])
    #
    # plt.xticks(rotation=55,
    #             ha="right",
    #             rotation_mode="anchor",
    #             fontsize=12)
    # axes.set_title('bootstrapping')
    # y_label = "mcc"
    # plt.show()
