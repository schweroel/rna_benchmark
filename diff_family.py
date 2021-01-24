#!python
###=========================import modules ======================###
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
from numpy import genfromtxt
import matplotlib.pyplot as plt

###========================settings for pandas and numpy========###
pd.set_option('display.max_rows', 10000) #some default pandas settings will chop data
pd.set_option('display.max_columns', 10000)
pd.set_option('display.width', 10000)
np.set_printoptions(threshold=sys.maxsize) #if printed out, matrix won't be truncated


parser = argparse.ArgumentParser(description='Compare ct data with RNAfold')
parser.add_argument('--test','-t', action='store_true',
                    default=False,
                    help='run with testseq')

parser.add_argument('--input','-i', action='store',
                    default="comp_m_n_fin_f_mean.csv",
                    help='input filename')

parser.add_argument('--label','-l', action='store',
                    default="label missing",
                    help='input a label')

parser.add_argument('--graphics','-g', action='store_true',
                    default=False,
                    help='produce plot')

args=parser.parse_args()


### ==================================== INPUT DATA ====================================== ###
row_names = ["fold_mfe_ppv", "fold_bpp_ppv", "drtr_mfe_ppv", "drtr_prob_ppv", "drtr_bpp_ppv",
             "fold_mfe_tpr", "fold_bpp_tpr", "drtr_mfe_tpr", "drtr_prob_tpr", "drtr_bpp_tpr",
             "fold_mfe_f1", "fold_bpp_f1", "drtr_mfe_f1", "drtr_prob_f1", "drtr_bpp_f1",
             "fold_mfe_mcc", "fold_bpp_mcc", "drtr_mfe_mcc", "drtr_prob_mcc", "drtr_bpp_mcc"]

input_name_n = args.input
input_name_cp = input_name_n [:7] + "cp" + input_name_n[8:]
print(input_name_n)
print(input_name_n[6:7])
if input_name_n[6:7]=="u":
    input_name_cp = input_name_n [:8] + "cp" + input_name_n[9:]
print(input_name_cp)
data_raw = pd.read_csv(input_name_n)
data_raw_cp = pd.read_csv(input_name_cp)


# data_raw = genfromtxt(input_name_n, delimiter = ",")
# print("x")
# print(data_raw[1])


family_dict = (
            "16SrRNA",
            "23SrRNA",
            "5SrRNA",
            "Cili.Telo. RNA",
            "Cis-reg.element",
            "GIIIntron",
            "GIIntron",
            "Ham.Ribozyme",
            "HDVRibozyme",
            "IRES",
            "OtherRibozyme",
            "OtherRNA",
            "OtherrRNA",
            "RNAIII",
            "RNaseE5UTR",
            "RNaseMRPRNA",
            "RNasePRNA",
            "snRNA",
            "SRPRNA",
            "SyntheticRNA",
            "tmRNA",
            "tRNA",
            "Vert.Telo. RNA",
            "Viral&Phage",
            "YRNA" )

for i in range(0,len(family_dict)):
    create_names = list()

# print(data_raw.iloc[0])
# print(data_raw)
# fold_mfe_ppv = data_raw.iloc[0]


fold_mfe_ppv = data_raw.iloc[0][1:]
fold_bpp_ppv = data_raw.iloc[1][1:]
drtr_mfe_ppv = data_raw.iloc[2][1:]
drtr_prob_ppv = data_raw.iloc[3][1:]
drtr_bpp_ppv = data_raw.iloc[4][1:]


fold_mfe_tpr = data_raw.iloc[5][1:]
fold_bpp_tpr = data_raw.iloc[6][1:]
drtr_mfe_tpr = data_raw.iloc[7][1:]
drtr_prob_tpr = data_raw.iloc[8][1:]
drtr_bpp_tpr = data_raw.iloc[9][1:]


fold_mfe_f1 = data_raw.iloc[10][1:]
fold_bpp_f1 = data_raw.iloc[11][1:]
drtr_mfe_f1 = data_raw.iloc[12][1:]
drtr_prob_f1 = data_raw.iloc[13][1:]
drtr_bpp_f1 = data_raw.iloc[14][1:]


fold_mfe_mcc = data_raw.iloc[15][1:]
fold_bpp_mcc = data_raw.iloc[16][1:]
drtr_mfe_mcc = data_raw.iloc[17][1:]
drtr_prob_mcc = data_raw.iloc[18][1:]
drtr_bpp_mcc = data_raw.iloc[19][1:]



fold_mfe_ppv_cp = data_raw_cp.iloc[0][1:]
fold_bpp_ppv_cp = data_raw_cp.iloc[1][1:]
drtr_mfe_ppv_cp = data_raw_cp.iloc[2][1:]
drtr_prob_ppv_cp = data_raw_cp.iloc[3][1:]
drtr_bpp_ppv_cp = data_raw_cp.iloc[4][1:]


fold_mfe_tpr_cp = data_raw_cp.iloc[5][1:]
fold_bpp_tpr_cp = data_raw_cp.iloc[6][1:]
drtr_mfe_tpr_cp = data_raw_cp.iloc[7][1:]
drtr_prob_tpr_cp = data_raw_cp.iloc[8][1:]
drtr_bpp_tpr_cp = data_raw_cp.iloc[9][1:]


fold_mfe_f1_cp = data_raw_cp.iloc[10][1:]
fold_bpp_f1_cp = data_raw_cp.iloc[11][1:]
drtr_mfe_f1_cp = data_raw_cp.iloc[12][1:]
drtr_prob_f1_cp = data_raw_cp.iloc[13][1:]
drtr_bpp_f1_cp = data_raw_cp.iloc[14][1:]


fold_mfe_mcc_cp = data_raw_cp.iloc[15][1:]
fold_bpp_mcc_cp = data_raw_cp.iloc[16][1:]
drtr_mfe_mcc_cp = data_raw_cp.iloc[17][1:]
drtr_prob_mcc_cp = data_raw_cp.iloc[18][1:]
drtr_bpp_mcc_cp = data_raw_cp.iloc[19][1:]

if args.graphics==True:

    titel = args.label
    barwidth = 0.18

    r_1 = np.arange(len(fold_mfe_mcc))
    r_2 = [x + barwidth for x in r_1]
    r_3 = [x + barwidth for x in r_2]
    r_4 = [x + barwidth for x in r_3]
    r_5 = [x + barwidth for x in r_4]

    plt.bar(r_1, fold_mfe_mcc_cp, width=barwidth, color='#d9d9d9', label='non canonical base pairs and pseudoknots removed from reference structure')
    plt.bar(r_2, drtr_mfe_mcc_cp, width=barwidth, color='#d9d9d9')
    plt.bar(r_3, drtr_prob_mcc_cp, width=barwidth, color='#d9d9d9')
    plt.bar(r_4, fold_bpp_mcc_cp, width=barwidth, color='#d9d9d9')
    plt.bar(r_5, drtr_bpp_mcc_cp, width=barwidth, color='#d9d9d9')

#0F7A26
    plt.bar(r_1, fold_mfe_mcc, width=barwidth, color='#173F5F', label='RNAfold: mfe')
    plt.bar(r_2, drtr_mfe_mcc, width=barwidth, color='#3CAEA3', label='DrTransformer: mfe')
    plt.bar(r_3, drtr_prob_mcc, width=barwidth, color='#0F7A26', label='DrTransformer: prob')
    plt.bar(r_4, fold_bpp_mcc, width=barwidth, color='#FF7A26', label='RNAfold: bpp')
    plt.bar(r_5, drtr_bpp_mcc, width=barwidth, color='#FFFF38', label='DrTransformer: bpp')


    plt.title(titel, fontsize=20, y=1.03)
    plt.xlabel('Family', fontsize=16)
    plt.ylabel('MCC',  fontsize=16)
    plt.axis([-0.4, 25, 0, 1])
    plt.xticks([r + barwidth for r in range (len(fold_mfe_mcc))],
                                                   [family_dict[0] + "(n=" + "42" + ")",
                                                   family_dict[5] + "(n=" + "12" + ")",
                                                   family_dict[15] + "(n=" + "5" + ")",
                                                   family_dict[8] + "(n=" + "7" + ")",
                                                   family_dict[20] + "(n=" + "657" + ")",
                                                   family_dict[6] + "(n=" + "112" + ")",
                                                   family_dict[10] + "(n=" + "14" + ")",
                                                   family_dict[3] + "(n=" + "18" + ")",
                                                   family_dict[16] + "(n=" + "437" + ")",
                                                   family_dict[22] + "(n=" + "6" + ")",
                                                   family_dict[7] + "(n=" + "136" + ")",
                                                   family_dict[18] + "(n=" + "375" + ")",
                                                   family_dict[24] + "(n=" + "15" + ")",
                                                   family_dict[4] + "(n=" + "40" + ")",
                                                   family_dict[21] + "(n=" + "613" + ")",
                                                   family_dict[2] + "(n=" + "136" + ")",
                                                   family_dict[11] + "(n=" + "123" + ")",
                                                   family_dict[13] + "(n=" + "4" + ")",
                                                   family_dict[14] + "(n=" + "5" + ")",
                                                   family_dict[12] + "(n=" + "15" + ")",
                                                   family_dict[19] + "(n=" + "148" + ")",
                                                   family_dict[9] + "(n=" + "5" + ")",
                                                   family_dict[1] + "(n=" + "4" + ")",
                                                   family_dict[23] + "(n=" + "13" + ")",
                                                   family_dict[17] + "(n=" + "4" + ")"],
                                                   rotation=56,
                                                   ha="right",
                                                   rotation_mode="anchor",
                                                   fontsize=12)

    plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    plt.gcf().subplots_adjust(bottom=0.31, left=0.06, right=0.94, top=0.91)
    # plt.set_facecolor('xkcd:salmon')
    plt.legend()
    plt.show()


    # fold_mfe_ppv_16SrRNA = float(fold_mfe_ppv.iloc[1])
    # fold_bpp_ppv_16SrRNA = float(fold_bpp_ppv.iloc[1])
    # drtr_mfe_ppv_16SrRNA = float(drtr_mfe_ppv.iloc[1])
    # drtr_prob_ppv_16SrRNA = float(drtr_prob_ppv.iloc[1])
    # drtr_bpp_ppv_16SrRNA = float(drtr_bpp_ppv.iloc[1])
    #
    #
    # fold_mfe_tpr_16SrRNA = float(fold_mfe_tpr.iloc[1])
    # fold_bpp_tpr_16SrRNA = float(fold_bpp_tpr.iloc[1])
    # drtr_mfe_tpr_16SrRNA = float(drtr_mfe_tpr.iloc[1])
    # drtr_prob_tpr_16SrRNA = float(drtr_prob_tpr.iloc[1])
    # drtr_bpp_tpr_16SrRNA = float(drtr_bpp_tpr.iloc[1])
    #
    #
    # fold_mfe_f1_16SrRNA = float(fold_mfe_f1.iloc[1])
    # fold_bpp_f1_16SrRNA = float(fold_bpp_f1.iloc[1])
    # drtr_mfe_f1_16SrRNA = float(drtr_mfe_f1.iloc[1])
    # drtr_prob_f1_16SrRNA = float(drtr_prob_f1.iloc[1])
    # drtr_bpp_f1_16SrRNA = float(drtr_bpp_f1.iloc[1])
    #
    #
    # fold_mfe_mcc_16SrRNA = float(fold_mfe_mcc.iloc[1])
    # fold_bpp_mcc_16SrRNA = float(fold_bpp_mcc.iloc[1])
    # drtr_mfe_mcc_16SrRNA = float(drtr_mfe_mcc.iloc[1])
    # drtr_prob_mcc_16SrRNA = float(drtr_prob_mcc.iloc[1])
    # drtr_bpp_mcc_16SrRNA = float(drtr_bpp_mcc.iloc[1])
