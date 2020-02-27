#!python
# content: script that calulates tpr/ppv/f-measure and mcc from raw results;
# calculates values based on RNA families

### ================= import modules ======================================= ###
import pandas as pd
import math
import sys
import matplotlib.pyplot as plt
import argparse
from matplotlib.pyplot import figure

### ================= settings for pandas and numpy ======================== ###
# prevents output from bein truncated when printed
pd.set_option('display.max_rows', 10000) #some default pandas settings will chop data
pd.set_option('display.max_columns', 10000)
pd.set_option('display.width', 10000)
np.set_printoptions(threshold=sys.maxsize) #if printed out, matrix won't be truncated


### ================= argparse arguments =================================== ###
parser = argparse.ArgumentParser(description='Compare ct data with RNAfold')

parser.add_argument('--printout','-p', action='store_true',
                    default=False,
                    help='print all values')

parser.add_argument('--output','-o', action='store_true',
                    default=False,
                    help='produce csv file with arithmetic means')

parser.add_argument('--test', '-t', action='store_true',
                    default=False,
                    help='run with test csv')

parser.add_argument('--verbose', '-v', action='store_true',
                    default=False,
                    help='show raw data')


parser.add_argument('--label','-l', action='store',
                    default="xyz",
                    help='parse label to graph')

parser.add_argument('--input','-i', action='store',
                    default="comp_n_n_fin.csv",
                    help='specify cdv input file')

parser.add_argument('--graphics','-g', action='store_true',
                    default=False,
                    help='print vector graphics plot, default:OFF')

args=parser.parse_args()


if args.test == False:
    input_name = args.input

if args.test == True:
    input_name = "comp_mu_cp_fin.csv"

output_name = input_name[:-4] + "_f_mean.csv"


### ================= functions ============================================ ###
def average(column):
    n=len(column)
    sum = 0
    avg = 0
    for i in range(0, n):
        sum += column.iloc[i]
    avg = sum/n
    return avg



def tpr(tp_data,fn_data):
    n = len(tp_data)
    value = 0
    result = 0
    tp_sum = 0
    fn_sum = 0
    for i in range(0,n):
        tp = float(tp_data.iloc[i])
        fn = float(fn_data.iloc[i])


        tp_sum += tp
        fn_sum += fn
    result = tp_sum / (tp_sum + fn_sum)
    return result

def ppv(tp_data,fp_data):
    n = len(tp_data)
    result = 0
    value = 0
    tp_sum = 0
    fp_sum = 0
    for i in range(0,n):
        tp = float(tp_data.iloc[i])
        fp = float(fp_data.iloc[i])
        tp_sum += tp
        fp_sum += fp
    result = tp_sum / (tp_sum + fp_sum)
    return result


def f_measure(tp_data, fp_data, fn_data):
    ppv = 0
    tpr = 0
    tp = 0
    fp = 0
    fn = 0
    tp_sum = 0
    fp_sum = 0
    fn_sum = 0
    f = 0
    result = 0
    n = len(tp_data)
    # print(n)
    for i in range (0, n):
        tp = float(tp_data.iloc[i])
        fp = float(fp_data.iloc[i])
        fn = float(fn_data.iloc[i])
        tp_sum += tp
        fp_sum += fp
        fn_sum += fn
    ppv = tp_sum / (tp_sum + fp_sum)
    tpr = tp_sum / (tp_sum + fn_sum)
    f = 2 * (ppv * tpr)/(ppv + tpr)
    result = f
    return result


def mcc(tp_data,fp_data,tn_data,fn_data):
    n = len(tp_data)
    value = 0
    mcc = 0
    result = 0
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    tp_sum = 0
    fp_sum = 0
    tn_sum = 0
    fn_sum = 0
    for i in range(0, n):
        tp = float(tp_data.iloc[i])
        fp = float(fp_data.iloc[i])
        tn = float(tn_data.iloc[i])
        fn = float(fn_data.iloc[i])
        tp_sum += tp
        fp_sum += fp
        tn_sum += tn
        fn_sum += fn
    root = math.sqrt(((tp_sum + fp_sum) * (tp_sum + fn_sum) * (tn_sum + fp_sum) * (tn_sum + fn_sum)))
    mcc = (tp_sum * tn_sum - fp_sum * fn_sum) / root
    result = mcc
    return result


### =============================== file-import ============================ ###
raw_data = pd.read_csv( open(input_name),float_precision = "high")
data_length = len(raw_data)

seql = raw_data['sequencelength']
pps_fold = raw_data['RNAfold_pps']
pps_drtr = raw_data['drtr_pps']
drtr_tp_mfe = raw_data['drtr_tp_mfe']
drtr_fp_mfe = raw_data['drtr_fp_mfe']
drtr_tn_mfe = raw_data['drtr_tn_mfe']
drtr_fn_mfe = raw_data['drtr_fn_mfe']
drtr_tp = raw_data['drtr_tp']
drtr_fp = raw_data['drtr_fp']
drtr_tn = raw_data['drtr_tn']
drtr_fn = raw_data['drtr_fn']
drtr_tp_prob = raw_data['drtr_tp_prob']
drtr_fp_prob = raw_data['drtr_fp_prob']
drtr_tn_prob = raw_data['drtr_tn_prob']
drtr_fn_prob = raw_data['drtr_fn_prob']
RNAfold_mfe_tp = raw_data['RNAfold_mfe_tp']
RNAfold_mfe_fp = raw_data['RNAfold_mfe_fp']
RNAfold_mfe_tn = raw_data['RNAfold_mfe_tn']
RNAfold_mfe_fn = raw_data['RNAfold_mfe_fn']
RNAfold_bpp_tp = raw_data['RNAfold_bpp_tp']
RNAfold_bpp_fp = raw_data['RNAfold_bpp_fp']
RNAfold_bpp_tn = raw_data['RNAfold_bpp_tn']
RNAfold_bpp_fn = raw_data['RNAfold_bpp_fn']
lendata = len(pps_fold)


### ========================== split based on RNA-family =================== ###
family_data = raw_data
family_data.sort_values(by='family', inplace=True, axis=0)
family_data.set_index(keys=['family'], drop=False, inplace=True)
family=family_data['family'].unique().tolist()


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

# ======================= data for all RNA-families ======================== ###

#16SrRNA
sixteen_SrRNA=family_data.loc[family_data.family=='16SrRNA']

tp_sixteen_SrRNA_fold_mfe = (sixteen_SrRNA[['RNAfold_mfe_tp']])
fp_sixteen_SrRNA_fold_mfe = (sixteen_SrRNA[['RNAfold_mfe_fp']])
tn_sixteen_SrRNA_fold_mfe = (sixteen_SrRNA[['RNAfold_mfe_tn']])
fn_sixteen_SrRNA_fold_mfe = (sixteen_SrRNA[['RNAfold_mfe_fn']])

tp_sixteen_SrRNA_fold_bpp = (sixteen_SrRNA[['RNAfold_bpp_tp']])
fp_sixteen_SrRNA_fold_bpp = (sixteen_SrRNA[['RNAfold_bpp_fp']])
tn_sixteen_SrRNA_fold_bpp = (sixteen_SrRNA[['RNAfold_bpp_tn']])
fn_sixteen_SrRNA_fold_bpp = (sixteen_SrRNA[['RNAfold_bpp_fn']])

tp_sixteen_SrRNA_drtr = (sixteen_SrRNA[['drtr_tp']])
fp_sixteen_SrRNA_drtr = (sixteen_SrRNA[['drtr_fp']])
tn_sixteen_SrRNA_drtr = (sixteen_SrRNA[['drtr_tn']])
fn_sixteen_SrRNA_drtr = (sixteen_SrRNA[['drtr_fn']])

tp_sixteen_SrRNA_drtr_mfe = (sixteen_SrRNA[['drtr_tp_mfe']])
fp_sixteen_SrRNA_drtr_mfe = (sixteen_SrRNA[['drtr_fp_mfe']])
tn_sixteen_SrRNA_drtr_mfe = (sixteen_SrRNA[['drtr_tn_mfe']])
fn_sixteen_SrRNA_drtr_mfe = (sixteen_SrRNA[['drtr_fn_mfe']])

tp_sixteen_SrRNA_drtr_prob = (sixteen_SrRNA[['drtr_tp_prob']])
fp_sixteen_SrRNA_drtr_prob = (sixteen_SrRNA[['drtr_fp_prob']])
tn_sixteen_SrRNA_drtr_prob = (sixteen_SrRNA[['drtr_tn_prob']])
fn_sixteen_SrRNA_drtr_prob = (sixteen_SrRNA[['drtr_fn_prob']])

#TPR
tpr_sixteen_SrRNA_fold_mfe = tpr(tp_sixteen_SrRNA_fold_mfe, fn_sixteen_SrRNA_fold_mfe)
#print(tpr_sixteen_SrRNA_fold_mfe)

tpr_sixteen_SrRNA_fold_bpp = tpr(tp_sixteen_SrRNA_fold_bpp, fn_sixteen_SrRNA_fold_bpp)
#print(tpr_sixteen_SrRNA_fold_bpp)

tpr_sixteen_SrRNA_drtr = tpr(tp_sixteen_SrRNA_drtr, fn_sixteen_SrRNA_drtr)
#print(tpr_sixteen_SrRNA_drtr)

tpr_sixteen_SrRNA_drtr_mfe = tpr(tp_sixteen_SrRNA_drtr_mfe, fn_sixteen_SrRNA_drtr_mfe)

tpr_sixteen_SrRNA_drtr_prob = tpr(tp_sixteen_SrRNA_drtr_prob, fn_sixteen_SrRNA_drtr_prob)


#PPV
ppv_sixteen_SrRNA_fold_mfe = ppv(tp_sixteen_SrRNA_fold_mfe, fp_sixteen_SrRNA_fold_mfe)
#print(ppv_sixteen_SrRNA_fold_mfe)

ppv_sixteen_SrRNA_fold_bpp = ppv(tp_sixteen_SrRNA_fold_bpp, fp_sixteen_SrRNA_fold_bpp)
#print(ppv_sixteen_SrRNA_fold_bpp)

ppv_sixteen_SrRNA_drtr = ppv(tp_sixteen_SrRNA_drtr, fp_sixteen_SrRNA_drtr)
#print(ppv_sixteen_SrRNA_drtr)

ppv_sixteen_SrRNA_drtr_mfe = ppv(tp_sixteen_SrRNA_drtr_mfe, fp_sixteen_SrRNA_drtr_mfe)

ppv_sixteen_SrRNA_drtr_prob = ppv(tp_sixteen_SrRNA_drtr_prob, fp_sixteen_SrRNA_drtr_prob)


#F1
f1_sixteen_SrRNA_fold_mfe = f_measure(tp_sixteen_SrRNA_fold_mfe, fp_sixteen_SrRNA_fold_mfe, fn_sixteen_SrRNA_fold_mfe)
#print(f1_sixteen_SrRNA_fold_mfe)

f1_sixteen_SrRNA_fold_bpp = f_measure(tp_sixteen_SrRNA_fold_bpp, fp_sixteen_SrRNA_fold_bpp, fn_sixteen_SrRNA_fold_bpp)
#print(f1_sixteen_SrRNA_fold_bpp)

f1_sixteen_SrRNA_drtr = f_measure(tp_sixteen_SrRNA_drtr, fp_sixteen_SrRNA_drtr, fn_sixteen_SrRNA_drtr)
#print(f1_sixteen_SrRNA_drtr)

f1_sixteen_SrRNA_drtr_mfe = f_measure(tp_sixteen_SrRNA_drtr_mfe, fp_sixteen_SrRNA_drtr_mfe, fn_sixteen_SrRNA_drtr_mfe)

f1_sixteen_SrRNA_drtr_prob = f_measure(tp_sixteen_SrRNA_drtr_prob, fp_sixteen_SrRNA_drtr_prob, fn_sixteen_SrRNA_drtr_prob)


#MCC
mcc_sixteen_SrRNA_fold_mfe = mcc (tp_sixteen_SrRNA_fold_mfe, fp_sixteen_SrRNA_fold_mfe, tn_sixteen_SrRNA_fold_mfe, fn_sixteen_SrRNA_fold_mfe)
#print(mcc_sixteen_SrRNA_fold_mfe)

mcc_sixteen_SrRNA_fold_bpp = mcc (tp_sixteen_SrRNA_fold_bpp, fp_sixteen_SrRNA_fold_bpp, tn_sixteen_SrRNA_fold_bpp, fn_sixteen_SrRNA_fold_bpp)
#print(mcc_sixteen_SrRNA_fold_bpp)

mcc_sixteen_SrRNA_drtr = mcc (tp_sixteen_SrRNA_drtr, fp_sixteen_SrRNA_drtr, tn_sixteen_SrRNA_drtr, fn_sixteen_SrRNA_drtr)
#print(mcc_sixteen_SrRNA_drtr)

mcc_sixteen_SrRNA_drtr_mfe = mcc (tp_sixteen_SrRNA_drtr_mfe, fp_sixteen_SrRNA_drtr_mfe, tn_sixteen_SrRNA_drtr_mfe, fn_sixteen_SrRNA_drtr_mfe)

mcc_sixteen_SrRNA_drtr_prob = mcc (tp_sixteen_SrRNA_drtr_prob, fp_sixteen_SrRNA_drtr_prob, tn_sixteen_SrRNA_drtr_prob, fn_sixteen_SrRNA_drtr_prob)



#23SrRNA
twentythree_SrRNA=family_data.loc[family_data.family=='23SrRNA']

tp_twentythree_SrRNA_fold_mfe = (twentythree_SrRNA[['RNAfold_mfe_tp']])
fp_twentythree_SrRNA_fold_mfe = (twentythree_SrRNA[['RNAfold_mfe_fp']])
tn_twentythree_SrRNA_fold_mfe = (twentythree_SrRNA[['RNAfold_mfe_tn']])
fn_twentythree_SrRNA_fold_mfe = (twentythree_SrRNA[['RNAfold_mfe_fn']])

tp_twentythree_SrRNA_fold_bpp = (twentythree_SrRNA[['RNAfold_bpp_tp']])
fp_twentythree_SrRNA_fold_bpp = (twentythree_SrRNA[['RNAfold_bpp_fp']])
tn_twentythree_SrRNA_fold_bpp = (twentythree_SrRNA[['RNAfold_bpp_tn']])
fn_twentythree_SrRNA_fold_bpp = (twentythree_SrRNA[['RNAfold_bpp_fn']])

tp_twentythree_SrRNA_drtr = (twentythree_SrRNA[['drtr_tp']])
fp_twentythree_SrRNA_drtr = (twentythree_SrRNA[['drtr_fp']])
tn_twentythree_SrRNA_drtr = (twentythree_SrRNA[['drtr_tn']])
fn_twentythree_SrRNA_drtr = (twentythree_SrRNA[['drtr_fn']])

tp_twentythree_SrRNA_drtr_mfe = (twentythree_SrRNA[['drtr_tp_mfe']])
fp_twentythree_SrRNA_drtr_mfe = (twentythree_SrRNA[['drtr_fp_mfe']])
tn_twentythree_SrRNA_drtr_mfe = (twentythree_SrRNA[['drtr_tn_mfe']])
fn_twentythree_SrRNA_drtr_mfe = (twentythree_SrRNA[['drtr_fn_mfe']])

tp_twentythree_SrRNA_drtr_prob = (twentythree_SrRNA[['drtr_tp_prob']])
fp_twentythree_SrRNA_drtr_prob = (twentythree_SrRNA[['drtr_fp_prob']])
tn_twentythree_SrRNA_drtr_prob = (twentythree_SrRNA[['drtr_tn_prob']])
fn_twentythree_SrRNA_drtr_prob = (twentythree_SrRNA[['drtr_fn_prob']])


pps_twentythree_SrRNA_fold=(twentythree_SrRNA[['RNAfold_pps']])
avg_pps_twentythree_SrRNA_fold=average(pps_twentythree_SrRNA_fold)
#print(*avg_pps_twentythree_SrRNA_fold)

twentythree_SrRNA_drtr=family_data.loc[family_data.family=='23SrRNA']
pps_twentythree_SrRNA_drtr=(twentythree_SrRNA_drtr[['drtr_pps']])
avg_pps_twentythree_SrRNA_drtr=average(pps_twentythree_SrRNA_drtr)
#print(*avg_pps_twentythree_SrRNA_drtr)

#TPR
tpr_twentythree_SrRNA_fold_mfe = tpr(tp_twentythree_SrRNA_fold_mfe, fn_twentythree_SrRNA_fold_mfe)
#print(tpr_twentythree_SrRNA_fold_mfe)

tpr_twentythree_SrRNA_fold_bpp = tpr(tp_twentythree_SrRNA_fold_bpp, fn_twentythree_SrRNA_fold_bpp)
#print(tpr_twentythree_SrRNA_fold_bpp)

tpr_twentythree_SrRNA_drtr = tpr(tp_twentythree_SrRNA_drtr, fn_twentythree_SrRNA_drtr)
#print(tpr_twentythree_SrRNA_drtr)

tpr_twentythree_SrRNA_drtr_mfe = tpr(tp_twentythree_SrRNA_drtr_mfe, fn_twentythree_SrRNA_drtr_mfe)

tpr_twentythree_SrRNA_drtr_prob = tpr(tp_twentythree_SrRNA_drtr_prob, fn_twentythree_SrRNA_drtr_prob)


#PPV
ppv_twentythree_SrRNA_fold_mfe = ppv(tp_twentythree_SrRNA_fold_mfe, fp_twentythree_SrRNA_fold_mfe)
#print(ppv_twentythree_SrRNA_fold_mfe)

ppv_twentythree_SrRNA_fold_bpp = ppv(tp_twentythree_SrRNA_fold_bpp, fp_twentythree_SrRNA_fold_bpp)
#print(ppv_twentythree_SrRNA_fold_bpp)

ppv_twentythree_SrRNA_drtr = ppv(tp_twentythree_SrRNA_drtr, fp_twentythree_SrRNA_drtr)
#print(ppv_twentythree_SrRNA_drtr)

ppv_twentythree_SrRNA_drtr_mfe = ppv(tp_twentythree_SrRNA_drtr_mfe, fp_twentythree_SrRNA_drtr_mfe)

ppv_twentythree_SrRNA_drtr_prob = ppv(tp_twentythree_SrRNA_drtr_prob, fp_twentythree_SrRNA_drtr_prob)


#F1
f1_twentythree_SrRNA_fold_mfe = f_measure(tp_twentythree_SrRNA_fold_mfe, fp_twentythree_SrRNA_fold_mfe, fn_twentythree_SrRNA_fold_mfe)
#print(f1_twentythree_SrRNA_fold_mfe)

f1_twentythree_SrRNA_fold_bpp = f_measure(tp_twentythree_SrRNA_fold_bpp, fp_twentythree_SrRNA_fold_bpp, fn_twentythree_SrRNA_fold_bpp)
#print(f1_twentythree_SrRNA_fold_bpp)

f1_twentythree_SrRNA_drtr = f_measure(tp_twentythree_SrRNA_drtr, fp_twentythree_SrRNA_drtr, fn_twentythree_SrRNA_drtr)
#print(f1_twentythree_SrRNA_drtr)

f1_twentythree_SrRNA_drtr_mfe = f_measure(tp_twentythree_SrRNA_drtr_mfe, fp_twentythree_SrRNA_drtr_mfe, fn_twentythree_SrRNA_drtr_mfe)

f1_twentythree_SrRNA_drtr_prob = f_measure(tp_twentythree_SrRNA_drtr_prob, fp_twentythree_SrRNA_drtr_prob, fn_twentythree_SrRNA_drtr_prob)


#MCC
mcc_twentythree_SrRNA_fold_mfe = mcc (tp_twentythree_SrRNA_fold_mfe, fp_twentythree_SrRNA_fold_mfe, tn_twentythree_SrRNA_fold_mfe, fn_twentythree_SrRNA_fold_mfe)
#print(mcc_twentythree_SrRNA_fold_mfe)

mcc_twentythree_SrRNA_fold_bpp = mcc (tp_twentythree_SrRNA_fold_bpp, fp_twentythree_SrRNA_fold_bpp, tn_twentythree_SrRNA_fold_bpp, fn_twentythree_SrRNA_fold_bpp)
#print(mcc_twentythree_SrRNA_fold_bpp)

mcc_twentythree_SrRNA_drtr = mcc (tp_twentythree_SrRNA_drtr, fp_twentythree_SrRNA_drtr, tn_twentythree_SrRNA_drtr, fn_twentythree_SrRNA_drtr)
#print(mcc_twentythree_SrRNA_drtr)
mcc_twentythree_SrRNA_drtr_mfe = mcc (tp_twentythree_SrRNA_drtr_mfe, fp_twentythree_SrRNA_drtr_mfe, tn_twentythree_SrRNA_drtr_mfe, fn_twentythree_SrRNA_drtr_mfe)

mcc_twentythree_SrRNA_drtr_prob = mcc (tp_twentythree_SrRNA_drtr_prob, fp_twentythree_SrRNA_drtr_prob, tn_twentythree_SrRNA_drtr_prob, fn_twentythree_SrRNA_drtr_prob)



#5SrRNA
five_SrRNA=family_data.loc[family_data.family=='5SrRNA']

tp_five_SrRNA_fold_mfe = (five_SrRNA[['RNAfold_mfe_tp']])
fp_five_SrRNA_fold_mfe = (five_SrRNA[['RNAfold_mfe_fp']])
tn_five_SrRNA_fold_mfe = (five_SrRNA[['RNAfold_mfe_tn']])
fn_five_SrRNA_fold_mfe = (five_SrRNA[['RNAfold_mfe_fn']])

tp_five_SrRNA_fold_bpp = (five_SrRNA[['RNAfold_bpp_tp']])
fp_five_SrRNA_fold_bpp = (five_SrRNA[['RNAfold_bpp_fp']])
tn_five_SrRNA_fold_bpp = (five_SrRNA[['RNAfold_bpp_tn']])
fn_five_SrRNA_fold_bpp = (five_SrRNA[['RNAfold_bpp_fn']])

tp_five_SrRNA_drtr = (five_SrRNA[['drtr_tp']])
fp_five_SrRNA_drtr = (five_SrRNA[['drtr_fp']])
tn_five_SrRNA_drtr = (five_SrRNA[['drtr_tn']])
fn_five_SrRNA_drtr = (five_SrRNA[['drtr_fn']])

tp_five_SrRNA_drtr_mfe = (five_SrRNA[['drtr_tp_mfe']])
fp_five_SrRNA_drtr_mfe = (five_SrRNA[['drtr_fp_mfe']])
tn_five_SrRNA_drtr_mfe = (five_SrRNA[['drtr_tn_mfe']])
fn_five_SrRNA_drtr_mfe = (five_SrRNA[['drtr_fn_mfe']])

tp_five_SrRNA_drtr_prob = (five_SrRNA[['drtr_tp_prob']])
fp_five_SrRNA_drtr_prob = (five_SrRNA[['drtr_fp_prob']])
tn_five_SrRNA_drtr_prob = (five_SrRNA[['drtr_tn_prob']])
fn_five_SrRNA_drtr_prob = (five_SrRNA[['drtr_fn_prob']])

#TPR
tpr_five_SrRNA_fold_mfe = tpr(tp_five_SrRNA_fold_mfe, fn_five_SrRNA_fold_mfe)
#print(tpr_five_SrRNA_fold_mfe)

tpr_five_SrRNA_fold_bpp = tpr(tp_five_SrRNA_fold_bpp, fn_five_SrRNA_fold_bpp)
#print(tpr_five_SrRNA_fold_bpp)

tpr_five_SrRNA_drtr = tpr(tp_five_SrRNA_drtr, fn_five_SrRNA_drtr)
#print(tpr_five_SrRNA_drtr)

tpr_five_SrRNA_drtr_mfe = tpr(tp_five_SrRNA_drtr_mfe, fn_five_SrRNA_drtr_mfe)

tpr_five_SrRNA_drtr_prob = tpr(tp_five_SrRNA_drtr_prob, fn_five_SrRNA_drtr_prob)


#PPV
ppv_five_SrRNA_fold_mfe = ppv(tp_five_SrRNA_fold_mfe, fp_five_SrRNA_fold_mfe)
#print(ppv_five_SrRNA_fold_mfe)

ppv_five_SrRNA_fold_bpp = ppv(tp_five_SrRNA_fold_bpp, fp_five_SrRNA_fold_bpp)
#print(ppv_five_SrRNA_fold_bpp)

ppv_five_SrRNA_drtr = ppv(tp_five_SrRNA_drtr, fp_five_SrRNA_drtr)
#print(ppv_five_SrRNA_drtr)

ppv_five_SrRNA_drtr_mfe = ppv(tp_five_SrRNA_drtr_mfe, fp_five_SrRNA_drtr_mfe)

ppv_five_SrRNA_drtr_prob = ppv(tp_five_SrRNA_drtr_prob, fp_five_SrRNA_drtr_prob)


#F1
f1_five_SrRNA_fold_mfe = f_measure(tp_five_SrRNA_fold_mfe, fp_five_SrRNA_fold_mfe, fn_five_SrRNA_fold_mfe)
#print(f1_five_SrRNA_fold_mfe)

f1_five_SrRNA_fold_bpp = f_measure(tp_five_SrRNA_fold_bpp, fp_five_SrRNA_fold_bpp, fn_five_SrRNA_fold_bpp)
#print(f1_five_SrRNA_fold_bpp)

f1_five_SrRNA_drtr = f_measure(tp_five_SrRNA_drtr, fp_five_SrRNA_drtr, fn_five_SrRNA_drtr)
#print(f1_five_SrRNA_drtr)

f1_five_SrRNA_drtr_mfe = f_measure(tp_five_SrRNA_drtr_mfe, fp_five_SrRNA_drtr_mfe, fn_five_SrRNA_drtr_mfe)

f1_five_SrRNA_drtr_prob = f_measure(tp_five_SrRNA_drtr_prob, fp_five_SrRNA_drtr_prob, fn_five_SrRNA_drtr_prob)


#MCC
mcc_five_SrRNA_fold_mfe = mcc (tp_five_SrRNA_fold_mfe, fp_five_SrRNA_fold_mfe, tn_five_SrRNA_fold_mfe, fn_five_SrRNA_fold_mfe)
#print(mcc_five_SrRNA_fold_mfe)

mcc_five_SrRNA_fold_bpp = mcc (tp_five_SrRNA_fold_bpp, fp_five_SrRNA_fold_bpp, tn_five_SrRNA_fold_bpp, fn_five_SrRNA_fold_bpp)
#print(mcc_five_SrRNA_fold_bpp)

mcc_five_SrRNA_drtr = mcc (tp_five_SrRNA_drtr, fp_five_SrRNA_drtr, tn_five_SrRNA_drtr, fn_five_SrRNA_drtr)
#print(mcc_five_SrRNA_drtr)

mcc_five_SrRNA_drtr_mfe = mcc (tp_five_SrRNA_drtr_mfe, fp_five_SrRNA_drtr_mfe, tn_five_SrRNA_drtr_mfe, fn_five_SrRNA_drtr_mfe)

mcc_five_SrRNA_drtr_prob = mcc (tp_five_SrRNA_drtr_prob, fp_five_SrRNA_drtr_prob, tn_five_SrRNA_drtr_prob, fn_five_SrRNA_drtr_prob)



#Cili.Telo. RNA
cili_telo_RNA=family_data.loc[family_data.family=='Cili.Telo. RNA']

tp_cili_telo_RNA_fold_mfe = (cili_telo_RNA[['RNAfold_mfe_tp']])
fp_cili_telo_RNA_fold_mfe = (cili_telo_RNA[['RNAfold_mfe_fp']])
tn_cili_telo_RNA_fold_mfe = (cili_telo_RNA[['RNAfold_mfe_tn']])
fn_cili_telo_RNA_fold_mfe = (cili_telo_RNA[['RNAfold_mfe_fn']])

tp_cili_telo_RNA_fold_bpp = (cili_telo_RNA[['RNAfold_bpp_tp']])
fp_cili_telo_RNA_fold_bpp = (cili_telo_RNA[['RNAfold_bpp_fp']])
tn_cili_telo_RNA_fold_bpp = (cili_telo_RNA[['RNAfold_bpp_tn']])
fn_cili_telo_RNA_fold_bpp = (cili_telo_RNA[['RNAfold_bpp_fn']])

tp_cili_telo_RNA_drtr = (cili_telo_RNA[['drtr_tp']])
fp_cili_telo_RNA_drtr = (cili_telo_RNA[['drtr_fp']])
tn_cili_telo_RNA_drtr = (cili_telo_RNA[['drtr_tn']])
fn_cili_telo_RNA_drtr = (cili_telo_RNA[['drtr_fn']])

tp_cili_telo_RNA_drtr_mfe = (cili_telo_RNA[['drtr_tp_mfe']])
fp_cili_telo_RNA_drtr_mfe = (cili_telo_RNA[['drtr_fp_mfe']])
tn_cili_telo_RNA_drtr_mfe = (cili_telo_RNA[['drtr_tn_mfe']])
fn_cili_telo_RNA_drtr_mfe = (cili_telo_RNA[['drtr_fn_mfe']])

tp_cili_telo_RNA_drtr_prob = (cili_telo_RNA[['drtr_tp_prob']])
fp_cili_telo_RNA_drtr_prob = (cili_telo_RNA[['drtr_fp_prob']])
tn_cili_telo_RNA_drtr_prob = (cili_telo_RNA[['drtr_tn_prob']])
fn_cili_telo_RNA_drtr_prob = (cili_telo_RNA[['drtr_fn_prob']])


#TPR
tpr_cili_telo_RNA_fold_mfe = tpr(tp_cili_telo_RNA_fold_mfe, fn_cili_telo_RNA_fold_mfe)
#print(tpr_cili_telo_RNA_fold_mfe)

tpr_cili_telo_RNA_fold_bpp = tpr(tp_cili_telo_RNA_fold_bpp, fn_cili_telo_RNA_fold_bpp)
#print(tpr_cili_telo_RNA_fold_bpp)

tpr_cili_telo_RNA_drtr = tpr(tp_cili_telo_RNA_drtr, fn_cili_telo_RNA_drtr)
#print(tpr_cili_telo_RNA_drtr)

tpr_cili_telo_RNA_drtr_mfe = tpr(tp_cili_telo_RNA_drtr_mfe, fn_cili_telo_RNA_drtr_mfe)

tpr_cili_telo_RNA_drtr_prob = tpr(tp_cili_telo_RNA_drtr_prob, fn_cili_telo_RNA_drtr_prob)


#PPV
ppv_cili_telo_RNA_fold_mfe = ppv(tp_cili_telo_RNA_fold_mfe, fp_cili_telo_RNA_fold_mfe)
#print(ppv_cili_telo_RNA_fold_mfe)

ppv_cili_telo_RNA_fold_bpp = ppv(tp_cili_telo_RNA_fold_bpp, fp_cili_telo_RNA_fold_bpp)
#print(ppv_cili_telo_RNA_fold_bpp)

ppv_cili_telo_RNA_drtr = ppv(tp_cili_telo_RNA_drtr, fp_cili_telo_RNA_drtr)
#print(ppv_cili_telo_RNA_drtr)

ppv_cili_telo_RNA_drtr_mfe = ppv(tp_cili_telo_RNA_drtr_mfe, fp_cili_telo_RNA_drtr_mfe)

ppv_cili_telo_RNA_drtr_prob = ppv(tp_cili_telo_RNA_drtr_prob, fp_cili_telo_RNA_drtr_prob)


#F1
f1_cili_telo_RNA_fold_mfe = f_measure(tp_cili_telo_RNA_fold_mfe, fp_cili_telo_RNA_fold_mfe, fn_cili_telo_RNA_fold_mfe)
#print(f1_cili_telo_RNA_fold_mfe)

f1_cili_telo_RNA_fold_bpp = f_measure(tp_cili_telo_RNA_fold_bpp, fp_cili_telo_RNA_fold_bpp, fn_cili_telo_RNA_fold_bpp)
#print(f1_cili_telo_RNA_fold_bpp)

f1_cili_telo_RNA_drtr = f_measure(tp_cili_telo_RNA_drtr, fp_cili_telo_RNA_drtr, fn_cili_telo_RNA_drtr)
#print(f1_cili_telo_RNA_drtr)

f1_cili_telo_RNA_drtr_mfe = f_measure(tp_cili_telo_RNA_drtr_mfe, fp_cili_telo_RNA_drtr_mfe, fn_cili_telo_RNA_drtr_mfe)

f1_cili_telo_RNA_drtr_prob = f_measure(tp_cili_telo_RNA_drtr_prob, fp_cili_telo_RNA_drtr_prob, fn_cili_telo_RNA_drtr_prob)


#MCC
mcc_cili_telo_RNA_fold_mfe = mcc (tp_cili_telo_RNA_fold_mfe, fp_cili_telo_RNA_fold_mfe, tn_cili_telo_RNA_fold_mfe, fn_cili_telo_RNA_fold_mfe)
#print(mcc_cili_telo_RNA_fold_mfe)

mcc_cili_telo_RNA_fold_bpp = mcc (tp_cili_telo_RNA_fold_bpp, fp_cili_telo_RNA_fold_bpp, tn_cili_telo_RNA_fold_bpp, fn_cili_telo_RNA_fold_bpp)
#print(mcc_cili_telo_RNA_fold_bpp)

mcc_cili_telo_RNA_drtr = mcc (tp_cili_telo_RNA_drtr, fp_cili_telo_RNA_drtr, tn_cili_telo_RNA_drtr, fn_cili_telo_RNA_drtr)
#print(mcc_cili_telo_RNA_drtr)

mcc_cili_telo_RNA_drtr_mfe = mcc (tp_cili_telo_RNA_drtr_mfe, fp_cili_telo_RNA_drtr_mfe, tn_cili_telo_RNA_drtr_mfe, fn_cili_telo_RNA_drtr_mfe)

mcc_cili_telo_RNA_drtr_prob = mcc (tp_cili_telo_RNA_drtr_prob, fp_cili_telo_RNA_drtr_prob, tn_cili_telo_RNA_drtr_prob, fn_cili_telo_RNA_drtr_prob)


#Cis-reg.element
cis_regulatory_element=family_data.loc[family_data.family=='Cis-reg.element']

tp_cis_regulatory_element_fold_mfe = (cis_regulatory_element[['RNAfold_mfe_tp']])
fp_cis_regulatory_element_fold_mfe = (cis_regulatory_element[['RNAfold_mfe_fp']])
tn_cis_regulatory_element_fold_mfe = (cis_regulatory_element[['RNAfold_mfe_tn']])
fn_cis_regulatory_element_fold_mfe = (cis_regulatory_element[['RNAfold_mfe_fn']])

tp_cis_regulatory_element_fold_bpp = (cis_regulatory_element[['RNAfold_bpp_tp']])
fp_cis_regulatory_element_fold_bpp = (cis_regulatory_element[['RNAfold_bpp_fp']])
tn_cis_regulatory_element_fold_bpp = (cis_regulatory_element[['RNAfold_bpp_tn']])
fn_cis_regulatory_element_fold_bpp = (cis_regulatory_element[['RNAfold_bpp_fn']])

tp_cis_regulatory_element_drtr = (cis_regulatory_element[['drtr_tp']])
fp_cis_regulatory_element_drtr = (cis_regulatory_element[['drtr_fp']])
tn_cis_regulatory_element_drtr = (cis_regulatory_element[['drtr_tn']])
fn_cis_regulatory_element_drtr = (cis_regulatory_element[['drtr_fn']])

tp_cis_regulatory_element_drtr_mfe = (cis_regulatory_element[['drtr_tp_mfe']])
fp_cis_regulatory_element_drtr_mfe = (cis_regulatory_element[['drtr_fp_mfe']])
tn_cis_regulatory_element_drtr_mfe = (cis_regulatory_element[['drtr_tn_mfe']])
fn_cis_regulatory_element_drtr_mfe = (cis_regulatory_element[['drtr_fn_mfe']])

tp_cis_regulatory_element_drtr_prob = (cis_regulatory_element[['drtr_tp_prob']])
fp_cis_regulatory_element_drtr_prob = (cis_regulatory_element[['drtr_fp_prob']])
tn_cis_regulatory_element_drtr_prob = (cis_regulatory_element[['drtr_tn_prob']])
fn_cis_regulatory_element_drtr_prob = (cis_regulatory_element[['drtr_fn_prob']])

#TPR
tpr_cis_regulatory_element_fold_mfe = tpr(tp_cis_regulatory_element_fold_mfe, fn_cis_regulatory_element_fold_mfe)
#print(tpr_cis_regulatory_element_fold_mfe)

tpr_cis_regulatory_element_fold_bpp = tpr(tp_cis_regulatory_element_fold_bpp, fn_cis_regulatory_element_fold_bpp)
#print(tpr_cis_regulatory_element_fold_bpp)

tpr_cis_regulatory_element_drtr = tpr(tp_cis_regulatory_element_drtr, fn_cis_regulatory_element_drtr)
#print(tpr_cis_regulatory_element_drtr)

tpr_cis_regulatory_element_drtr_mfe = tpr(tp_cis_regulatory_element_drtr_mfe, fn_cis_regulatory_element_drtr_mfe)

tpr_cis_regulatory_element_drtr_prob = tpr(tp_cis_regulatory_element_drtr_prob, fn_cis_regulatory_element_drtr_prob)


#PPV
ppv_cis_regulatory_element_fold_mfe = ppv(tp_cis_regulatory_element_fold_mfe, fp_cis_regulatory_element_fold_mfe)
#print(ppv_cis_regulatory_element_fold_mfe)

ppv_cis_regulatory_element_fold_bpp = ppv(tp_cis_regulatory_element_fold_bpp, fp_cis_regulatory_element_fold_bpp)
#print(ppv_cis_regulatory_element_fold_bpp)

ppv_cis_regulatory_element_drtr = ppv(tp_cis_regulatory_element_drtr, fp_cis_regulatory_element_drtr)
#print(ppv_cis_regulatory_element_drtr)

ppv_cis_regulatory_element_drtr_mfe = ppv(tp_cis_regulatory_element_drtr_mfe, fp_cis_regulatory_element_drtr_mfe)

ppv_cis_regulatory_element_drtr_prob = ppv(tp_cis_regulatory_element_drtr_prob, fp_cis_regulatory_element_drtr_prob)


#F1
f1_cis_regulatory_element_fold_mfe = f_measure(tp_cis_regulatory_element_fold_mfe, fp_cis_regulatory_element_fold_mfe, fn_cis_regulatory_element_fold_mfe)
#print(f1_cis_regulatory_element_fold_mfe)

f1_cis_regulatory_element_fold_bpp = f_measure(tp_cis_regulatory_element_fold_bpp, fp_cis_regulatory_element_fold_bpp, fn_cis_regulatory_element_fold_bpp)
#print(f1_cis_regulatory_element_fold_bpp)

f1_cis_regulatory_element_drtr = f_measure(tp_cis_regulatory_element_drtr, fp_cis_regulatory_element_drtr, fn_cis_regulatory_element_drtr)
#print(f1_cis_regulatory_element_drtr)

f1_cis_regulatory_element_drtr_mfe = f_measure(tp_cis_regulatory_element_drtr_mfe, fp_cis_regulatory_element_drtr_mfe, fn_cis_regulatory_element_drtr_mfe)

f1_cis_regulatory_element_drtr_prob = f_measure(tp_cis_regulatory_element_drtr_prob, fp_cis_regulatory_element_drtr_prob, fn_cis_regulatory_element_drtr_prob)


#MCC
mcc_cis_regulatory_element_fold_mfe = mcc (tp_cis_regulatory_element_fold_mfe, fp_cis_regulatory_element_fold_mfe, tn_cis_regulatory_element_fold_mfe, fn_cis_regulatory_element_fold_mfe)
#print(mcc_cis_regulatory_element_fold_mfe)

mcc_cis_regulatory_element_fold_bpp = mcc (tp_cis_regulatory_element_fold_bpp, fp_cis_regulatory_element_fold_bpp, tn_cis_regulatory_element_fold_bpp, fn_cis_regulatory_element_fold_bpp)
#print(mcc_cis_regulatory_element_fold_bpp)

mcc_cis_regulatory_element_drtr = mcc (tp_cis_regulatory_element_drtr, fp_cis_regulatory_element_drtr, tn_cis_regulatory_element_drtr, fn_cis_regulatory_element_drtr)
#print(mcc_cis_regulatory_element_drtr)

mcc_cis_regulatory_element_drtr_mfe = mcc (tp_cis_regulatory_element_drtr_mfe, fp_cis_regulatory_element_drtr_mfe, tn_cis_regulatory_element_drtr_mfe, fn_cis_regulatory_element_drtr_mfe)

mcc_cis_regulatory_element_drtr_prob = mcc (tp_cis_regulatory_element_drtr_prob, fp_cis_regulatory_element_drtr_prob, tn_cis_regulatory_element_drtr_prob, fn_cis_regulatory_element_drtr_prob)


#GIIIntron
GIIIntron=family_data.loc[family_data.family=='GIIIntron']

tp_GIIIntron_fold_mfe = (GIIIntron[['RNAfold_mfe_tp']])
fp_GIIIntron_fold_mfe = (GIIIntron[['RNAfold_mfe_fp']])
tn_GIIIntron_fold_mfe = (GIIIntron[['RNAfold_mfe_tn']])
fn_GIIIntron_fold_mfe = (GIIIntron[['RNAfold_mfe_fn']])

tp_GIIIntron_fold_bpp = (GIIIntron[['RNAfold_bpp_tp']])
fp_GIIIntron_fold_bpp = (GIIIntron[['RNAfold_bpp_fp']])
tn_GIIIntron_fold_bpp = (GIIIntron[['RNAfold_bpp_tn']])
fn_GIIIntron_fold_bpp = (GIIIntron[['RNAfold_bpp_fn']])

tp_GIIIntron_drtr = (GIIIntron[['drtr_tp']])
fp_GIIIntron_drtr = (GIIIntron[['drtr_fp']])
tn_GIIIntron_drtr = (GIIIntron[['drtr_tn']])
fn_GIIIntron_drtr = (GIIIntron[['drtr_fn']])

tp_GIIIntron_drtr_mfe = (GIIIntron[['drtr_tp_mfe']])
fp_GIIIntron_drtr_mfe = (GIIIntron[['drtr_fp_mfe']])
tn_GIIIntron_drtr_mfe = (GIIIntron[['drtr_tn_mfe']])
fn_GIIIntron_drtr_mfe = (GIIIntron[['drtr_fn_mfe']])

tp_GIIIntron_drtr_prob = (GIIIntron[['drtr_tp_prob']])
fp_GIIIntron_drtr_prob = (GIIIntron[['drtr_fp_prob']])
tn_GIIIntron_drtr_prob = (GIIIntron[['drtr_tn_prob']])
fn_GIIIntron_drtr_prob = (GIIIntron[['drtr_fn_prob']])

#TPR
tpr_GIIIntron_fold_mfe = tpr(tp_GIIIntron_fold_mfe, fn_GIIIntron_fold_mfe)
#print(tpr_GIIIntron_fold_mfe)

tpr_GIIIntron_fold_bpp = tpr(tp_GIIIntron_fold_bpp, fn_GIIIntron_fold_bpp)
#print(tpr_GIIIntron_fold_bpp)

tpr_GIIIntron_drtr = tpr(tp_GIIIntron_drtr, fn_GIIIntron_drtr)
#print(tpr_GIIIntron_drtr)

tpr_GIIIntron_drtr_mfe = tpr(tp_GIIIntron_drtr_mfe, fn_GIIIntron_drtr_mfe)

tpr_GIIIntron_drtr_prob = tpr(tp_GIIIntron_drtr_prob, fn_GIIIntron_drtr_prob)


#PPV
ppv_GIIIntron_fold_mfe = ppv(tp_GIIIntron_fold_mfe, fp_GIIIntron_fold_mfe)
#print(ppv_GIIIntron_fold_mfe)

ppv_GIIIntron_fold_bpp = ppv(tp_GIIIntron_fold_bpp, fp_GIIIntron_fold_bpp)
#print(ppv_GIIIntron_fold_bpp)

ppv_GIIIntron_drtr = ppv(tp_GIIIntron_drtr, fp_GIIIntron_drtr)
#print(ppv_GIIIntron_drtr)

ppv_GIIIntron_drtr_mfe = ppv(tp_GIIIntron_drtr_mfe, fp_GIIIntron_drtr_mfe)

ppv_GIIIntron_drtr_prob = ppv(tp_GIIIntron_drtr_prob, fp_GIIIntron_drtr_prob)


#F1
f1_GIIIntron_fold_mfe = f_measure(tp_GIIIntron_fold_mfe, fp_GIIIntron_fold_mfe, fn_GIIIntron_fold_mfe)
#print(f1_GIIIntron_fold_mfe)

f1_GIIIntron_fold_bpp = f_measure(tp_GIIIntron_fold_bpp, fp_GIIIntron_fold_bpp, fn_GIIIntron_fold_bpp)
#print(f1_GIIIntron_fold_bpp)

f1_GIIIntron_drtr = f_measure(tp_GIIIntron_drtr, fp_GIIIntron_drtr, fn_GIIIntron_drtr)
#print(f1_GIIIntron_drtr)

f1_GIIIntron_drtr_mfe = f_measure(tp_GIIIntron_drtr_mfe, fp_GIIIntron_drtr_mfe, fn_GIIIntron_drtr_mfe)

f1_GIIIntron_drtr_prob = f_measure(tp_GIIIntron_drtr_prob, fp_GIIIntron_drtr_prob, fn_GIIIntron_drtr_prob)


#MCC
mcc_GIIIntron_fold_mfe = mcc (tp_GIIIntron_fold_mfe, fp_GIIIntron_fold_mfe, tn_GIIIntron_fold_mfe, fn_GIIIntron_fold_mfe)
#print(mcc_GIIIntron_fold_mfe)

mcc_GIIIntron_fold_bpp = mcc (tp_GIIIntron_fold_bpp, fp_GIIIntron_fold_bpp, tn_GIIIntron_fold_bpp, fn_GIIIntron_fold_bpp)
#print(mcc_GIIIntron_fold_bpp)

mcc_GIIIntron_drtr = mcc (tp_GIIIntron_drtr, fp_GIIIntron_drtr, tn_GIIIntron_drtr, fn_GIIIntron_drtr)
#print(mcc_GIIIntron_drtr)

mcc_GIIIntron_drtr_mfe = mcc (tp_GIIIntron_drtr_mfe, fp_GIIIntron_drtr_mfe, tn_GIIIntron_drtr_mfe, fn_GIIIntron_drtr_mfe)

mcc_GIIIntron_drtr_prob = mcc (tp_GIIIntron_drtr_prob, fp_GIIIntron_drtr_prob, tn_GIIIntron_drtr_prob, fn_GIIIntron_drtr_prob)



#GIIntron
GIIntron=family_data.loc[family_data.family=='GIIntron']

tp_GIIntron_fold_mfe = (GIIntron[['RNAfold_mfe_tp']])
fp_GIIntron_fold_mfe = (GIIntron[['RNAfold_mfe_fp']])
tn_GIIntron_fold_mfe = (GIIntron[['RNAfold_mfe_tn']])
fn_GIIntron_fold_mfe = (GIIntron[['RNAfold_mfe_fn']])

tp_GIIntron_fold_bpp = (GIIntron[['RNAfold_bpp_tp']])
fp_GIIntron_fold_bpp = (GIIntron[['RNAfold_bpp_fp']])
tn_GIIntron_fold_bpp = (GIIntron[['RNAfold_bpp_tn']])
fn_GIIntron_fold_bpp = (GIIntron[['RNAfold_bpp_fn']])

tp_GIIntron_drtr = (GIIntron[['drtr_tp']])
fp_GIIntron_drtr = (GIIntron[['drtr_fp']])
tn_GIIntron_drtr = (GIIntron[['drtr_tn']])
fn_GIIntron_drtr = (GIIntron[['drtr_fn']])

tp_GIIntron_drtr_mfe = (GIIntron[['drtr_tp_mfe']])
fp_GIIntron_drtr_mfe = (GIIntron[['drtr_fp_mfe']])
tn_GIIntron_drtr_mfe = (GIIntron[['drtr_tn_mfe']])
fn_GIIntron_drtr_mfe = (GIIntron[['drtr_fn_mfe']])

tp_GIIntron_drtr_prob = (GIIntron[['drtr_tp_prob']])
fp_GIIntron_drtr_prob = (GIIntron[['drtr_fp_prob']])
tn_GIIntron_drtr_prob = (GIIntron[['drtr_tn_prob']])
fn_GIIntron_drtr_prob = (GIIntron[['drtr_fn_prob']])


#TPR
tpr_GIIntron_fold_mfe = tpr(tp_GIIntron_fold_mfe, fn_GIIntron_fold_mfe)
#print(tpr_GIIntron_fold_mfe)

tpr_GIIntron_fold_bpp = tpr(tp_GIIntron_fold_bpp, fn_GIIntron_fold_bpp)
#print(tpr_GIIntron_fold_bpp)

tpr_GIIntron_drtr = tpr(tp_GIIntron_drtr, fn_GIIntron_drtr)
#print(tpr_GIIntron_drtr)

tpr_GIIntron_drtr_mfe = tpr(tp_GIIntron_drtr_mfe, fn_GIIntron_drtr_mfe)

tpr_GIIntron_drtr_prob = tpr(tp_GIIntron_drtr_prob, fn_GIIntron_drtr_prob)



#PPV
ppv_GIIntron_fold_mfe = ppv(tp_GIIntron_fold_mfe, fp_GIIntron_fold_mfe)
#print(ppv_GIIntron_fold_mfe)

ppv_GIIntron_fold_bpp = ppv(tp_GIIntron_fold_bpp, fp_GIIntron_fold_bpp)
#print(ppv_GIIntron_fold_bpp)

ppv_GIIntron_drtr = ppv(tp_GIIntron_drtr, fp_GIIntron_drtr)
#print(ppv_GIIntron_drtr)

ppv_GIIntron_drtr_mfe = ppv(tp_GIIntron_drtr_mfe, fp_GIIntron_drtr_mfe)

ppv_GIIntron_drtr_prob = ppv(tp_GIIntron_drtr_prob, fp_GIIntron_drtr_prob)



#F1
f1_GIIntron_fold_mfe = f_measure(tp_GIIntron_fold_mfe, fp_GIIntron_fold_mfe, fn_GIIntron_fold_mfe)
#print(f1_GIIntron_fold_mfe)

f1_GIIntron_fold_bpp = f_measure(tp_GIIntron_fold_bpp, fp_GIIntron_fold_bpp, fn_GIIntron_fold_bpp)
#print(f1_GIIntron_fold_bpp)

#print(f1_GIIntron_drtr)
f1_GIIntron_drtr = f_measure(tp_GIIntron_drtr, fp_GIIntron_drtr, fn_GIIntron_drtr)


f1_GIIntron_drtr_mfe = f_measure(tp_GIIntron_drtr_mfe, fp_GIIntron_drtr_mfe, fn_GIIntron_drtr_mfe)

f1_GIIntron_drtr_prob = f_measure(tp_GIIntron_drtr_prob, fp_GIIntron_drtr_prob, fn_GIIntron_drtr_prob)



#MCC
mcc_GIIntron_fold_mfe = mcc (tp_GIIntron_fold_mfe, fp_GIIntron_fold_mfe, tn_GIIntron_fold_mfe, fn_GIIntron_fold_mfe)
#print(mcc_GIIntron_fold_mfe)

mcc_GIIntron_fold_bpp = mcc (tp_GIIntron_fold_bpp, fp_GIIntron_fold_bpp, tn_GIIntron_fold_bpp, fn_GIIntron_fold_bpp)
#print(mcc_GIIntron_fold_bpp)

mcc_GIIntron_drtr = mcc (tp_GIIntron_drtr, fp_GIIntron_drtr, tn_GIIntron_drtr, fn_GIIntron_drtr)
#print(mcc_GIIntron_drtr)

mcc_GIIntron_drtr_mfe = mcc (tp_GIIntron_drtr_mfe, fp_GIIntron_drtr_mfe, tn_GIIntron_drtr_mfe, fn_GIIntron_drtr_mfe)

mcc_GIIntron_drtr_prob = mcc (tp_GIIntron_drtr_prob, fp_GIIntron_drtr_prob, tn_GIIntron_drtr_prob, fn_GIIntron_drtr_prob)



#Ham.Ribozyme
ham_ribozyme = family_data.loc[family_data.family=='Ham.Ribozyme']

tp_ham_ribozyme_fold_mfe = (ham_ribozyme[['RNAfold_mfe_tp']])
fp_ham_ribozyme_fold_mfe = (ham_ribozyme[['RNAfold_mfe_fp']])
tn_ham_ribozyme_fold_mfe = (ham_ribozyme[['RNAfold_mfe_tn']])
fn_ham_ribozyme_fold_mfe = (ham_ribozyme[['RNAfold_mfe_fn']])

tp_ham_ribozyme_fold_bpp = (ham_ribozyme[['RNAfold_bpp_tp']])
fp_ham_ribozyme_fold_bpp = (ham_ribozyme[['RNAfold_bpp_fp']])
tn_ham_ribozyme_fold_bpp = (ham_ribozyme[['RNAfold_bpp_tn']])
fn_ham_ribozyme_fold_bpp = (ham_ribozyme[['RNAfold_bpp_fn']])

tp_ham_ribozyme_drtr = (ham_ribozyme[['drtr_tp']])
fp_ham_ribozyme_drtr = (ham_ribozyme[['drtr_fp']])
tn_ham_ribozyme_drtr = (ham_ribozyme[['drtr_tn']])
fn_ham_ribozyme_drtr = (ham_ribozyme[['drtr_fn']])

tp_ham_ribozyme_drtr_mfe = (ham_ribozyme[['drtr_tp_mfe']])
fp_ham_ribozyme_drtr_mfe = (ham_ribozyme[['drtr_fp_mfe']])
tn_ham_ribozyme_drtr_mfe = (ham_ribozyme[['drtr_tn_mfe']])
fn_ham_ribozyme_drtr_mfe = (ham_ribozyme[['drtr_fn_mfe']])

tp_ham_ribozyme_drtr_prob = (ham_ribozyme[['drtr_tp_prob']])
fp_ham_ribozyme_drtr_prob = (ham_ribozyme[['drtr_fp_prob']])
tn_ham_ribozyme_drtr_prob = (ham_ribozyme[['drtr_tn_prob']])
fn_ham_ribozyme_drtr_prob = (ham_ribozyme[['drtr_fn_prob']])

#TPR
tpr_ham_ribozyme_fold_mfe = tpr(tp_ham_ribozyme_fold_mfe, fn_ham_ribozyme_fold_mfe)
#print(tpr_ham_ribozyme_fold_mfe)

tpr_ham_ribozyme_fold_bpp = tpr(tp_ham_ribozyme_fold_bpp, fn_ham_ribozyme_fold_bpp)
#print(tpr_ham_ribozyme_fold_bpp)

tpr_ham_ribozyme_drtr = tpr(tp_ham_ribozyme_drtr, fn_ham_ribozyme_drtr)
#print(tpr_ham_ribozyme_drtr)

tpr_ham_ribozyme_drtr_mfe = tpr(tp_ham_ribozyme_drtr_mfe, fn_ham_ribozyme_drtr_mfe)

tpr_ham_ribozyme_drtr_prob = tpr(tp_ham_ribozyme_drtr_prob, fn_ham_ribozyme_drtr_prob)



#PPV
ppv_ham_ribozyme_fold_mfe = ppv(tp_ham_ribozyme_fold_mfe, fp_ham_ribozyme_fold_mfe)
#print(ppv_ham_ribozyme_fold_mfe)

ppv_ham_ribozyme_fold_bpp = ppv(tp_ham_ribozyme_fold_bpp, fp_ham_ribozyme_fold_bpp)
#print(ppv_ham_ribozyme_fold_bpp)

ppv_ham_ribozyme_drtr = ppv(tp_ham_ribozyme_drtr, fp_ham_ribozyme_drtr)
#print(ppv_ham_ribozyme_drtr)

ppv_ham_ribozyme_drtr_mfe = ppv(tp_ham_ribozyme_drtr_mfe, fp_ham_ribozyme_drtr_mfe)

ppv_ham_ribozyme_drtr_prob = ppv(tp_ham_ribozyme_drtr_prob, fp_ham_ribozyme_drtr_prob)



#F1
f1_ham_ribozyme_fold_mfe = f_measure(tp_ham_ribozyme_fold_mfe, fp_ham_ribozyme_fold_mfe, fn_ham_ribozyme_fold_mfe)
#print(f1_ham_ribozyme_fold_mfe)

f1_ham_ribozyme_fold_bpp = f_measure(tp_ham_ribozyme_fold_bpp, fp_ham_ribozyme_fold_bpp, fn_ham_ribozyme_fold_bpp)
#print(f1_ham_ribozyme_fold_bpp)

f1_ham_ribozyme_drtr = f_measure(tp_ham_ribozyme_drtr, fp_ham_ribozyme_drtr, fn_ham_ribozyme_drtr)
#print(f1_ham_ribozyme_drtr)

f1_ham_ribozyme_drtr_mfe = f_measure(tp_ham_ribozyme_drtr_mfe, fp_ham_ribozyme_drtr_mfe, fn_ham_ribozyme_drtr_mfe)

f1_ham_ribozyme_drtr_prob = f_measure(tp_ham_ribozyme_drtr_prob, fp_ham_ribozyme_drtr_prob, fn_ham_ribozyme_drtr_prob)



#MCC
mcc_ham_ribozyme_fold_mfe = mcc (tp_ham_ribozyme_fold_mfe, fp_ham_ribozyme_fold_mfe, tn_ham_ribozyme_fold_mfe, fn_ham_ribozyme_fold_mfe)
#print(mcc_ham_ribozyme_fold_mfe)

mcc_ham_ribozyme_fold_bpp = mcc (tp_ham_ribozyme_fold_bpp, fp_ham_ribozyme_fold_bpp, tn_ham_ribozyme_fold_bpp, fn_ham_ribozyme_fold_bpp)
#print(mcc_ham_ribozyme_fold_bpp)

mcc_ham_ribozyme_drtr = mcc (tp_ham_ribozyme_drtr, fp_ham_ribozyme_drtr, tn_ham_ribozyme_drtr, fn_ham_ribozyme_drtr)
#print(mcc_ham_ribozyme_drtr)

mcc_ham_ribozyme_drtr_mfe = mcc (tp_ham_ribozyme_drtr_mfe, fp_ham_ribozyme_drtr_mfe, tn_ham_ribozyme_drtr_mfe, fn_ham_ribozyme_drtr_mfe)

mcc_ham_ribozyme_drtr_prob = mcc (tp_ham_ribozyme_drtr_prob, fp_ham_ribozyme_drtr_prob, tn_ham_ribozyme_drtr_prob, fn_ham_ribozyme_drtr_prob)



#HDVRibozyme
hdvr_ribozyme=family_data.loc[family_data.family=='HDVRibozyme']

tp_hdvr_ribozyme_fold_mfe = (hdvr_ribozyme[['RNAfold_mfe_tp']])
fp_hdvr_ribozyme_fold_mfe = (hdvr_ribozyme[['RNAfold_mfe_fp']])
tn_hdvr_ribozyme_fold_mfe = (hdvr_ribozyme[['RNAfold_mfe_tn']])
fn_hdvr_ribozyme_fold_mfe = (hdvr_ribozyme[['RNAfold_mfe_fn']])

tp_hdvr_ribozyme_fold_bpp = (hdvr_ribozyme[['RNAfold_bpp_tp']])
fp_hdvr_ribozyme_fold_bpp = (hdvr_ribozyme[['RNAfold_bpp_fp']])
tn_hdvr_ribozyme_fold_bpp = (hdvr_ribozyme[['RNAfold_bpp_tn']])
fn_hdvr_ribozyme_fold_bpp = (hdvr_ribozyme[['RNAfold_bpp_fn']])

tp_hdvr_ribozyme_drtr = (hdvr_ribozyme[['drtr_tp']])
fp_hdvr_ribozyme_drtr = (hdvr_ribozyme[['drtr_fp']])
tn_hdvr_ribozyme_drtr = (hdvr_ribozyme[['drtr_tn']])
fn_hdvr_ribozyme_drtr = (hdvr_ribozyme[['drtr_fn']])

tp_hdvr_ribozyme_drtr_mfe = (hdvr_ribozyme[['drtr_tp_mfe']])
fp_hdvr_ribozyme_drtr_mfe = (hdvr_ribozyme[['drtr_fp_mfe']])
tn_hdvr_ribozyme_drtr_mfe = (hdvr_ribozyme[['drtr_tn_mfe']])
fn_hdvr_ribozyme_drtr_mfe = (hdvr_ribozyme[['drtr_fn_mfe']])

tp_hdvr_ribozyme_drtr_prob = (hdvr_ribozyme[['drtr_tp_prob']])
fp_hdvr_ribozyme_drtr_prob = (hdvr_ribozyme[['drtr_fp_prob']])
tn_hdvr_ribozyme_drtr_prob = (hdvr_ribozyme[['drtr_tn_prob']])
fn_hdvr_ribozyme_drtr_prob = (hdvr_ribozyme[['drtr_fn_prob']])


#TPR
tpr_hdvr_ribozyme_fold_mfe = tpr(tp_hdvr_ribozyme_fold_mfe, fn_hdvr_ribozyme_fold_mfe)
#print(tpr_hdvr_ribozyme_fold_mfe)

tpr_hdvr_ribozyme_fold_bpp = tpr(tp_hdvr_ribozyme_fold_bpp, fn_hdvr_ribozyme_fold_bpp)
#print(tpr_hdvr_ribozyme_fold_bpp)

tpr_hdvr_ribozyme_drtr = tpr(tp_hdvr_ribozyme_drtr, fn_hdvr_ribozyme_drtr)
#print(tpr_hdvr_ribozyme_drtr)

tpr_hdvr_ribozyme_drtr_mfe = tpr(tp_hdvr_ribozyme_drtr_mfe, fn_hdvr_ribozyme_drtr_mfe)

tpr_hdvr_ribozyme_drtr_prob = tpr(tp_hdvr_ribozyme_drtr_prob, fn_hdvr_ribozyme_drtr_prob)


#PPV
ppv_hdvr_ribozyme_fold_mfe = ppv(tp_hdvr_ribozyme_fold_mfe, fp_hdvr_ribozyme_fold_mfe)
#print(ppv_hdvr_ribozyme_fold_mfe)

ppv_hdvr_ribozyme_fold_bpp = ppv(tp_hdvr_ribozyme_fold_bpp, fp_hdvr_ribozyme_fold_bpp)
#print(ppv_hdvr_ribozyme_fold_bpp)

ppv_hdvr_ribozyme_drtr = ppv(tp_hdvr_ribozyme_drtr, fp_hdvr_ribozyme_drtr)
#print(ppv_hdvr_ribozyme_drtr)

ppv_hdvr_ribozyme_drtr_mfe = ppv(tp_hdvr_ribozyme_drtr_mfe, fp_hdvr_ribozyme_drtr_mfe)

ppv_hdvr_ribozyme_drtr_prob = ppv(tp_hdvr_ribozyme_drtr_prob, fp_hdvr_ribozyme_drtr_prob)


#F1
f1_hdvr_ribozyme_fold_mfe = f_measure(tp_hdvr_ribozyme_fold_mfe, fp_hdvr_ribozyme_fold_mfe, fn_hdvr_ribozyme_fold_mfe)
#print(f1_hdvr_ribozyme_fold_mfe)

f1_hdvr_ribozyme_fold_bpp = f_measure(tp_hdvr_ribozyme_fold_bpp, fp_hdvr_ribozyme_fold_bpp, fn_hdvr_ribozyme_fold_bpp)
#print(f1_hdvr_ribozyme_fold_bpp)

f1_hdvr_ribozyme_drtr = f_measure(tp_hdvr_ribozyme_drtr, fp_hdvr_ribozyme_drtr, fn_hdvr_ribozyme_drtr)
#print(f1_hdvr_ribozyme_drtr)

f1_hdvr_ribozyme_drtr_mfe = f_measure(tp_hdvr_ribozyme_drtr_mfe, fp_hdvr_ribozyme_drtr_mfe, fn_hdvr_ribozyme_drtr_mfe)

f1_hdvr_ribozyme_drtr_prob = f_measure(tp_hdvr_ribozyme_drtr_prob, fp_hdvr_ribozyme_drtr_prob, fn_hdvr_ribozyme_drtr_prob)



#MCC
mcc_hdvr_ribozyme_fold_mfe = mcc (tp_hdvr_ribozyme_fold_mfe, fp_hdvr_ribozyme_fold_mfe, tn_hdvr_ribozyme_fold_mfe, fn_hdvr_ribozyme_fold_mfe)
#print(mcc_hdvr_ribozyme_fold_mfe)

mcc_hdvr_ribozyme_fold_bpp = mcc (tp_hdvr_ribozyme_fold_bpp, fp_hdvr_ribozyme_fold_bpp, tn_hdvr_ribozyme_fold_bpp, fn_hdvr_ribozyme_fold_bpp)
#print(mcc_hdvr_ribozyme_fold_bpp)

mcc_hdvr_ribozyme_drtr = mcc (tp_hdvr_ribozyme_drtr, fp_hdvr_ribozyme_drtr, tn_hdvr_ribozyme_drtr, fn_hdvr_ribozyme_drtr)
#print(mcc_hdvr_ribozyme_drtr)

mcc_hdvr_ribozyme_drtr_mfe = mcc (tp_hdvr_ribozyme_drtr_mfe, fp_hdvr_ribozyme_drtr_mfe, tn_hdvr_ribozyme_drtr_mfe, fn_hdvr_ribozyme_drtr_mfe)

mcc_hdvr_ribozyme_drtr_prob = mcc (tp_hdvr_ribozyme_drtr_prob, fp_hdvr_ribozyme_drtr_prob, tn_hdvr_ribozyme_drtr_prob, fn_hdvr_ribozyme_drtr_prob)



#IRES
ires=family_data.loc[family_data.family=='IRES']

tp_ires_fold_mfe = (ires[['RNAfold_mfe_tp']])
fp_ires_fold_mfe = (ires[['RNAfold_mfe_fp']])
tn_ires_fold_mfe = (ires[['RNAfold_mfe_tn']])
fn_ires_fold_mfe = (ires[['RNAfold_mfe_fn']])

tp_ires_fold_bpp = (ires[['RNAfold_bpp_tp']])
fp_ires_fold_bpp = (ires[['RNAfold_bpp_fp']])
tn_ires_fold_bpp = (ires[['RNAfold_bpp_tn']])
fn_ires_fold_bpp = (ires[['RNAfold_bpp_fn']])

tp_ires_drtr = (ires[['drtr_tp']])
fp_ires_drtr = (ires[['drtr_fp']])
tn_ires_drtr = (ires[['drtr_tn']])
fn_ires_drtr = (ires[['drtr_fn']])

tp_ires_drtr_mfe = (ires[['drtr_tp_mfe']])
fp_ires_drtr_mfe = (ires[['drtr_fp_mfe']])
tn_ires_drtr_mfe = (ires[['drtr_tn_mfe']])
fn_ires_drtr_mfe = (ires[['drtr_fn_mfe']])

tp_ires_drtr_prob = (ires[['drtr_tp_prob']])
fp_ires_drtr_prob = (ires[['drtr_fp_prob']])
tn_ires_drtr_prob = (ires[['drtr_tn_prob']])
fn_ires_drtr_prob = (ires[['drtr_fn_prob']])


#TPR
tpr_ires_fold_mfe = tpr(tp_ires_fold_mfe, fn_ires_fold_mfe)
#print(tpr_ires_fold_mfe)

tpr_ires_fold_bpp = tpr(tp_ires_fold_bpp, fn_ires_fold_bpp)
#print(tpr_ires_fold_bpp)

tpr_ires_drtr = tpr(tp_ires_drtr, fn_ires_drtr)
#print(tpr_ires_drtr)

tpr_ires_drtr_mfe = tpr(tp_ires_drtr_mfe, fn_ires_drtr_mfe)

tpr_ires_drtr_prob = tpr(tp_ires_drtr_prob, fn_ires_drtr_prob)



#PPV
ppv_ires_fold_mfe = ppv(tp_ires_fold_mfe, fp_ires_fold_mfe)
#print(ppv_ires_fold_mfe)

ppv_ires_fold_bpp = ppv(tp_ires_fold_bpp, fp_ires_fold_bpp)
#print(ppv_ires_fold_bpp)

ppv_ires_drtr = ppv(tp_ires_drtr, fp_ires_drtr)
#print(ppv_ires_drtr)

ppv_ires_drtr_mfe = ppv(tp_ires_drtr_mfe, fp_ires_drtr_mfe)

ppv_ires_drtr_prob = ppv(tp_ires_drtr_prob, fp_ires_drtr_prob)



#F1
f1_ires_fold_mfe = f_measure(tp_ires_fold_mfe, fp_ires_fold_mfe, fn_ires_fold_mfe)
#print(f1_ires_fold_mfe)

f1_ires_fold_bpp = f_measure(tp_ires_fold_bpp, fp_ires_fold_bpp, fn_ires_fold_bpp)
#print(f1_ires_fold_bpp)

f1_ires_drtr = f_measure(tp_ires_drtr, fp_ires_drtr, fn_ires_drtr)
#print(f1_ires_drtr)

f1_ires_drtr_mfe = f_measure(tp_ires_drtr_mfe, fp_ires_drtr_mfe, fn_ires_drtr_mfe)

f1_ires_drtr_prob = f_measure(tp_ires_drtr_prob, fp_ires_drtr_prob, fn_ires_drtr_prob)



#MCC
mcc_ires_fold_mfe = mcc (tp_ires_fold_mfe, fp_ires_fold_mfe, tn_ires_fold_mfe, fn_ires_fold_mfe)
#print(mcc_ires_fold_mfe)

mcc_ires_fold_bpp = mcc (tp_ires_fold_bpp, fp_ires_fold_bpp, tn_ires_fold_bpp, fn_ires_fold_bpp)
#print(mcc_ires_fold_bpp)

mcc_ires_drtr = mcc (tp_ires_drtr, fp_ires_drtr, tn_ires_drtr, fn_ires_drtr)
#print(mcc_ires_drtr)

mcc_ires_drtr_mfe = mcc (tp_ires_drtr_mfe, fp_ires_drtr_mfe, tn_ires_drtr_mfe, fn_ires_drtr_mfe)

mcc_ires_drtr_prob = mcc (tp_ires_drtr_prob, fp_ires_drtr_prob, tn_ires_drtr_prob, fn_ires_drtr_prob)



#OtherRibozyme
other_ribozyme=family_data.loc[family_data.family=='OtherRibozyme']

tp_other_ribozyme_fold_mfe = (other_ribozyme[['RNAfold_mfe_tp']])
fp_other_ribozyme_fold_mfe = (other_ribozyme[['RNAfold_mfe_fp']])
tn_other_ribozyme_fold_mfe = (other_ribozyme[['RNAfold_mfe_tn']])
fn_other_ribozyme_fold_mfe = (other_ribozyme[['RNAfold_mfe_fn']])

tp_other_ribozyme_fold_bpp = (other_ribozyme[['RNAfold_bpp_tp']])
fp_other_ribozyme_fold_bpp = (other_ribozyme[['RNAfold_bpp_fp']])
tn_other_ribozyme_fold_bpp = (other_ribozyme[['RNAfold_bpp_tn']])
fn_other_ribozyme_fold_bpp = (other_ribozyme[['RNAfold_bpp_fn']])

tp_other_ribozyme_drtr = (other_ribozyme[['drtr_tp']])
fp_other_ribozyme_drtr = (other_ribozyme[['drtr_fp']])
tn_other_ribozyme_drtr = (other_ribozyme[['drtr_tn']])
fn_other_ribozyme_drtr = (other_ribozyme[['drtr_fn']])

tp_other_ribozyme_drtr_mfe = (other_ribozyme[['drtr_tp_mfe']])
fp_other_ribozyme_drtr_mfe = (other_ribozyme[['drtr_fp_mfe']])
tn_other_ribozyme_drtr_mfe = (other_ribozyme[['drtr_tn_mfe']])
fn_other_ribozyme_drtr_mfe = (other_ribozyme[['drtr_fn_mfe']])

tp_other_ribozyme_drtr_prob = (other_ribozyme[['drtr_tp_prob']])
fp_other_ribozyme_drtr_prob = (other_ribozyme[['drtr_fp_prob']])
tn_other_ribozyme_drtr_prob = (other_ribozyme[['drtr_tn_prob']])
fn_other_ribozyme_drtr_prob = (other_ribozyme[['drtr_fn_prob']])



#TPR
tpr_other_ribozyme_fold_mfe = tpr(tp_other_ribozyme_fold_mfe, fn_other_ribozyme_fold_mfe)
#print(tpr_other_ribozyme_fold_mfe)

tpr_other_ribozyme_fold_bpp = tpr(tp_other_ribozyme_fold_bpp, fn_other_ribozyme_fold_bpp)
#print(tpr_other_ribozyme_fold_bpp)

tpr_other_ribozyme_drtr = tpr(tp_other_ribozyme_drtr, fn_other_ribozyme_drtr)
#print(tpr_other_ribozyme_drtr)

tpr_other_ribozyme_drtr_mfe = tpr(tp_other_ribozyme_drtr_mfe, fn_other_ribozyme_drtr_mfe)

tpr_other_ribozyme_drtr_prob = tpr(tp_other_ribozyme_drtr_prob, fn_other_ribozyme_drtr_prob)



#PPV
ppv_other_ribozyme_fold_mfe = ppv(tp_other_ribozyme_fold_mfe, fp_other_ribozyme_fold_mfe)
#print(ppv_other_ribozyme_fold_mfe)

ppv_other_ribozyme_fold_bpp = ppv(tp_other_ribozyme_fold_bpp, fp_other_ribozyme_fold_bpp)
#print(ppv_other_ribozyme_fold_bpp)

ppv_other_ribozyme_drtr = ppv(tp_other_ribozyme_drtr, fp_other_ribozyme_drtr)
#print(ppv_other_ribozyme_drtr)

ppv_other_ribozyme_drtr_mfe = ppv(tp_other_ribozyme_drtr_mfe, fp_other_ribozyme_drtr_mfe)

ppv_other_ribozyme_drtr_prob = ppv(tp_other_ribozyme_drtr_prob, fp_other_ribozyme_drtr_prob)



#F1
f1_other_ribozyme_fold_mfe = f_measure(tp_other_ribozyme_fold_mfe, fp_other_ribozyme_fold_mfe, fn_other_ribozyme_fold_mfe)
#print(f1_other_ribozyme_fold_mfe)

f1_other_ribozyme_fold_bpp = f_measure(tp_other_ribozyme_fold_bpp, fp_other_ribozyme_fold_bpp, fn_other_ribozyme_fold_bpp)
#print(f1_other_ribozyme_fold_bpp)

f1_other_ribozyme_drtr = f_measure(tp_other_ribozyme_drtr, fp_other_ribozyme_drtr, fn_other_ribozyme_drtr)
#print(f1_other_ribozyme_drtr)

f1_other_ribozyme_drtr_mfe = f_measure(tp_other_ribozyme_drtr_mfe, fp_other_ribozyme_drtr_mfe, fn_other_ribozyme_drtr_mfe)

f1_other_ribozyme_drtr_prob = f_measure(tp_other_ribozyme_drtr_prob, fp_other_ribozyme_drtr_prob, fn_other_ribozyme_drtr_prob)


#MCC
mcc_other_ribozyme_fold_mfe = mcc (tp_other_ribozyme_fold_mfe, fp_other_ribozyme_fold_mfe, tn_other_ribozyme_fold_mfe, fn_other_ribozyme_fold_mfe)
#print(mcc_other_ribozyme_fold_mfe)

mcc_other_ribozyme_fold_bpp = mcc (tp_other_ribozyme_fold_bpp, fp_other_ribozyme_fold_bpp, tn_other_ribozyme_fold_bpp, fn_other_ribozyme_fold_bpp)
#print(mcc_other_ribozyme_fold_bpp)

mcc_other_ribozyme_drtr = mcc (tp_other_ribozyme_drtr, fp_other_ribozyme_drtr, tn_other_ribozyme_drtr, fn_other_ribozyme_drtr)
#print(mcc_other_ribozyme_drtr)

mcc_other_ribozyme_drtr_mfe = mcc (tp_other_ribozyme_drtr_mfe, fp_other_ribozyme_drtr_mfe, tn_other_ribozyme_drtr_mfe, fn_other_ribozyme_drtr_mfe)

mcc_other_ribozyme_drtr_prob = mcc (tp_other_ribozyme_drtr_prob, fp_other_ribozyme_drtr_prob, tn_other_ribozyme_drtr_prob, fn_other_ribozyme_drtr_prob)


#OtherRNA
other_RNA=family_data.loc[family_data.family=='OtherRNA']

tp_other_RNA_fold_mfe = (other_RNA[['RNAfold_mfe_tp']])
fp_other_RNA_fold_mfe = (other_RNA[['RNAfold_mfe_fp']])
tn_other_RNA_fold_mfe = (other_RNA[['RNAfold_mfe_tn']])
fn_other_RNA_fold_mfe = (other_RNA[['RNAfold_mfe_fn']])

tp_other_RNA_fold_bpp = (other_RNA[['RNAfold_bpp_tp']])
fp_other_RNA_fold_bpp = (other_RNA[['RNAfold_bpp_fp']])
tn_other_RNA_fold_bpp = (other_RNA[['RNAfold_bpp_tn']])
fn_other_RNA_fold_bpp = (other_RNA[['RNAfold_bpp_fn']])

tp_other_RNA_drtr = (other_RNA[['drtr_tp']])
fp_other_RNA_drtr = (other_RNA[['drtr_fp']])
tn_other_RNA_drtr = (other_RNA[['drtr_tn']])
fn_other_RNA_drtr = (other_RNA[['drtr_fn']])

tp_other_RNA_drtr_mfe = (other_RNA[['drtr_tp_mfe']])
fp_other_RNA_drtr_mfe = (other_RNA[['drtr_fp_mfe']])
tn_other_RNA_drtr_mfe = (other_RNA[['drtr_tn_mfe']])
fn_other_RNA_drtr_mfe = (other_RNA[['drtr_fn_mfe']])

tp_other_RNA_drtr_prob = (other_RNA[['drtr_tp_prob']])
fp_other_RNA_drtr_prob = (other_RNA[['drtr_fp_prob']])
tn_other_RNA_drtr_prob = (other_RNA[['drtr_tn_prob']])
fn_other_RNA_drtr_prob = (other_RNA[['drtr_fn_prob']])


#TPR
tpr_other_RNA_fold_mfe = tpr(tp_other_RNA_fold_mfe, fn_other_RNA_fold_mfe)
#print(tpr_other_RNA_fold_mfe)

tpr_other_RNA_fold_bpp = tpr(tp_other_RNA_fold_bpp, fn_other_RNA_fold_bpp)
#print(tpr_other_RNA_fold_bpp)

tpr_other_RNA_drtr = tpr(tp_other_RNA_drtr, fn_other_RNA_drtr)
#print(tpr_other_RNA_drtr)

tpr_other_RNA_drtr_mfe = tpr(tp_other_RNA_drtr_mfe, fn_other_RNA_drtr_mfe)

tpr_other_RNA_drtr_prob = tpr(tp_other_RNA_drtr_prob, fn_other_RNA_drtr_prob)



#PPV
ppv_other_RNA_fold_mfe = ppv(tp_other_RNA_fold_mfe, fp_other_RNA_fold_mfe)
#print(ppv_other_RNA_fold_mfe)

ppv_other_RNA_fold_bpp = ppv(tp_other_RNA_fold_bpp, fp_other_RNA_fold_bpp)
#print(ppv_other_RNA_fold_bpp)

ppv_other_RNA_drtr = ppv(tp_other_RNA_drtr, fp_other_RNA_drtr)
#print(ppv_other_RNA_drtr)

ppv_other_RNA_drtr_mfe = ppv(tp_other_RNA_drtr_mfe, fp_other_RNA_drtr_mfe)

ppv_other_RNA_drtr_prob = ppv(tp_other_RNA_drtr_prob, fp_other_RNA_drtr_prob)



#F1
f1_other_RNA_fold_mfe = f_measure(tp_other_RNA_fold_mfe, fp_other_RNA_fold_mfe, fn_other_RNA_fold_mfe)
#print(f1_other_RNA_fold_mfe)

f1_other_RNA_fold_bpp = f_measure(tp_other_RNA_fold_bpp, fp_other_RNA_fold_bpp, fn_other_RNA_fold_bpp)
#print(f1_other_RNA_fold_bpp)

f1_other_RNA_drtr = f_measure(tp_other_RNA_drtr, fp_other_RNA_drtr, fn_other_RNA_drtr)
#print(f1_other_RNA_drtr)

f1_other_RNA_drtr_mfe = f_measure(tp_other_RNA_drtr_mfe, fp_other_RNA_drtr_mfe, fn_other_RNA_drtr_mfe)

f1_other_RNA_drtr_prob = f_measure(tp_other_RNA_drtr_prob, fp_other_RNA_drtr_prob, fn_other_RNA_drtr_prob)



#MCC
mcc_other_RNA_fold_mfe = mcc (tp_other_RNA_fold_mfe, fp_other_RNA_fold_mfe, tn_other_RNA_fold_mfe, fn_other_RNA_fold_mfe)
#print(mcc_other_RNA_fold_mfe)

mcc_other_RNA_fold_bpp = mcc (tp_other_RNA_fold_bpp, fp_other_RNA_fold_bpp, tn_other_RNA_fold_bpp, fn_other_RNA_fold_bpp)
#print(mcc_other_RNA_fold_bpp)

mcc_other_RNA_drtr = mcc (tp_other_RNA_drtr, fp_other_RNA_drtr, tn_other_RNA_drtr, fn_other_RNA_drtr)
#print(mcc_other_RNA_drtr)

mcc_other_RNA_drtr_mfe = mcc (tp_other_RNA_drtr_mfe, fp_other_RNA_drtr_mfe, tn_other_RNA_drtr_mfe, fn_other_RNA_drtr_mfe)

mcc_other_RNA_drtr_prob = mcc (tp_other_RNA_drtr_prob, fp_other_RNA_drtr_prob, tn_other_RNA_drtr_prob, fn_other_RNA_drtr_prob)


#OtherrRNA
other_rRNA=family_data.loc[family_data.family=='OtherrRNA']

tp_other_rRNA_fold_mfe = (other_rRNA[['RNAfold_mfe_tp']])
fp_other_rRNA_fold_mfe = (other_rRNA[['RNAfold_mfe_fp']])
tn_other_rRNA_fold_mfe = (other_rRNA[['RNAfold_mfe_tn']])
fn_other_rRNA_fold_mfe = (other_rRNA[['RNAfold_mfe_fn']])

tp_other_rRNA_fold_bpp = (other_rRNA[['RNAfold_bpp_tp']])
fp_other_rRNA_fold_bpp = (other_rRNA[['RNAfold_bpp_fp']])
tn_other_rRNA_fold_bpp = (other_rRNA[['RNAfold_bpp_tn']])
fn_other_rRNA_fold_bpp = (other_rRNA[['RNAfold_bpp_fn']])

tp_other_rRNA_drtr = (other_rRNA[['drtr_tp']])
fp_other_rRNA_drtr = (other_rRNA[['drtr_fp']])
tn_other_rRNA_drtr = (other_rRNA[['drtr_tn']])
fn_other_rRNA_drtr = (other_rRNA[['drtr_fn']])

tp_other_rRNA_drtr_mfe = (other_rRNA[['drtr_tp_mfe']])
fp_other_rRNA_drtr_mfe = (other_rRNA[['drtr_fp_mfe']])
tn_other_rRNA_drtr_mfe = (other_rRNA[['drtr_tn_mfe']])
fn_other_rRNA_drtr_mfe = (other_rRNA[['drtr_fn_mfe']])

tp_other_rRNA_drtr_prob = (other_rRNA[['drtr_tp_prob']])
fp_other_rRNA_drtr_prob = (other_rRNA[['drtr_fp_prob']])
tn_other_rRNA_drtr_prob = (other_rRNA[['drtr_tn_prob']])
fn_other_rRNA_drtr_prob = (other_rRNA[['drtr_fn_prob']])


#TPR
tpr_other_rRNA_fold_mfe = tpr(tp_other_rRNA_fold_mfe, fn_other_rRNA_fold_mfe)
#print(tpr_other_rRNA_fold_mfe)

tpr_other_rRNA_fold_bpp = tpr(tp_other_rRNA_fold_bpp, fn_other_rRNA_fold_bpp)
#print(tpr_other_rRNA_fold_bpp)

tpr_other_rRNA_drtr = tpr(tp_other_rRNA_drtr, fn_other_rRNA_drtr)
#print(tpr_other_rRNA_drtr)

tpr_other_rRNA_drtr_mfe = tpr(tp_other_rRNA_drtr_mfe, fn_other_rRNA_drtr_mfe)

tpr_other_rRNA_drtr_prob = tpr(tp_other_rRNA_drtr_prob, fn_other_rRNA_drtr_prob)



#PPV
ppv_other_rRNA_fold_mfe = ppv(tp_other_rRNA_fold_mfe, fp_other_rRNA_fold_mfe)
#print(ppv_other_rRNA_fold_mfe)

ppv_other_rRNA_fold_bpp = ppv(tp_other_rRNA_fold_bpp, fp_other_rRNA_fold_bpp)
#print(ppv_other_rRNA_fold_bpp)

ppv_other_rRNA_drtr = ppv(tp_other_rRNA_drtr, fp_other_rRNA_drtr)
#print(ppv_other_rRNA_drtr)

ppv_other_rRNA_drtr_mfe = ppv(tp_other_rRNA_drtr_mfe, fp_other_rRNA_drtr_mfe)

ppv_other_rRNA_drtr_prob = ppv(tp_other_rRNA_drtr_prob, fp_other_rRNA_drtr_prob)



#F1
f1_other_rRNA_fold_mfe = f_measure(tp_other_rRNA_fold_mfe, fp_other_rRNA_fold_mfe, fn_other_rRNA_fold_mfe)
#print(f1_other_rRNA_fold_mfe)

f1_other_rRNA_fold_bpp = f_measure(tp_other_rRNA_fold_bpp, fp_other_rRNA_fold_bpp, fn_other_rRNA_fold_bpp)
#print(f1_other_rRNA_fold_bpp)

f1_other_rRNA_drtr = f_measure(tp_other_rRNA_drtr, fp_other_rRNA_drtr, fn_other_rRNA_drtr)
#print(f1_other_rRNA_drtr)

f1_other_rRNA_drtr_mfe = f_measure(tp_other_rRNA_drtr_mfe, fp_other_rRNA_drtr_mfe, fn_other_rRNA_drtr_mfe)

f1_other_rRNA_drtr_prob = f_measure(tp_other_rRNA_drtr_prob, fp_other_rRNA_drtr_prob, fn_other_rRNA_drtr_prob)



#MCC
mcc_other_rRNA_fold_mfe = mcc (tp_other_rRNA_fold_mfe, fp_other_rRNA_fold_mfe, tn_other_rRNA_fold_mfe, fn_other_rRNA_fold_mfe)
#print(mcc_other_rRNA_fold_mfe)

mcc_other_rRNA_fold_bpp = mcc (tp_other_rRNA_fold_bpp, fp_other_rRNA_fold_bpp, tn_other_rRNA_fold_bpp, fn_other_rRNA_fold_bpp)
#print(mcc_other_rRNA_fold_bpp)

mcc_other_rRNA_drtr = mcc (tp_other_rRNA_drtr, fp_other_rRNA_drtr, tn_other_rRNA_drtr, fn_other_rRNA_drtr)
#print(mcc_other_rRNA_drtr)

mcc_other_rRNA_drtr_mfe = mcc (tp_other_rRNA_drtr_mfe, fp_other_rRNA_drtr_mfe, tn_other_rRNA_drtr_mfe, fn_other_rRNA_drtr_mfe)

mcc_other_rRNA_drtr_prob = mcc (tp_other_rRNA_drtr_prob, fp_other_rRNA_drtr_prob, tn_other_rRNA_drtr_prob, fn_other_rRNA_drtr_prob)



#RNAIII
RNAIII=family_data.loc[family_data.family=='RNAIII']

tp_RNAIII_fold_mfe = (RNAIII[['RNAfold_mfe_tp']])
fp_RNAIII_fold_mfe = (RNAIII[['RNAfold_mfe_fp']])
tn_RNAIII_fold_mfe = (RNAIII[['RNAfold_mfe_tn']])
fn_RNAIII_fold_mfe = (RNAIII[['RNAfold_mfe_fn']])

tp_RNAIII_fold_bpp = (RNAIII[['RNAfold_bpp_tp']])
fp_RNAIII_fold_bpp = (RNAIII[['RNAfold_bpp_fp']])
tn_RNAIII_fold_bpp = (RNAIII[['RNAfold_bpp_tn']])
fn_RNAIII_fold_bpp = (RNAIII[['RNAfold_bpp_fn']])

tp_RNAIII_drtr = (RNAIII[['drtr_tp']])
fp_RNAIII_drtr = (RNAIII[['drtr_fp']])
tn_RNAIII_drtr = (RNAIII[['drtr_tn']])
fn_RNAIII_drtr = (RNAIII[['drtr_fn']])

tp_RNAIII_drtr_mfe = (RNAIII[['drtr_tp_mfe']])
fp_RNAIII_drtr_mfe = (RNAIII[['drtr_fp_mfe']])
tn_RNAIII_drtr_mfe = (RNAIII[['drtr_tn_mfe']])
fn_RNAIII_drtr_mfe = (RNAIII[['drtr_fn_mfe']])

tp_RNAIII_drtr_prob = (RNAIII[['drtr_tp_prob']])
fp_RNAIII_drtr_prob = (RNAIII[['drtr_fp_prob']])
tn_RNAIII_drtr_prob = (RNAIII[['drtr_tn_prob']])
fn_RNAIII_drtr_prob = (RNAIII[['drtr_fn_prob']])


#TPR
tpr_RNAIII_fold_mfe = tpr(tp_RNAIII_fold_mfe, fn_RNAIII_fold_mfe)
#print(tpr_RNAIII_fold_mfe)

tpr_RNAIII_fold_bpp = tpr(tp_RNAIII_fold_bpp, fn_RNAIII_fold_bpp)
#print(tpr_RNAIII_fold_bpp)

tpr_RNAIII_drtr = tpr(tp_RNAIII_drtr, fn_RNAIII_drtr)
#print(tpr_RNAIII_drtr)

tpr_RNAIII_drtr_mfe = tpr(tp_RNAIII_drtr_mfe, fn_RNAIII_drtr_mfe)

tpr_RNAIII_drtr_prob = tpr(tp_RNAIII_drtr_prob, fn_RNAIII_drtr_prob)



#PPV
ppv_RNAIII_fold_mfe = ppv(tp_RNAIII_fold_mfe, fp_RNAIII_fold_mfe)
#print(ppv_RNAIII_fold_mfe)

ppv_RNAIII_fold_bpp = ppv(tp_RNAIII_fold_bpp, fp_RNAIII_fold_bpp)
#print(ppv_RNAIII_fold_bpp)

ppv_RNAIII_drtr = ppv(tp_RNAIII_drtr, fp_RNAIII_drtr)
#print(ppv_RNAIII_drtr)

ppv_RNAIII_drtr_mfe = ppv(tp_RNAIII_drtr_mfe, fp_RNAIII_drtr_mfe)

ppv_RNAIII_drtr_prob = ppv(tp_RNAIII_drtr_prob, fp_RNAIII_drtr_prob)



#F1
f1_RNAIII_fold_mfe = f_measure(tp_RNAIII_fold_mfe, fp_RNAIII_fold_mfe, fn_RNAIII_fold_mfe)
#print(f1_RNAIII_fold_mfe)

f1_RNAIII_fold_bpp = f_measure(tp_RNAIII_fold_bpp, fp_RNAIII_fold_bpp, fn_RNAIII_fold_bpp)
#print(f1_RNAIII_fold_bpp)

f1_RNAIII_drtr = f_measure(tp_RNAIII_drtr, fp_RNAIII_drtr, fn_RNAIII_drtr)
#print(f1_RNAIII_drtr)

f1_RNAIII_drtr_mfe = f_measure(tp_RNAIII_drtr_mfe, fp_RNAIII_drtr_mfe, fn_RNAIII_drtr_mfe)

f1_RNAIII_drtr_prob = f_measure(tp_RNAIII_drtr_prob, fp_RNAIII_drtr_prob, fn_RNAIII_drtr_prob)



#MCC
mcc_RNAIII_fold_mfe = mcc (tp_RNAIII_fold_mfe, fp_RNAIII_fold_mfe, tn_RNAIII_fold_mfe, fn_RNAIII_fold_mfe)
#print(mcc_RNAIII_fold_mfe)

mcc_RNAIII_fold_bpp = mcc (tp_RNAIII_fold_bpp, fp_RNAIII_fold_bpp, tn_RNAIII_fold_bpp, fn_RNAIII_fold_bpp)
#print(mcc_RNAIII_fold_bpp)

mcc_RNAIII_drtr = mcc (tp_RNAIII_drtr, fp_RNAIII_drtr, tn_RNAIII_drtr, fn_RNAIII_drtr)
#print(mcc_RNAIII_drtr)

mcc_RNAIII_drtr_mfe = mcc (tp_RNAIII_drtr_mfe, fp_RNAIII_drtr_mfe, tn_RNAIII_drtr_mfe, fn_RNAIII_drtr_mfe)

mcc_RNAIII_drtr_prob = mcc (tp_RNAIII_drtr_prob, fp_RNAIII_drtr_prob, tn_RNAIII_drtr_prob, fn_RNAIII_drtr_prob)



#RNaseE5UTR
RNaseE5UTR=family_data.loc[family_data.family=='RNaseE5UTR']

tp_RNaseE5UTR_fold_mfe = (RNaseE5UTR[['RNAfold_mfe_tp']])
fp_RNaseE5UTR_fold_mfe = (RNaseE5UTR[['RNAfold_mfe_fp']])
tn_RNaseE5UTR_fold_mfe = (RNaseE5UTR[['RNAfold_mfe_tn']])
fn_RNaseE5UTR_fold_mfe = (RNaseE5UTR[['RNAfold_mfe_fn']])

tp_RNaseE5UTR_fold_bpp = (RNaseE5UTR[['RNAfold_bpp_tp']])
fp_RNaseE5UTR_fold_bpp = (RNaseE5UTR[['RNAfold_bpp_fp']])
tn_RNaseE5UTR_fold_bpp = (RNaseE5UTR[['RNAfold_bpp_tn']])
fn_RNaseE5UTR_fold_bpp = (RNaseE5UTR[['RNAfold_bpp_fn']])

tp_RNaseE5UTR_drtr = (RNaseE5UTR[['drtr_tp']])
fp_RNaseE5UTR_drtr = (RNaseE5UTR[['drtr_fp']])
tn_RNaseE5UTR_drtr = (RNaseE5UTR[['drtr_tn']])
fn_RNaseE5UTR_drtr = (RNaseE5UTR[['drtr_fn']])

tp_RNaseE5UTR_drtr_mfe = (RNaseE5UTR[['drtr_tp_mfe']])
fp_RNaseE5UTR_drtr_mfe = (RNaseE5UTR[['drtr_fp_mfe']])
tn_RNaseE5UTR_drtr_mfe = (RNaseE5UTR[['drtr_tn_mfe']])
fn_RNaseE5UTR_drtr_mfe = (RNaseE5UTR[['drtr_fn_mfe']])

tp_RNaseE5UTR_drtr_prob = (RNaseE5UTR[['drtr_tp_prob']])
fp_RNaseE5UTR_drtr_prob = (RNaseE5UTR[['drtr_fp_prob']])
tn_RNaseE5UTR_drtr_prob = (RNaseE5UTR[['drtr_tn_prob']])
fn_RNaseE5UTR_drtr_prob = (RNaseE5UTR[['drtr_fn_prob']])


#TPR
tpr_RNaseE5UTR_fold_mfe = tpr(tp_RNaseE5UTR_fold_mfe, fn_RNaseE5UTR_fold_mfe)
#print(tpr_RNaseE5UTR_fold_mfe)

tpr_RNaseE5UTR_fold_bpp = tpr(tp_RNaseE5UTR_fold_bpp, fn_RNaseE5UTR_fold_bpp)
#print(tpr_RNaseE5UTR_fold_bpp)

tpr_RNaseE5UTR_drtr = tpr(tp_RNaseE5UTR_drtr, fn_RNaseE5UTR_drtr)
#print(tpr_RNaseE5UTR_drtr)

tpr_RNaseE5UTR_drtr_mfe = tpr(tp_RNaseE5UTR_drtr_mfe, fn_RNaseE5UTR_drtr_mfe)

tpr_RNaseE5UTR_drtr_prob = tpr(tp_RNaseE5UTR_drtr_prob, fn_RNaseE5UTR_drtr_prob)



#PPV
ppv_RNaseE5UTR_fold_mfe = ppv(tp_RNaseE5UTR_fold_mfe, fp_RNaseE5UTR_fold_mfe)
#print(ppv_RNaseE5UTR_fold_mfe)

ppv_RNaseE5UTR_fold_bpp = ppv(tp_RNaseE5UTR_fold_bpp, fp_RNaseE5UTR_fold_bpp)
#print(ppv_RNaseE5UTR_fold_bpp)

ppv_RNaseE5UTR_drtr = ppv(tp_RNaseE5UTR_drtr, fp_RNaseE5UTR_drtr)
#print(ppv_RNaseE5UTR_drtr)

ppv_RNaseE5UTR_drtr_mfe = ppv(tp_RNaseE5UTR_drtr_mfe, fp_RNaseE5UTR_drtr_mfe)

ppv_RNaseE5UTR_drtr_prob = ppv(tp_RNaseE5UTR_drtr_prob, fp_RNaseE5UTR_drtr_prob)



#F1
f1_RNaseE5UTR_fold_mfe = f_measure(tp_RNaseE5UTR_fold_mfe, fp_RNaseE5UTR_fold_mfe, fn_RNaseE5UTR_fold_mfe)
#print(f1_RNaseE5UTR_fold_mfe)

f1_RNaseE5UTR_fold_bpp = f_measure(tp_RNaseE5UTR_fold_bpp, fp_RNaseE5UTR_fold_bpp, fn_RNaseE5UTR_fold_bpp)
#print(f1_RNaseE5UTR_fold_bpp)

f1_RNaseE5UTR_drtr = f_measure(tp_RNaseE5UTR_drtr, fp_RNaseE5UTR_drtr, fn_RNaseE5UTR_drtr)
#print(f1_RNaseE5UTR_drtr)

f1_RNaseE5UTR_drtr_mfe = f_measure(tp_RNaseE5UTR_drtr_mfe, fp_RNaseE5UTR_drtr_mfe, fn_RNaseE5UTR_drtr_mfe)

f1_RNaseE5UTR_drtr_prob = f_measure(tp_RNaseE5UTR_drtr_prob, fp_RNaseE5UTR_drtr_prob, fn_RNaseE5UTR_drtr_prob)



#MCC
mcc_RNaseE5UTR_fold_mfe = mcc (tp_RNaseE5UTR_fold_mfe, fp_RNaseE5UTR_fold_mfe, tn_RNaseE5UTR_fold_mfe, fn_RNaseE5UTR_fold_mfe)
#print(mcc_RNaseE5UTR_fold_mfe)

mcc_RNaseE5UTR_fold_bpp = mcc (tp_RNaseE5UTR_fold_bpp, fp_RNaseE5UTR_fold_bpp, tn_RNaseE5UTR_fold_bpp, fn_RNaseE5UTR_fold_bpp)
#print(mcc_RNaseE5UTR_fold_bpp)

mcc_RNaseE5UTR_drtr = mcc (tp_RNaseE5UTR_drtr, fp_RNaseE5UTR_drtr, tn_RNaseE5UTR_drtr, fn_RNaseE5UTR_drtr)
#print(mcc_RNaseE5UTR_drtr)

mcc_RNaseE5UTR_drtr_mfe = mcc (tp_RNaseE5UTR_drtr_mfe, fp_RNaseE5UTR_drtr_mfe, tn_RNaseE5UTR_drtr_mfe, fn_RNaseE5UTR_drtr_mfe)

mcc_RNaseE5UTR_drtr_prob = mcc (tp_RNaseE5UTR_drtr_prob, fp_RNaseE5UTR_drtr_prob, tn_RNaseE5UTR_drtr_prob, fn_RNaseE5UTR_drtr_prob)



#RNaseMRPRNA
RNaseMRPRNA=family_data.loc[family_data.family=='RNaseMRPRNA']

tp_RNaseMRPRNA_fold_mfe = (RNaseMRPRNA[['RNAfold_mfe_tp']])
fp_RNaseMRPRNA_fold_mfe = (RNaseMRPRNA[['RNAfold_mfe_fp']])
tn_RNaseMRPRNA_fold_mfe = (RNaseMRPRNA[['RNAfold_mfe_tn']])
fn_RNaseMRPRNA_fold_mfe = (RNaseMRPRNA[['RNAfold_mfe_fn']])

tp_RNaseMRPRNA_fold_bpp = (RNaseMRPRNA[['RNAfold_bpp_tp']])
fp_RNaseMRPRNA_fold_bpp = (RNaseMRPRNA[['RNAfold_bpp_fp']])
tn_RNaseMRPRNA_fold_bpp = (RNaseMRPRNA[['RNAfold_bpp_tn']])
fn_RNaseMRPRNA_fold_bpp = (RNaseMRPRNA[['RNAfold_bpp_fn']])

tp_RNaseMRPRNA_drtr = (RNaseMRPRNA[['drtr_tp']])
fp_RNaseMRPRNA_drtr = (RNaseMRPRNA[['drtr_fp']])
tn_RNaseMRPRNA_drtr = (RNaseMRPRNA[['drtr_tn']])
fn_RNaseMRPRNA_drtr = (RNaseMRPRNA[['drtr_fn']])

tp_RNaseMRPRNA_drtr_mfe = (RNaseMRPRNA[['drtr_tp_mfe']])
fp_RNaseMRPRNA_drtr_mfe = (RNaseMRPRNA[['drtr_fp_mfe']])
tn_RNaseMRPRNA_drtr_mfe = (RNaseMRPRNA[['drtr_tn_mfe']])
fn_RNaseMRPRNA_drtr_mfe = (RNaseMRPRNA[['drtr_fn_mfe']])

tp_RNaseMRPRNA_drtr_prob = (RNaseMRPRNA[['drtr_tp_prob']])
fp_RNaseMRPRNA_drtr_prob = (RNaseMRPRNA[['drtr_fp_prob']])
tn_RNaseMRPRNA_drtr_prob = (RNaseMRPRNA[['drtr_tn_prob']])
fn_RNaseMRPRNA_drtr_prob = (RNaseMRPRNA[['drtr_fn_prob']])


#TPR
tpr_RNaseMRPRNA_fold_mfe = tpr(tp_RNaseMRPRNA_fold_mfe, fn_RNaseMRPRNA_fold_mfe)
#print(tpr_RNaseMRPRNA_fold_mfe)

tpr_RNaseMRPRNA_fold_bpp = tpr(tp_RNaseMRPRNA_fold_bpp, fn_RNaseMRPRNA_fold_bpp)
#print(tpr_RNaseMRPRNA_fold_bpp)

tpr_RNaseMRPRNA_drtr = tpr(tp_RNaseMRPRNA_drtr, fn_RNaseMRPRNA_drtr)
#print(tpr_RNaseMRPRNA_drtr)

tpr_RNaseMRPRNA_drtr_mfe = tpr(tp_RNaseMRPRNA_drtr_mfe, fn_RNaseMRPRNA_drtr_mfe)

tpr_RNaseMRPRNA_drtr_prob = tpr(tp_RNaseMRPRNA_drtr_prob, fn_RNaseMRPRNA_drtr_prob)



#PPV
ppv_RNaseMRPRNA_fold_mfe = ppv(tp_RNaseMRPRNA_fold_mfe, fp_RNaseMRPRNA_fold_mfe)
#print(ppv_RNaseMRPRNA_fold_mfe)

ppv_RNaseMRPRNA_fold_bpp = ppv(tp_RNaseMRPRNA_fold_bpp, fp_RNaseMRPRNA_fold_bpp)
#print(ppv_RNaseMRPRNA_fold_bpp)

ppv_RNaseMRPRNA_drtr = ppv(tp_RNaseMRPRNA_drtr, fp_RNaseMRPRNA_drtr)
#print(ppv_RNaseMRPRNA_drtr)

ppv_RNaseMRPRNA_drtr_mfe = ppv(tp_RNaseMRPRNA_drtr_mfe, fp_RNaseMRPRNA_drtr_mfe)

ppv_RNaseMRPRNA_drtr_prob = ppv(tp_RNaseMRPRNA_drtr_prob, fp_RNaseMRPRNA_drtr_prob)



#F1
f1_RNaseMRPRNA_fold_mfe = f_measure(tp_RNaseMRPRNA_fold_mfe, fp_RNaseMRPRNA_fold_mfe, fn_RNaseMRPRNA_fold_mfe)
#print(f1_RNaseMRPRNA_fold_mfe)

f1_RNaseMRPRNA_fold_bpp = f_measure(tp_RNaseMRPRNA_fold_bpp, fp_RNaseMRPRNA_fold_bpp, fn_RNaseMRPRNA_fold_bpp)
#print(f1_RNaseMRPRNA_fold_bpp)

f1_RNaseMRPRNA_drtr = f_measure(tp_RNaseMRPRNA_drtr, fp_RNaseMRPRNA_drtr, fn_RNaseMRPRNA_drtr)
#print(f1_RNaseMRPRNA_drtr)

f1_RNaseMRPRNA_drtr_mfe = f_measure(tp_RNaseMRPRNA_drtr_mfe, fp_RNaseMRPRNA_drtr_mfe, fn_RNaseMRPRNA_drtr_mfe)

f1_RNaseMRPRNA_drtr_prob = f_measure(tp_RNaseMRPRNA_drtr_prob, fp_RNaseMRPRNA_drtr_prob, fn_RNaseMRPRNA_drtr_prob)



#MCC
mcc_RNaseMRPRNA_fold_mfe = mcc (tp_RNaseMRPRNA_fold_mfe, fp_RNaseMRPRNA_fold_mfe, tn_RNaseMRPRNA_fold_mfe, fn_RNaseMRPRNA_fold_mfe)
#print(mcc_RNaseMRPRNA_fold_mfe)

mcc_RNaseMRPRNA_fold_bpp = mcc (tp_RNaseMRPRNA_fold_bpp, fp_RNaseMRPRNA_fold_bpp, tn_RNaseMRPRNA_fold_bpp, fn_RNaseMRPRNA_fold_bpp)
#print(mcc_RNaseMRPRNA_fold_bpp)

mcc_RNaseMRPRNA_drtr = mcc (tp_RNaseMRPRNA_drtr, fp_RNaseMRPRNA_drtr, tn_RNaseMRPRNA_drtr, fn_RNaseMRPRNA_drtr)
#print(mcc_RNaseMRPRNA_drtr)

mcc_RNaseMRPRNA_drtr_mfe = mcc (tp_RNaseMRPRNA_drtr_mfe, fp_RNaseMRPRNA_drtr_mfe, tn_RNaseMRPRNA_drtr_mfe, fn_RNaseMRPRNA_drtr_mfe)

mcc_RNaseMRPRNA_drtr_prob = mcc (tp_RNaseMRPRNA_drtr_prob, fp_RNaseMRPRNA_drtr_prob, tn_RNaseMRPRNA_drtr_prob, fn_RNaseMRPRNA_drtr_prob)



#RNasePRNA
RNasePRNA=family_data.loc[family_data.family=='RNasePRNA']

tp_RNasePRNA_fold_mfe = (RNasePRNA[['RNAfold_mfe_tp']])
fp_RNasePRNA_fold_mfe = (RNasePRNA[['RNAfold_mfe_fp']])
tn_RNasePRNA_fold_mfe = (RNasePRNA[['RNAfold_mfe_tn']])
fn_RNasePRNA_fold_mfe = (RNasePRNA[['RNAfold_mfe_fn']])

tp_RNasePRNA_fold_bpp = (RNasePRNA[['RNAfold_bpp_tp']])
fp_RNasePRNA_fold_bpp = (RNasePRNA[['RNAfold_bpp_fp']])
tn_RNasePRNA_fold_bpp = (RNasePRNA[['RNAfold_bpp_tn']])
fn_RNasePRNA_fold_bpp = (RNasePRNA[['RNAfold_bpp_fn']])

tp_RNasePRNA_drtr = (RNasePRNA[['drtr_tp']])
fp_RNasePRNA_drtr = (RNasePRNA[['drtr_fp']])
tn_RNasePRNA_drtr = (RNasePRNA[['drtr_tn']])
fn_RNasePRNA_drtr = (RNasePRNA[['drtr_fn']])

tp_RNasePRNA_drtr_mfe = (RNasePRNA[['drtr_tp_mfe']])
fp_RNasePRNA_drtr_mfe = (RNasePRNA[['drtr_fp_mfe']])
tn_RNasePRNA_drtr_mfe = (RNasePRNA[['drtr_tn_mfe']])
fn_RNasePRNA_drtr_mfe = (RNasePRNA[['drtr_fn_mfe']])

tp_RNasePRNA_drtr_prob = (RNasePRNA[['drtr_tp_prob']])
fp_RNasePRNA_drtr_prob = (RNasePRNA[['drtr_fp_prob']])
tn_RNasePRNA_drtr_prob = (RNasePRNA[['drtr_tn_prob']])
fn_RNasePRNA_drtr_prob = (RNasePRNA[['drtr_fn_prob']])


#TPR
tpr_RNasePRNA_fold_mfe = tpr(tp_RNasePRNA_fold_mfe, fn_RNasePRNA_fold_mfe)
#print(tpr_RNasePRNA_fold_mfe)

tpr_RNasePRNA_fold_bpp = tpr(tp_RNasePRNA_fold_bpp, fn_RNasePRNA_fold_bpp)
#print(tpr_RNasePRNA_fold_bpp)

tpr_RNasePRNA_drtr = tpr(tp_RNasePRNA_drtr, fn_RNasePRNA_drtr)
#print(tpr_RNasePRNA_drtr)

tpr_RNasePRNA_drtr_mfe = tpr(tp_RNasePRNA_drtr_mfe, fn_RNasePRNA_drtr_mfe)

tpr_RNasePRNA_drtr_prob = tpr(tp_RNasePRNA_drtr_prob, fn_RNasePRNA_drtr_prob)


#PPV
ppv_RNasePRNA_fold_mfe = ppv(tp_RNasePRNA_fold_mfe, fp_RNasePRNA_fold_mfe)
#print(ppv_RNasePRNA_fold_mfe)

ppv_RNasePRNA_fold_bpp = ppv(tp_RNasePRNA_fold_bpp, fp_RNasePRNA_fold_bpp)
#print(ppv_RNasePRNA_fold_bpp)

ppv_RNasePRNA_drtr = ppv(tp_RNasePRNA_drtr, fp_RNasePRNA_drtr)
#print(ppv_RNasePRNA_drtr)

ppv_RNasePRNA_drtr_mfe = ppv(tp_RNasePRNA_drtr_mfe, fp_RNasePRNA_drtr_mfe)

ppv_RNasePRNA_drtr_prob = ppv(tp_RNasePRNA_drtr_prob, fp_RNasePRNA_drtr_prob)


#F1
f1_RNasePRNA_fold_mfe = f_measure(tp_RNasePRNA_fold_mfe, fp_RNasePRNA_fold_mfe, fn_RNasePRNA_fold_mfe)
#print(f1_RNasePRNA_fold_mfe)

f1_RNasePRNA_fold_bpp = f_measure(tp_RNasePRNA_fold_bpp, fp_RNasePRNA_fold_bpp, fn_RNasePRNA_fold_bpp)
#print(f1_RNasePRNA_fold_bpp)

f1_RNasePRNA_drtr = f_measure(tp_RNasePRNA_drtr, fp_RNasePRNA_drtr, fn_RNasePRNA_drtr)
#print(f1_RNasePRNA_drtr)

f1_RNasePRNA_drtr_mfe = f_measure(tp_RNasePRNA_drtr_mfe, fp_RNasePRNA_drtr_mfe, fn_RNasePRNA_drtr_mfe)

f1_RNasePRNA_drtr_prob = f_measure(tp_RNasePRNA_drtr_prob, fp_RNasePRNA_drtr_prob, fn_RNasePRNA_drtr_prob)


#MCC
mcc_RNasePRNA_fold_mfe = mcc (tp_RNasePRNA_fold_mfe, fp_RNasePRNA_fold_mfe, tn_RNasePRNA_fold_mfe, fn_RNasePRNA_fold_mfe)
#print(mcc_RNasePRNA_fold_mfe)

mcc_RNasePRNA_fold_bpp = mcc (tp_RNasePRNA_fold_bpp, fp_RNasePRNA_fold_bpp, tn_RNasePRNA_fold_bpp, fn_RNasePRNA_fold_bpp)
#print(mcc_RNasePRNA_fold_bpp)

mcc_RNasePRNA_drtr = mcc (tp_RNasePRNA_drtr, fp_RNasePRNA_drtr, tn_RNasePRNA_drtr, fn_RNasePRNA_drtr)
#print(mcc_RNasePRNA_drtr)

mcc_RNasePRNA_drtr_mfe = mcc (tp_RNasePRNA_drtr_mfe, fp_RNasePRNA_drtr_mfe, tn_RNasePRNA_drtr_mfe, fn_RNasePRNA_drtr_mfe)

mcc_RNasePRNA_drtr_prob = mcc (tp_RNasePRNA_drtr_prob, fp_RNasePRNA_drtr_prob, tn_RNasePRNA_drtr_prob, fn_RNasePRNA_drtr_prob)




#snRNA
snRNA=family_data.loc[family_data.family=='snRNA']

tp_snRNA_fold_mfe = (snRNA[['RNAfold_mfe_tp']])
fp_snRNA_fold_mfe = (snRNA[['RNAfold_mfe_fp']])
tn_snRNA_fold_mfe = (snRNA[['RNAfold_mfe_tn']])
fn_snRNA_fold_mfe = (snRNA[['RNAfold_mfe_fn']])

tp_snRNA_fold_bpp = (snRNA[['RNAfold_bpp_tp']])
fp_snRNA_fold_bpp = (snRNA[['RNAfold_bpp_fp']])
tn_snRNA_fold_bpp = (snRNA[['RNAfold_bpp_tn']])
fn_snRNA_fold_bpp = (snRNA[['RNAfold_bpp_fn']])

tp_snRNA_drtr = (snRNA[['drtr_tp']])
fp_snRNA_drtr = (snRNA[['drtr_fp']])
tn_snRNA_drtr = (snRNA[['drtr_tn']])
fn_snRNA_drtr = (snRNA[['drtr_fn']])

tp_snRNA_drtr_mfe = (snRNA[['drtr_tp_mfe']])
fp_snRNA_drtr_mfe = (snRNA[['drtr_fp_mfe']])
tn_snRNA_drtr_mfe = (snRNA[['drtr_tn_mfe']])
fn_snRNA_drtr_mfe = (snRNA[['drtr_fn_mfe']])

tp_snRNA_drtr_prob = (snRNA[['drtr_tp_prob']])
fp_snRNA_drtr_prob = (snRNA[['drtr_fp_prob']])
tn_snRNA_drtr_prob = (snRNA[['drtr_tn_prob']])
fn_snRNA_drtr_prob = (snRNA[['drtr_fn_prob']])


#TPR
tpr_snRNA_fold_mfe = tpr(tp_snRNA_fold_mfe, fn_snRNA_fold_mfe)
#print(tpr_snRNA_fold_mfe)

tpr_snRNA_fold_bpp = tpr(tp_snRNA_fold_bpp, fn_snRNA_fold_bpp)
#print(tpr_snRNA_fold_bpp)

tpr_snRNA_drtr = tpr(tp_snRNA_drtr, fn_snRNA_drtr)
#print(tpr_snRNA_drtr)

tpr_snRNA_drtr_mfe = tpr(tp_snRNA_drtr_mfe, fn_snRNA_drtr_mfe)

tpr_snRNA_drtr_prob = tpr(tp_snRNA_drtr_prob, fn_snRNA_drtr_prob)


#PPV
ppv_snRNA_fold_mfe = ppv(tp_snRNA_fold_mfe, fp_snRNA_fold_mfe)
#print(ppv_snRNA_fold_mfe)

ppv_snRNA_fold_bpp = ppv(tp_snRNA_fold_bpp, fp_snRNA_fold_bpp)
#print(ppv_snRNA_fold_bpp)

ppv_snRNA_drtr = ppv(tp_snRNA_drtr, fp_snRNA_drtr)
#print(ppv_snRNA_drtr)

ppv_snRNA_drtr_mfe = ppv(tp_snRNA_drtr_mfe, fp_snRNA_drtr_mfe)

ppv_snRNA_drtr_prob = ppv(tp_snRNA_drtr_prob, fp_snRNA_drtr_prob)


#F1
f1_snRNA_fold_mfe = f_measure(tp_snRNA_fold_mfe, fp_snRNA_fold_mfe, fn_snRNA_fold_mfe)
#print(f1_snRNA_fold_mfe)

f1_snRNA_fold_bpp = f_measure(tp_snRNA_fold_bpp, fp_snRNA_fold_bpp, fn_snRNA_fold_bpp)
#print(f1_snRNA_fold_bpp)

f1_snRNA_drtr = f_measure(tp_snRNA_drtr, fp_snRNA_drtr, fn_snRNA_drtr)
#print(f1_snRNA_drtr)

f1_snRNA_drtr_mfe = f_measure(tp_snRNA_drtr_mfe, fp_snRNA_drtr_mfe, fn_snRNA_drtr_mfe)

f1_snRNA_drtr_prob = f_measure(tp_snRNA_drtr_prob, fp_snRNA_drtr_prob, fn_snRNA_drtr_prob)


#MCC
mcc_snRNA_fold_mfe = mcc (tp_snRNA_fold_mfe, fp_snRNA_fold_mfe, tn_snRNA_fold_mfe, fn_snRNA_fold_mfe)
#print(mcc_snRNA_fold_mfe)

mcc_snRNA_fold_bpp = mcc (tp_snRNA_fold_bpp, fp_snRNA_fold_bpp, tn_snRNA_fold_bpp, fn_snRNA_fold_bpp)
#print(mcc_snRNA_fold_bpp)

mcc_snRNA_drtr = mcc (tp_snRNA_drtr, fp_snRNA_drtr, tn_snRNA_drtr, fn_snRNA_drtr)
#print(mcc_snRNA_drtr)

mcc_snRNA_drtr_mfe = mcc (tp_snRNA_drtr_mfe, fp_snRNA_drtr_mfe, tn_snRNA_drtr_mfe, fn_snRNA_drtr_mfe)

mcc_snRNA_drtr_prob = mcc (tp_snRNA_drtr_prob, fp_snRNA_drtr_prob, tn_snRNA_drtr_prob, fn_snRNA_drtr_prob)



#SRPRNA
SRPRNA=family_data.loc[family_data.family=='SRPRNA']

tp_SRPRNA_fold_mfe = (SRPRNA[['RNAfold_mfe_tp']])
fp_SRPRNA_fold_mfe = (SRPRNA[['RNAfold_mfe_fp']])
tn_SRPRNA_fold_mfe = (SRPRNA[['RNAfold_mfe_tn']])
fn_SRPRNA_fold_mfe = (SRPRNA[['RNAfold_mfe_fn']])

tp_SRPRNA_fold_bpp = (SRPRNA[['RNAfold_bpp_tp']])
fp_SRPRNA_fold_bpp = (SRPRNA[['RNAfold_bpp_fp']])
tn_SRPRNA_fold_bpp = (SRPRNA[['RNAfold_bpp_tn']])
fn_SRPRNA_fold_bpp = (SRPRNA[['RNAfold_bpp_fn']])

tp_SRPRNA_drtr = (SRPRNA[['drtr_tp']])
fp_SRPRNA_drtr = (SRPRNA[['drtr_fp']])
tn_SRPRNA_drtr = (SRPRNA[['drtr_tn']])
fn_SRPRNA_drtr = (SRPRNA[['drtr_fn']])

tp_SRPRNA_drtr_mfe = (SRPRNA[['drtr_tp_mfe']])
fp_SRPRNA_drtr_mfe = (SRPRNA[['drtr_fp_mfe']])
tn_SRPRNA_drtr_mfe = (SRPRNA[['drtr_tn_mfe']])
fn_SRPRNA_drtr_mfe = (SRPRNA[['drtr_fn_mfe']])

tp_SRPRNA_drtr_prob = (SRPRNA[['drtr_tp_prob']])
fp_SRPRNA_drtr_prob = (SRPRNA[['drtr_fp_prob']])
tn_SRPRNA_drtr_prob = (SRPRNA[['drtr_tn_prob']])
fn_SRPRNA_drtr_prob = (SRPRNA[['drtr_fn_prob']])


#TPR
tpr_SRPRNA_fold_mfe = tpr(tp_SRPRNA_fold_mfe, fn_SRPRNA_fold_mfe)
#print(tpr_SRPRNA_fold_mfe)

tpr_SRPRNA_fold_bpp = tpr(tp_SRPRNA_fold_bpp, fn_SRPRNA_fold_bpp)
#print(tpr_SRPRNA_fold_bpp)

tpr_SRPRNA_drtr = tpr(tp_SRPRNA_drtr, fn_SRPRNA_drtr)
#print(tpr_SRPRNA_drtr)

tpr_SRPRNA_drtr_mfe = tpr(tp_SRPRNA_drtr_mfe, fn_SRPRNA_drtr_mfe)

tpr_SRPRNA_drtr_prob = tpr(tp_SRPRNA_drtr_prob, fn_SRPRNA_drtr_prob)



#PPV
ppv_SRPRNA_fold_mfe = ppv(tp_SRPRNA_fold_mfe, fp_SRPRNA_fold_mfe)
#print(ppv_SRPRNA_fold_mfe)

ppv_SRPRNA_fold_bpp = ppv(tp_SRPRNA_fold_bpp, fp_SRPRNA_fold_bpp)
#print(ppv_SRPRNA_fold_bpp)

ppv_SRPRNA_drtr = ppv(tp_SRPRNA_drtr, fp_SRPRNA_drtr)
#print(ppv_SRPRNA_drtr)

ppv_SRPRNA_drtr_mfe = ppv(tp_SRPRNA_drtr_mfe, fp_SRPRNA_drtr_mfe)

ppv_SRPRNA_drtr_prob = ppv(tp_SRPRNA_drtr_prob, fp_SRPRNA_drtr_prob)



#F1
f1_SRPRNA_fold_mfe = f_measure(tp_SRPRNA_fold_mfe, fp_SRPRNA_fold_mfe, fn_SRPRNA_fold_mfe)
#print(f1_SRPRNA_fold_mfe)

f1_SRPRNA_fold_bpp = f_measure(tp_SRPRNA_fold_bpp, fp_SRPRNA_fold_bpp, fn_SRPRNA_fold_bpp)
#print(f1_SRPRNA_fold_bpp)

f1_SRPRNA_drtr = f_measure(tp_SRPRNA_drtr, fp_SRPRNA_drtr, fn_SRPRNA_drtr)
#print(f1_SRPRNA_drtr)

f1_SRPRNA_drtr_mfe = f_measure(tp_SRPRNA_drtr_mfe, fp_SRPRNA_drtr_mfe, fn_SRPRNA_drtr_mfe)

f1_SRPRNA_drtr_prob = f_measure(tp_SRPRNA_drtr_prob, fp_SRPRNA_drtr_prob, fn_SRPRNA_drtr_prob)


#MCC
mcc_SRPRNA_fold_mfe = mcc (tp_SRPRNA_fold_mfe, fp_SRPRNA_fold_mfe, tn_SRPRNA_fold_mfe, fn_SRPRNA_fold_mfe)
#print(mcc_SRPRNA_fold_mfe)

mcc_SRPRNA_fold_bpp = mcc (tp_SRPRNA_fold_bpp, fp_SRPRNA_fold_bpp, tn_SRPRNA_fold_bpp, fn_SRPRNA_fold_bpp)
#print(mcc_SRPRNA_fold_bpp)

mcc_SRPRNA_drtr = mcc (tp_SRPRNA_drtr, fp_SRPRNA_drtr, tn_SRPRNA_drtr, fn_SRPRNA_drtr)
#print(mcc_SRPRNA_drtr)

mcc_SRPRNA_drtr_mfe = mcc (tp_SRPRNA_drtr_mfe, fp_SRPRNA_drtr_mfe, tn_SRPRNA_drtr_mfe, fn_SRPRNA_drtr_mfe)

mcc_SRPRNA_drtr_prob = mcc (tp_SRPRNA_drtr_prob, fp_SRPRNA_drtr_prob, tn_SRPRNA_drtr_prob, fn_SRPRNA_drtr_prob)



#SyntheticRNA
SyntheticRNA=family_data.loc[family_data.family=='SyntheticRNA']

tp_SyntheticRNA_fold_mfe = (SyntheticRNA[['RNAfold_mfe_tp']])
fp_SyntheticRNA_fold_mfe = (SyntheticRNA[['RNAfold_mfe_fp']])
tn_SyntheticRNA_fold_mfe = (SyntheticRNA[['RNAfold_mfe_tn']])
fn_SyntheticRNA_fold_mfe = (SyntheticRNA[['RNAfold_mfe_fn']])

tp_SyntheticRNA_fold_bpp = (SyntheticRNA[['RNAfold_bpp_tp']])
fp_SyntheticRNA_fold_bpp = (SyntheticRNA[['RNAfold_bpp_fp']])
tn_SyntheticRNA_fold_bpp = (SyntheticRNA[['RNAfold_bpp_tn']])
fn_SyntheticRNA_fold_bpp = (SyntheticRNA[['RNAfold_bpp_fn']])

tp_SyntheticRNA_drtr = (SyntheticRNA[['drtr_tp']])
fp_SyntheticRNA_drtr = (SyntheticRNA[['drtr_fp']])
tn_SyntheticRNA_drtr = (SyntheticRNA[['drtr_tn']])
fn_SyntheticRNA_drtr = (SyntheticRNA[['drtr_fn']])

tp_SyntheticRNA_drtr_mfe = (SyntheticRNA[['drtr_tp_mfe']])
fp_SyntheticRNA_drtr_mfe = (SyntheticRNA[['drtr_fp_mfe']])
tn_SyntheticRNA_drtr_mfe = (SyntheticRNA[['drtr_tn_mfe']])
fn_SyntheticRNA_drtr_mfe = (SyntheticRNA[['drtr_fn_mfe']])

tp_SyntheticRNA_drtr_prob = (SyntheticRNA[['drtr_tp_prob']])
fp_SyntheticRNA_drtr_prob = (SyntheticRNA[['drtr_fp_prob']])
tn_SyntheticRNA_drtr_prob = (SyntheticRNA[['drtr_tn_prob']])
fn_SyntheticRNA_drtr_prob = (SyntheticRNA[['drtr_fn_prob']])


#TPR
tpr_SyntheticRNA_fold_mfe = tpr(tp_SyntheticRNA_fold_mfe, fn_SyntheticRNA_fold_mfe)
#print(tpr_SyntheticRNA_fold_mfe)

tpr_SyntheticRNA_fold_bpp = tpr(tp_SyntheticRNA_fold_bpp, fn_SyntheticRNA_fold_bpp)
#print(tpr_SyntheticRNA_fold_bpp)

tpr_SyntheticRNA_drtr = tpr(tp_SyntheticRNA_drtr, fn_SyntheticRNA_drtr)
#print(tpr_SyntheticRNA_drtr)

tpr_SyntheticRNA_drtr_mfe = tpr(tp_SyntheticRNA_drtr_mfe, fn_SyntheticRNA_drtr_mfe)

tpr_SyntheticRNA_drtr_prob = tpr(tp_SyntheticRNA_drtr_prob, fn_SyntheticRNA_drtr_prob)


#PPV
ppv_SyntheticRNA_fold_mfe = ppv(tp_SyntheticRNA_fold_mfe, fp_SyntheticRNA_fold_mfe)
#print(ppv_SyntheticRNA_fold_mfe)

ppv_SyntheticRNA_fold_bpp = ppv(tp_SyntheticRNA_fold_bpp, fp_SyntheticRNA_fold_bpp)
#print(ppv_SyntheticRNA_fold_bpp)

ppv_SyntheticRNA_drtr = ppv(tp_SyntheticRNA_drtr, fp_SyntheticRNA_drtr)
#print(ppv_SyntheticRNA_drtr)

ppv_SyntheticRNA_drtr_mfe = ppv(tp_SyntheticRNA_drtr_mfe, fp_SyntheticRNA_drtr_mfe)

ppv_SyntheticRNA_drtr_prob = ppv(tp_SyntheticRNA_drtr_prob, fp_SyntheticRNA_drtr_prob)


#F1
f1_SyntheticRNA_fold_mfe = f_measure(tp_SyntheticRNA_fold_mfe, fp_SyntheticRNA_fold_mfe, fn_SyntheticRNA_fold_mfe)
#print(f1_SyntheticRNA_fold_mfe)

f1_SyntheticRNA_fold_bpp = f_measure(tp_SyntheticRNA_fold_bpp, fp_SyntheticRNA_fold_bpp, fn_SyntheticRNA_fold_bpp)
#print(f1_SyntheticRNA_fold_bpp)

f1_SyntheticRNA_drtr = f_measure(tp_SyntheticRNA_drtr, fp_SyntheticRNA_drtr, fn_SyntheticRNA_drtr)
#print(f1_SyntheticRNA_drtr)

f1_SyntheticRNA_drtr_mfe = f_measure(tp_SyntheticRNA_drtr_mfe, fp_SyntheticRNA_drtr_mfe, fn_SyntheticRNA_drtr_mfe)

f1_SyntheticRNA_drtr_prob = f_measure(tp_SyntheticRNA_drtr_prob, fp_SyntheticRNA_drtr_prob, fn_SyntheticRNA_drtr_prob)


#MCC
mcc_SyntheticRNA_fold_mfe = mcc (tp_SyntheticRNA_fold_mfe, fp_SyntheticRNA_fold_mfe, tn_SyntheticRNA_fold_mfe, fn_SyntheticRNA_fold_mfe)
#print(mcc_SyntheticRNA_fold_mfe)

mcc_SyntheticRNA_fold_bpp = mcc (tp_SyntheticRNA_fold_bpp, fp_SyntheticRNA_fold_bpp, tn_SyntheticRNA_fold_bpp, fn_SyntheticRNA_fold_bpp)
#print(mcc_SyntheticRNA_fold_bpp)

mcc_SyntheticRNA_drtr = mcc (tp_SyntheticRNA_drtr, fp_SyntheticRNA_drtr, tn_SyntheticRNA_drtr, fn_SyntheticRNA_drtr)
#print(mcc_SyntheticRNA_drtr)

mcc_SyntheticRNA_drtr_mfe = mcc (tp_SyntheticRNA_drtr_mfe, fp_SyntheticRNA_drtr_mfe, tn_SyntheticRNA_drtr_mfe, fn_SyntheticRNA_drtr_mfe)

mcc_SyntheticRNA_drtr_prob = mcc (tp_SyntheticRNA_drtr_prob, fp_SyntheticRNA_drtr_prob, tn_SyntheticRNA_drtr_prob, fn_SyntheticRNA_drtr_prob)


#tmRNA
tmRNA=family_data.loc[family_data.family=='tmRNA']

tp_tmRNA_fold_mfe = (tmRNA[['RNAfold_mfe_tp']])
fp_tmRNA_fold_mfe = (tmRNA[['RNAfold_mfe_fp']])
tn_tmRNA_fold_mfe = (tmRNA[['RNAfold_mfe_tn']])
fn_tmRNA_fold_mfe = (tmRNA[['RNAfold_mfe_fn']])

tp_tmRNA_fold_bpp = (tmRNA[['RNAfold_bpp_tp']])
fp_tmRNA_fold_bpp = (tmRNA[['RNAfold_bpp_fp']])
tn_tmRNA_fold_bpp = (tmRNA[['RNAfold_bpp_tn']])
fn_tmRNA_fold_bpp = (tmRNA[['RNAfold_bpp_fn']])

tp_tmRNA_drtr = (tmRNA[['drtr_tp']])
fp_tmRNA_drtr = (tmRNA[['drtr_fp']])
tn_tmRNA_drtr = (tmRNA[['drtr_tn']])
fn_tmRNA_drtr = (tmRNA[['drtr_fn']])

tp_tmRNA_drtr_mfe = (tmRNA[['drtr_tp_mfe']])
fp_tmRNA_drtr_mfe = (tmRNA[['drtr_fp_mfe']])
tn_tmRNA_drtr_mfe = (tmRNA[['drtr_tn_mfe']])
fn_tmRNA_drtr_mfe = (tmRNA[['drtr_fn_mfe']])

tp_tmRNA_drtr_prob = (tmRNA[['drtr_tp_prob']])
fp_tmRNA_drtr_prob = (tmRNA[['drtr_fp_prob']])
tn_tmRNA_drtr_prob = (tmRNA[['drtr_tn_prob']])
fn_tmRNA_drtr_prob = (tmRNA[['drtr_fn_prob']])


#TPR
tpr_tmRNA_fold_mfe = tpr(tp_tmRNA_fold_mfe, fn_tmRNA_fold_mfe)
#print(tpr_tmRNA_fold_mfe)

tpr_tmRNA_fold_bpp = tpr(tp_tmRNA_fold_bpp, fn_tmRNA_fold_bpp)
#print(tpr_tmRNA_fold_bpp)

tpr_tmRNA_drtr = tpr(tp_tmRNA_drtr, fn_tmRNA_drtr)
#print(tpr_tmRNA_drtr)

tpr_tmRNA_drtr_mfe = tpr(tp_tmRNA_drtr_mfe, fn_tmRNA_drtr_mfe)

tpr_tmRNA_drtr_prob = tpr(tp_tmRNA_drtr_prob, fn_tmRNA_drtr_prob)


#PPV
ppv_tmRNA_fold_mfe = ppv(tp_tmRNA_fold_mfe, fp_tmRNA_fold_mfe)
#print(ppv_tmRNA_fold_mfe)

ppv_tmRNA_fold_bpp = ppv(tp_tmRNA_fold_bpp, fp_tmRNA_fold_bpp)
#print(ppv_tmRNA_fold_bpp)

ppv_tmRNA_drtr = ppv(tp_tmRNA_drtr, fp_tmRNA_drtr)
#print(ppv_tmRNA_drtr)

ppv_tmRNA_drtr_mfe = ppv(tp_tmRNA_drtr_mfe, fp_tmRNA_drtr_mfe)

ppv_tmRNA_drtr_prob = ppv(tp_tmRNA_drtr_prob, fp_tmRNA_drtr_prob)


#F1
f1_tmRNA_fold_mfe = f_measure(tp_tmRNA_fold_mfe, fp_tmRNA_fold_mfe, fn_tmRNA_fold_mfe)
#print(f1_tmRNA_fold_mfe)

f1_tmRNA_fold_bpp = f_measure(tp_tmRNA_fold_bpp, fp_tmRNA_fold_bpp, fn_tmRNA_fold_bpp)
#print(f1_tmRNA_fold_bpp)

f1_tmRNA_drtr = f_measure(tp_tmRNA_drtr, fp_tmRNA_drtr, fn_tmRNA_drtr)
#print(f1_tmRNA_drtr)

f1_tmRNA_drtr_mfe = f_measure(tp_tmRNA_drtr_mfe, fp_tmRNA_drtr_mfe, fn_tmRNA_drtr_mfe)

f1_tmRNA_drtr_prob = f_measure(tp_tmRNA_drtr_prob, fp_tmRNA_drtr_prob, fn_tmRNA_drtr_prob)


#MCC
mcc_tmRNA_fold_mfe = mcc (tp_tmRNA_fold_mfe, fp_tmRNA_fold_mfe, tn_tmRNA_fold_mfe, fn_tmRNA_fold_mfe)
#print(mcc_tmRNA_fold_mfe)

mcc_tmRNA_fold_bpp = mcc (tp_tmRNA_fold_bpp, fp_tmRNA_fold_bpp, tn_tmRNA_fold_bpp, fn_tmRNA_fold_bpp)
#print(mcc_tmRNA_fold_bpp)

mcc_tmRNA_drtr = mcc (tp_tmRNA_drtr, fp_tmRNA_drtr, tn_tmRNA_drtr, fn_tmRNA_drtr)
#print(mcc_tmRNA_drtr)

mcc_tmRNA_drtr_mfe = mcc (tp_tmRNA_drtr_mfe, fp_tmRNA_drtr_mfe, tn_tmRNA_drtr_mfe, fn_tmRNA_drtr_mfe)

mcc_tmRNA_drtr_prob = mcc (tp_tmRNA_drtr_prob, fp_tmRNA_drtr_prob, tn_tmRNA_drtr_prob, fn_tmRNA_drtr_prob)


#tRNA
tRNA=family_data.loc[family_data.family=='tRNA']

tp_tRNA_fold_mfe = (tRNA[['RNAfold_mfe_tp']])
fp_tRNA_fold_mfe = (tRNA[['RNAfold_mfe_fp']])
tn_tRNA_fold_mfe = (tRNA[['RNAfold_mfe_tn']])
fn_tRNA_fold_mfe = (tRNA[['RNAfold_mfe_fn']])

tp_tRNA_fold_bpp = (tRNA[['RNAfold_bpp_tp']])
fp_tRNA_fold_bpp = (tRNA[['RNAfold_bpp_fp']])
tn_tRNA_fold_bpp = (tRNA[['RNAfold_bpp_tn']])
fn_tRNA_fold_bpp = (tRNA[['RNAfold_bpp_fn']])

tp_tRNA_drtr = (tRNA[['drtr_tp']])
fp_tRNA_drtr = (tRNA[['drtr_fp']])
tn_tRNA_drtr = (tRNA[['drtr_tn']])
fn_tRNA_drtr = (tRNA[['drtr_fn']])

tp_tRNA_drtr_mfe = (tRNA[['drtr_tp_mfe']])
fp_tRNA_drtr_mfe = (tRNA[['drtr_fp_mfe']])
tn_tRNA_drtr_mfe = (tRNA[['drtr_tn_mfe']])
fn_tRNA_drtr_mfe = (tRNA[['drtr_fn_mfe']])

tp_tRNA_drtr_prob = (tRNA[['drtr_tp_prob']])
fp_tRNA_drtr_prob = (tRNA[['drtr_fp_prob']])
tn_tRNA_drtr_prob = (tRNA[['drtr_tn_prob']])
fn_tRNA_drtr_prob = (tRNA[['drtr_fn_prob']])

#TPR
tpr_tRNA_fold_mfe = tpr(tp_tRNA_fold_mfe, fn_tRNA_fold_mfe)
#print(tpr_tRNA_fold_mfe)

tpr_tRNA_fold_bpp = tpr(tp_tRNA_fold_bpp, fn_tRNA_fold_bpp)
#print(tpr_tRNA_fold_bpp)

tpr_tRNA_drtr = tpr(tp_tRNA_drtr, fn_tRNA_drtr)
#print(tpr_tRNA_drtr)

tpr_tRNA_drtr_mfe = tpr(tp_tRNA_drtr_mfe, fn_tRNA_drtr_mfe)

tpr_tRNA_drtr_prob = tpr(tp_tRNA_drtr_prob, fn_tRNA_drtr_prob)

#PPV
ppv_tRNA_fold_mfe = ppv(tp_tRNA_fold_mfe, fp_tRNA_fold_mfe)
#print(ppv_tRNA_fold_mfe)

ppv_tRNA_fold_bpp = ppv(tp_tRNA_fold_bpp, fp_tRNA_fold_bpp)
#print(ppv_tRNA_fold_bpp)

ppv_tRNA_drtr = ppv(tp_tRNA_drtr, fp_tRNA_drtr)
#print(ppv_tRNA_drtr)

ppv_tRNA_drtr_mfe = ppv(tp_tRNA_drtr_mfe, fp_tRNA_drtr_mfe)

ppv_tRNA_drtr_prob = ppv(tp_tRNA_drtr_prob, fp_tRNA_drtr_prob)


#F1
f1_tRNA_fold_mfe = f_measure(tp_tRNA_fold_mfe, fp_tRNA_fold_mfe, fn_tRNA_fold_mfe)
#print(f1_tRNA_fold_mfe)

f1_tRNA_fold_bpp = f_measure(tp_tRNA_fold_bpp, fp_tRNA_fold_bpp, fn_tRNA_fold_bpp)
#print(f1_tRNA_fold_bpp)

f1_tRNA_drtr = f_measure(tp_tRNA_drtr, fp_tRNA_drtr, fn_tRNA_drtr)
#print(f1_tRNA_drtr)

f1_tRNA_drtr_mfe = f_measure(tp_tRNA_drtr_mfe, fp_tRNA_drtr_mfe, fn_tRNA_drtr_mfe)

f1_tRNA_drtr_prob = f_measure(tp_tRNA_drtr_prob, fp_tRNA_drtr_prob, fn_tRNA_drtr_prob)


#MCC
mcc_tRNA_fold_mfe = mcc (tp_tRNA_fold_mfe, fp_tRNA_fold_mfe, tn_tRNA_fold_mfe, fn_tRNA_fold_mfe)
#print(mcc_tRNA_fold_mfe)

mcc_tRNA_fold_bpp = mcc (tp_tRNA_fold_bpp, fp_tRNA_fold_bpp, tn_tRNA_fold_bpp, fn_tRNA_fold_bpp)
#print(mcc_tRNA_fold_bpp)

mcc_tRNA_drtr = mcc (tp_tRNA_drtr, fp_tRNA_drtr, tn_tRNA_drtr, fn_tRNA_drtr)
#print(mcc_tRNA_drtr)

mcc_tRNA_drtr_mfe = mcc (tp_tRNA_drtr_mfe, fp_tRNA_drtr_mfe, tn_tRNA_drtr_mfe, fn_tRNA_drtr_mfe)

mcc_tRNA_drtr_prob = mcc (tp_tRNA_drtr_prob, fp_tRNA_drtr_prob, tn_tRNA_drtr_prob, fn_tRNA_drtr_prob)


#Vert.Telo. RNA
Vert_Telo_RNA=family_data.loc[family_data.family=='Vert.Telo. RNA']

tp_Vert_Telo_RNA_fold_mfe = (Vert_Telo_RNA[['RNAfold_mfe_tp']])
fp_Vert_Telo_RNA_fold_mfe = (Vert_Telo_RNA[['RNAfold_mfe_fp']])
tn_Vert_Telo_RNA_fold_mfe = (Vert_Telo_RNA[['RNAfold_mfe_tn']])
fn_Vert_Telo_RNA_fold_mfe = (Vert_Telo_RNA[['RNAfold_mfe_fn']])

tp_Vert_Telo_RNA_fold_bpp = (Vert_Telo_RNA[['RNAfold_bpp_tp']])
fp_Vert_Telo_RNA_fold_bpp = (Vert_Telo_RNA[['RNAfold_bpp_fp']])
tn_Vert_Telo_RNA_fold_bpp = (Vert_Telo_RNA[['RNAfold_bpp_tn']])
fn_Vert_Telo_RNA_fold_bpp = (Vert_Telo_RNA[['RNAfold_bpp_fn']])

tp_Vert_Telo_RNA_drtr = (Vert_Telo_RNA[['drtr_tp']])
fp_Vert_Telo_RNA_drtr = (Vert_Telo_RNA[['drtr_fp']])
tn_Vert_Telo_RNA_drtr = (Vert_Telo_RNA[['drtr_tn']])
fn_Vert_Telo_RNA_drtr = (Vert_Telo_RNA[['drtr_fn']])

tp_Vert_Telo_RNA_drtr_mfe = (Vert_Telo_RNA[['drtr_tp_mfe']])
fp_Vert_Telo_RNA_drtr_mfe = (Vert_Telo_RNA[['drtr_fp_mfe']])
tn_Vert_Telo_RNA_drtr_mfe = (Vert_Telo_RNA[['drtr_tn_mfe']])
fn_Vert_Telo_RNA_drtr_mfe = (Vert_Telo_RNA[['drtr_fn_mfe']])

tp_Vert_Telo_RNA_drtr_prob = (Vert_Telo_RNA[['drtr_tp_prob']])
fp_Vert_Telo_RNA_drtr_prob = (Vert_Telo_RNA[['drtr_fp_prob']])
tn_Vert_Telo_RNA_drtr_prob = (Vert_Telo_RNA[['drtr_tn_prob']])
fn_Vert_Telo_RNA_drtr_prob = (Vert_Telo_RNA[['drtr_fn_prob']])


#TPR
tpr_Vert_Telo_RNA_fold_mfe = tpr(tp_Vert_Telo_RNA_fold_mfe, fn_Vert_Telo_RNA_fold_mfe)
#print(tpr_Vert_Telo_RNA_fold_mfe)

tpr_Vert_Telo_RNA_fold_bpp = tpr(tp_Vert_Telo_RNA_fold_bpp, fn_Vert_Telo_RNA_fold_bpp)
#print(tpr_Vert_Telo_RNA_fold_bpp)

tpr_Vert_Telo_RNA_drtr = tpr(tp_Vert_Telo_RNA_drtr, fn_Vert_Telo_RNA_drtr)
#print(tpr_Vert_Telo_RNA_drtr)

tpr_Vert_Telo_RNA_drtr_mfe = tpr(tp_Vert_Telo_RNA_drtr_mfe, fn_Vert_Telo_RNA_drtr_mfe)

tpr_Vert_Telo_RNA_drtr_prob = tpr(tp_Vert_Telo_RNA_drtr_prob, fn_Vert_Telo_RNA_drtr_prob)


#PPV
ppv_Vert_Telo_RNA_fold_mfe = ppv(tp_Vert_Telo_RNA_fold_mfe, fp_Vert_Telo_RNA_fold_mfe)
#print(ppv_Vert_Telo_RNA_fold_mfe)

ppv_Vert_Telo_RNA_fold_bpp = ppv(tp_Vert_Telo_RNA_fold_bpp, fp_Vert_Telo_RNA_fold_bpp)
#print(ppv_Vert_Telo_RNA_fold_bpp)

ppv_Vert_Telo_RNA_drtr = ppv(tp_Vert_Telo_RNA_drtr, fp_Vert_Telo_RNA_drtr)
#print(ppv_Vert_Telo_RNA_drtr)

ppv_Vert_Telo_RNA_drtr_mfe = ppv(tp_Vert_Telo_RNA_drtr_mfe, fp_Vert_Telo_RNA_drtr_mfe)

ppv_Vert_Telo_RNA_drtr_prob = ppv(tp_Vert_Telo_RNA_drtr_prob, fp_Vert_Telo_RNA_drtr_prob)


#F1
f1_Vert_Telo_RNA_fold_mfe = f_measure(tp_Vert_Telo_RNA_fold_mfe, fp_Vert_Telo_RNA_fold_mfe, fn_Vert_Telo_RNA_fold_mfe)
#print(f1_Vert_Telo_RNA_fold_mfe)

f1_Vert_Telo_RNA_fold_bpp = f_measure(tp_Vert_Telo_RNA_fold_bpp, fp_Vert_Telo_RNA_fold_bpp, fn_Vert_Telo_RNA_fold_bpp)
#print(f1_Vert_Telo_RNA_fold_bpp)

f1_Vert_Telo_RNA_drtr = f_measure(tp_Vert_Telo_RNA_drtr, fp_Vert_Telo_RNA_drtr, fn_Vert_Telo_RNA_drtr)
#print(f1_Vert_Telo_RNA_drtr)

f1_Vert_Telo_RNA_drtr_mfe = f_measure(tp_Vert_Telo_RNA_drtr_mfe, fp_Vert_Telo_RNA_drtr_mfe, fn_Vert_Telo_RNA_drtr_mfe)

f1_Vert_Telo_RNA_drtr_prob = f_measure(tp_Vert_Telo_RNA_drtr_prob, fp_Vert_Telo_RNA_drtr_prob, fn_Vert_Telo_RNA_drtr_prob)


#MCC
mcc_Vert_Telo_RNA_fold_mfe = mcc (tp_Vert_Telo_RNA_fold_mfe, fp_Vert_Telo_RNA_fold_mfe, tn_Vert_Telo_RNA_fold_mfe, fn_Vert_Telo_RNA_fold_mfe)
#print(mcc_Vert_Telo_RNA_fold_mfe)

mcc_Vert_Telo_RNA_fold_bpp = mcc (tp_Vert_Telo_RNA_fold_bpp, fp_Vert_Telo_RNA_fold_bpp, tn_Vert_Telo_RNA_fold_bpp, fn_Vert_Telo_RNA_fold_bpp)
#print(mcc_Vert_Telo_RNA_fold_bpp)

mcc_Vert_Telo_RNA_drtr = mcc (tp_Vert_Telo_RNA_drtr, fp_Vert_Telo_RNA_drtr, tn_Vert_Telo_RNA_drtr, fn_Vert_Telo_RNA_drtr)
#print(mcc_Vert_Telo_RNA_drtr)

mcc_Vert_Telo_RNA_drtr_mfe = mcc (tp_Vert_Telo_RNA_drtr_mfe, fp_Vert_Telo_RNA_drtr_mfe, tn_Vert_Telo_RNA_drtr_mfe, fn_Vert_Telo_RNA_drtr_mfe)

mcc_Vert_Telo_RNA_drtr_prob = mcc (tp_Vert_Telo_RNA_drtr_prob, fp_Vert_Telo_RNA_drtr_prob, tn_Vert_Telo_RNA_drtr_prob, fn_Vert_Telo_RNA_drtr_prob)


#Viral&Phage
Viral_Phage=family_data.loc[family_data.family=='Viral&Phage']

tp_Viral_Phage_fold_mfe = (Viral_Phage[['RNAfold_mfe_tp']])
fp_Viral_Phage_fold_mfe = (Viral_Phage[['RNAfold_mfe_fp']])
tn_Viral_Phage_fold_mfe = (Viral_Phage[['RNAfold_mfe_tn']])
fn_Viral_Phage_fold_mfe = (Viral_Phage[['RNAfold_mfe_fn']])

tp_Viral_Phage_fold_bpp = (Viral_Phage[['RNAfold_bpp_tp']])
fp_Viral_Phage_fold_bpp = (Viral_Phage[['RNAfold_bpp_fp']])
tn_Viral_Phage_fold_bpp = (Viral_Phage[['RNAfold_bpp_tn']])
fn_Viral_Phage_fold_bpp = (Viral_Phage[['RNAfold_bpp_fn']])

tp_Viral_Phage_drtr = (Viral_Phage[['drtr_tp']])
fp_Viral_Phage_drtr = (Viral_Phage[['drtr_fp']])
tn_Viral_Phage_drtr = (Viral_Phage[['drtr_tn']])
fn_Viral_Phage_drtr = (Viral_Phage[['drtr_fn']])

tp_Viral_Phage_drtr_mfe = (Viral_Phage[['drtr_tp_mfe']])
fp_Viral_Phage_drtr_mfe = (Viral_Phage[['drtr_fp_mfe']])
tn_Viral_Phage_drtr_mfe = (Viral_Phage[['drtr_tn_mfe']])
fn_Viral_Phage_drtr_mfe = (Viral_Phage[['drtr_fn_mfe']])

tp_Viral_Phage_drtr_prob = (Viral_Phage[['drtr_tp_prob']])
fp_Viral_Phage_drtr_prob = (Viral_Phage[['drtr_fp_prob']])
tn_Viral_Phage_drtr_prob = (Viral_Phage[['drtr_tn_prob']])
fn_Viral_Phage_drtr_prob = (Viral_Phage[['drtr_fn_prob']])


#TPR
tpr_Viral_Phage_fold_mfe = tpr(tp_Viral_Phage_fold_mfe, fn_Viral_Phage_fold_mfe)
#print(tpr_Viral_Phage_fold_mfe)

tpr_Viral_Phage_fold_bpp = tpr(tp_Viral_Phage_fold_bpp, fn_Viral_Phage_fold_bpp)
#print(tpr_Viral_Phage_fold_bpp)

tpr_Viral_Phage_drtr = tpr(tp_Viral_Phage_drtr, fn_Viral_Phage_drtr)
#print(tpr_Viral_Phage_drtr)

tpr_Viral_Phage_drtr_mfe = tpr(tp_Viral_Phage_drtr_mfe, fn_Viral_Phage_drtr_mfe)

tpr_Viral_Phage_drtr_prob = tpr(tp_Viral_Phage_drtr_prob, fn_Viral_Phage_drtr_prob)



#PPV
ppv_Viral_Phage_fold_mfe = ppv(tp_Viral_Phage_fold_mfe, fp_Viral_Phage_fold_mfe)
#print(ppv_Viral_Phage_fold_mfe)

ppv_Viral_Phage_fold_bpp = ppv(tp_Viral_Phage_fold_bpp, fp_Viral_Phage_fold_bpp)
#print(ppv_Viral_Phage_fold_bpp)

ppv_Viral_Phage_drtr = ppv(tp_Viral_Phage_drtr, fp_Viral_Phage_drtr)
#print(ppv_Viral_Phage_drtr)

ppv_Viral_Phage_drtr_mfe = ppv(tp_Viral_Phage_drtr_mfe, fp_Viral_Phage_drtr_mfe)

ppv_Viral_Phage_drtr_prob = ppv(tp_Viral_Phage_drtr_prob, fp_Viral_Phage_drtr_prob)


#F1
f1_Viral_Phage_fold_mfe = f_measure(tp_Viral_Phage_fold_mfe, fp_Viral_Phage_fold_mfe, fn_Viral_Phage_fold_mfe)
#print(f1_Viral_Phage_fold_mfe)

f1_Viral_Phage_fold_bpp = f_measure(tp_Viral_Phage_fold_bpp, fp_Viral_Phage_fold_bpp, fn_Viral_Phage_fold_bpp)
#print(f1_Viral_Phage_fold_bpp)

f1_Viral_Phage_drtr = f_measure(tp_Viral_Phage_drtr, fp_Viral_Phage_drtr, fn_Viral_Phage_drtr)
#print(f1_Viral_Phage_drtr)

f1_Viral_Phage_drtr_mfe = f_measure(tp_Viral_Phage_drtr_mfe, fp_Viral_Phage_drtr_mfe, fn_Viral_Phage_drtr_mfe)

f1_Viral_Phage_drtr_prob = f_measure(tp_Viral_Phage_drtr_prob, fp_Viral_Phage_drtr_prob, fn_Viral_Phage_drtr_prob)


#MCC
mcc_Viral_Phage_fold_mfe = mcc (tp_Viral_Phage_fold_mfe, fp_Viral_Phage_fold_mfe, tn_Viral_Phage_fold_mfe, fn_Viral_Phage_fold_mfe)
#print(mcc_Viral_Phage_fold_mfe)

mcc_Viral_Phage_fold_bpp = mcc (tp_Viral_Phage_fold_bpp, fp_Viral_Phage_fold_bpp, tn_Viral_Phage_fold_bpp, fn_Viral_Phage_fold_bpp)
#print(mcc_Viral_Phage_fold_bpp)

mcc_Viral_Phage_drtr = mcc (tp_Viral_Phage_drtr, fp_Viral_Phage_drtr, tn_Viral_Phage_drtr, fn_Viral_Phage_drtr)
#print(mcc_Viral_Phage_drtr)

mcc_Viral_Phage_drtr_mfe = mcc (tp_Viral_Phage_drtr_mfe, fp_Viral_Phage_drtr_mfe, tn_Viral_Phage_drtr_mfe, fn_Viral_Phage_drtr_mfe)

mcc_Viral_Phage_drtr_prob = mcc (tp_Viral_Phage_drtr_prob, fp_Viral_Phage_drtr_prob, tn_Viral_Phage_drtr_prob, fn_Viral_Phage_drtr_prob)


#YRNA
YRNA=family_data.loc[family_data.family=='YRNA']

tp_YRNA_fold_mfe = (YRNA[['RNAfold_mfe_tp']])
fp_YRNA_fold_mfe = (YRNA[['RNAfold_mfe_fp']])
tn_YRNA_fold_mfe = (YRNA[['RNAfold_mfe_tn']])
fn_YRNA_fold_mfe = (YRNA[['RNAfold_mfe_fn']])

tp_YRNA_fold_bpp = (YRNA[['RNAfold_bpp_tp']])
fp_YRNA_fold_bpp = (YRNA[['RNAfold_bpp_fp']])
tn_YRNA_fold_bpp = (YRNA[['RNAfold_bpp_tn']])
fn_YRNA_fold_bpp = (YRNA[['RNAfold_bpp_fn']])

tp_YRNA_drtr = (YRNA[['drtr_tp']])
fp_YRNA_drtr = (YRNA[['drtr_fp']])
tn_YRNA_drtr = (YRNA[['drtr_tn']])
fn_YRNA_drtr = (YRNA[['drtr_fn']])

tp_YRNA_drtr_mfe = (YRNA[['drtr_tp_mfe']])
fp_YRNA_drtr_mfe = (YRNA[['drtr_fp_mfe']])
tn_YRNA_drtr_mfe = (YRNA[['drtr_tn_mfe']])
fn_YRNA_drtr_mfe = (YRNA[['drtr_fn_mfe']])

tp_YRNA_drtr_prob = (YRNA[['drtr_tp_prob']])
fp_YRNA_drtr_prob = (YRNA[['drtr_fp_prob']])
tn_YRNA_drtr_prob = (YRNA[['drtr_tn_prob']])
fn_YRNA_drtr_prob = (YRNA[['drtr_fn_prob']])


#TPR
tpr_YRNA_fold_mfe = tpr(tp_YRNA_fold_mfe, fn_YRNA_fold_mfe)
#print(tpr_YRNA_fold_mfe)

tpr_YRNA_fold_bpp = tpr(tp_YRNA_fold_bpp, fn_YRNA_fold_bpp)
#print(tpr_YRNA_fold_bpp)

tpr_YRNA_drtr = tpr(tp_YRNA_drtr, fn_YRNA_drtr)
#print(tpr_YRNA_drtr)

tpr_YRNA_drtr_mfe = tpr(tp_YRNA_drtr_mfe, fn_YRNA_drtr_mfe)

tpr_YRNA_drtr_prob = tpr(tp_YRNA_drtr_prob, fn_YRNA_drtr_prob)


#PPV
ppv_YRNA_fold_mfe = ppv(tp_YRNA_fold_mfe, fp_YRNA_fold_mfe)
#print(ppv_YRNA_fold_mfe)

ppv_YRNA_fold_bpp = ppv(tp_YRNA_fold_bpp, fp_YRNA_fold_bpp)
#print(ppv_YRNA_fold_bpp)

ppv_YRNA_drtr = ppv(tp_YRNA_drtr, fp_YRNA_drtr)
#print(ppv_YRNA_drtr)

+ppv_YRNA_drtr_mfe = ppv(tp_YRNA_drtr_mfe, fp_YRNA_drtr_mfe)

ppv_YRNA_drtr_prob = ppv(tp_YRNA_drtr_prob, fp_YRNA_drtr_prob)


#F1
f1_YRNA_fold_mfe = f_measure(tp_YRNA_fold_mfe, fp_YRNA_fold_mfe, fn_YRNA_fold_mfe)
#print(f1_YRNA_fold_mfe)

f1_YRNA_fold_bpp = f_measure(tp_YRNA_fold_bpp, fp_YRNA_fold_bpp, fn_YRNA_fold_bpp)
#print(f1_YRNA_fold_bpp)

f1_YRNA_drtr = f_measure(tp_YRNA_drtr, fp_YRNA_drtr, fn_YRNA_drtr)
#print(f1_YRNA_drtr)

f1_YRNA_drtr_mfe = f_measure(tp_YRNA_drtr_mfe, fp_YRNA_drtr_mfe, fn_YRNA_drtr_mfe)

f1_YRNA_drtr_prob = f_measure(tp_YRNA_drtr_prob, fp_YRNA_drtr_prob, fn_YRNA_drtr_prob)


#MCC
mcc_YRNA_fold_mfe = mcc (tp_YRNA_fold_mfe, fp_YRNA_fold_mfe, tn_YRNA_fold_mfe, fn_YRNA_fold_mfe)
#print(mcc_YRNA_fold_mfe)

mcc_YRNA_fold_bpp = mcc (tp_YRNA_fold_bpp, fp_YRNA_fold_bpp, tn_YRNA_fold_bpp, fn_YRNA_fold_bpp)
#print(mcc_YRNA_fold_bpp)

mcc_YRNA_drtr = mcc (tp_YRNA_drtr, fp_YRNA_drtr, tn_YRNA_drtr, fn_YRNA_drtr)
#print(mcc_YRNA_drtr)

mcc_YRNA_drtr_mfe = mcc (tp_YRNA_drtr_mfe, fp_YRNA_drtr_mfe, tn_YRNA_drtr_mfe, fn_YRNA_drtr_mfe)

mcc_YRNA_drtr_prob = mcc (tp_YRNA_drtr_prob, fp_YRNA_drtr_prob, tn_YRNA_drtr_prob, fn_YRNA_drtr_prob)



if args.verbose == True:
    print(raw_data)

if args.graphics == True:

    fold_mcc_mfes = [
    mcc_sixteen_SrRNA_fold_mfe,                 #0
    mcc_GIIIntron_fold_mfe,                     #5
    mcc_RNaseMRPRNA_fold_mfe,                   #15
    mcc_hdvr_ribozyme_fold_mfe,                 #8
    mcc_tmRNA_fold_mfe,                         #20
    mcc_GIIntron_fold_mfe,                      #6
    mcc_other_ribozyme_fold_mfe,                #10
    mcc_cili_telo_RNA_fold_mfe,                 #3
    mcc_RNasePRNA_fold_mfe,                     #16
    mcc_Vert_Telo_RNA_fold_mfe,                 #22
    mcc_ham_ribozyme_fold_mfe,                  #7
    mcc_SRPRNA_fold_mfe,                        #18
    mcc_YRNA_fold_mfe,                          #24
    mcc_cis_regulatory_element_fold_mfe,        #4
    mcc_tRNA_fold_mfe,                          #21
    mcc_five_SrRNA_fold_mfe,                    #2
    mcc_other_RNA_fold_mfe,                     #11
    mcc_RNAIII_fold_mfe,                        #13
    mcc_RNaseE5UTR_fold_mfe,                    #14
    mcc_other_rRNA_fold_mfe,                    #12
    mcc_SyntheticRNA_fold_mfe,                  #19
    mcc_ires_fold_mfe,                          #9
    mcc_twentythree_SrRNA_fold_mfe,             #1
    mcc_Viral_Phage_fold_mfe,                   #23
    mcc_snRNA_fold_mfe                          #17
    ]

    fold_mcc_bpps = (
    mcc_sixteen_SrRNA_fold_bpp,
    mcc_GIIIntron_fold_bpp,
    mcc_RNaseMRPRNA_fold_bpp,
    mcc_hdvr_ribozyme_fold_bpp,
    mcc_tmRNA_fold_bpp,
    mcc_GIIntron_fold_bpp,
    mcc_other_ribozyme_fold_bpp,
    mcc_cili_telo_RNA_fold_bpp,
    mcc_RNasePRNA_fold_bpp,
    mcc_Vert_Telo_RNA_fold_bpp,
    mcc_ham_ribozyme_fold_bpp,
    mcc_SRPRNA_fold_bpp,
    mcc_YRNA_fold_bpp,
    mcc_cis_regulatory_element_fold_bpp,
    mcc_tRNA_fold_bpp,
    mcc_five_SrRNA_fold_bpp,
    mcc_other_RNA_fold_bpp,
    mcc_RNAIII_fold_bpp,
    mcc_RNaseE5UTR_fold_bpp,
    mcc_other_rRNA_fold_bpp,
    mcc_SyntheticRNA_fold_bpp,
    mcc_ires_fold_bpp,
    mcc_twentythree_SrRNA_fold_bpp,
    mcc_Viral_Phage_fold_bpp,
    mcc_snRNA_fold_bpp
    )

    fold_mcc_drtrs = (
    mcc_sixteen_SrRNA_drtr,
    mcc_GIIIntron_drtr,
    mcc_RNaseMRPRNA_drtr,
    mcc_hdvr_ribozyme_drtr,
    mcc_tmRNA_drtr,
    mcc_GIIntron_drtr,
    mcc_other_ribozyme_drtr,
    mcc_cili_telo_RNA_drtr,
    mcc_RNasePRNA_drtr,
    mcc_Vert_Telo_RNA_drtr,
    mcc_ham_ribozyme_drtr,
    mcc_SRPRNA_drtr,
    mcc_YRNA_drtr,
    mcc_cis_regulatory_element_drtr,
    mcc_tRNA_drtr,
    mcc_five_SrRNA_drtr,
    mcc_other_RNA_drtr,
    mcc_RNAIII_drtr,
    mcc_RNaseE5UTR_drtr,
    mcc_other_rRNA_drtr,
    mcc_SyntheticRNA_drtr,
    mcc_ires_drtr,
    mcc_twentythree_SrRNA_drtr,
    mcc_Viral_Phage_drtr,
    mcc_snRNA_drtr
    )

    fold_mcc_drtr_mfes = (
    mcc_sixteen_SrRNA_drtr_mfe,
    mcc_GIIIntron_drtr_mfe,
    mcc_RNaseMRPRNA_drtr_mfe,
    mcc_hdvr_ribozyme_drtr_mfe,
    mcc_tmRNA_drtr_mfe,
    mcc_GIIntron_drtr_mfe,
    mcc_other_ribozyme_drtr_mfe,
    mcc_cili_telo_RNA_drtr_mfe,
    mcc_RNasePRNA_drtr_mfe,
    mcc_Vert_Telo_RNA_drtr_mfe,
    mcc_ham_ribozyme_drtr_mfe,
    mcc_SRPRNA_drtr_mfe,
    mcc_YRNA_drtr_mfe,
    mcc_cis_regulatory_element_drtr_mfe,
    mcc_tRNA_drtr_mfe,
    mcc_five_SrRNA_drtr_mfe,
    mcc_other_RNA_drtr_mfe,
    mcc_RNAIII_drtr_mfe,
    mcc_RNaseE5UTR_drtr_mfe,
    mcc_other_rRNA_drtr_mfe,
    mcc_SyntheticRNA_drtr_mfe,
    mcc_ires_drtr_mfe,
    mcc_twentythree_SrRNA_drtr_mfe,
    mcc_Viral_Phage_drtr_mfe,
    mcc_snRNA_drtr_mfe
    )

    fold_mcc_drtr_probs = (
    mcc_sixteen_SrRNA_drtr_prob,
    mcc_GIIIntron_drtr_prob,
    mcc_RNaseMRPRNA_drtr_prob,
    mcc_hdvr_ribozyme_drtr_prob,
    mcc_tmRNA_drtr_prob,
    mcc_GIIntron_drtr_prob,
    mcc_other_ribozyme_drtr_prob,
    mcc_cili_telo_RNA_drtr_prob,
    mcc_RNasePRNA_drtr_prob,
    mcc_Vert_Telo_RNA_drtr_prob,
    mcc_ham_ribozyme_drtr_prob,
    mcc_SRPRNA_drtr_prob,
    mcc_YRNA_drtr_prob,
    mcc_cis_regulatory_element_drtr_prob,
    mcc_tRNA_drtr_prob,
    mcc_five_SrRNA_drtr_prob,
    mcc_other_RNA_drtr_prob,
    mcc_RNAIII_drtr_prob,
    mcc_RNaseE5UTR_drtr_prob,
    mcc_other_rRNA_drtr_prob,
    mcc_SyntheticRNA_drtr_prob,
    mcc_ires_drtr_prob,
    mcc_twentythree_SrRNA_drtr_prob,
    mcc_Viral_Phage_drtr_prob,
    mcc_snRNA_drtr_prob
    )

    fold_ppv_mfes = [
    ppv_sixteen_SrRNA_fold_mfe,                 #0
    ppv_GIIIntron_fold_mfe,                     #5
    ppv_RNaseMRPRNA_fold_mfe,                   #15
    ppv_hdvr_ribozyme_fold_mfe,                 #8
    ppv_tmRNA_fold_mfe,                         #20
    ppv_GIIntron_fold_mfe,                      #6
    ppv_other_ribozyme_fold_mfe,                #10
    ppv_cili_telo_RNA_fold_mfe,                 #3
    ppv_RNasePRNA_fold_mfe,                     #16
    ppv_Vert_Telo_RNA_fold_mfe,                 #22
    ppv_ham_ribozyme_fold_mfe,                  #7
    ppv_SRPRNA_fold_mfe,                        #18
    ppv_YRNA_fold_mfe,                          #24
    ppv_cis_regulatory_element_fold_mfe,        #4
    ppv_tRNA_fold_mfe,                          #21
    ppv_five_SrRNA_fold_mfe,                    #2
    ppv_other_RNA_fold_mfe,                     #11
    ppv_RNAIII_fold_mfe,                        #13
    ppv_RNaseE5UTR_fold_mfe,                    #14
    ppv_other_rRNA_fold_mfe,                    #12
    ppv_SyntheticRNA_fold_mfe,                  #19
    ppv_ires_fold_mfe,                          #9
    ppv_twentythree_SrRNA_fold_mfe,             #1
    ppv_Viral_Phage_fold_mfe,                   #23
    ppv_snRNA_fold_mfe                          #17
    ]

    fold_ppv_bpps = (
    ppv_sixteen_SrRNA_fold_bpp,
    ppv_GIIIntron_fold_bpp,
    ppv_RNaseMRPRNA_fold_bpp,
    ppv_hdvr_ribozyme_fold_bpp,
    ppv_tmRNA_fold_bpp,
    ppv_GIIntron_fold_bpp,
    ppv_other_ribozyme_fold_bpp,
    ppv_cili_telo_RNA_fold_bpp,
    ppv_RNasePRNA_fold_bpp,
    ppv_Vert_Telo_RNA_fold_bpp,
    ppv_ham_ribozyme_fold_bpp,
    ppv_SRPRNA_fold_bpp,
    ppv_YRNA_fold_bpp,
    ppv_cis_regulatory_element_fold_bpp,
    ppv_tRNA_fold_bpp,
    ppv_five_SrRNA_fold_bpp,
    ppv_other_RNA_fold_bpp,
    ppv_RNAIII_fold_bpp,
    ppv_RNaseE5UTR_fold_bpp,
    ppv_other_rRNA_fold_bpp,
    ppv_SyntheticRNA_fold_bpp,
    ppv_ires_fold_bpp,
    ppv_twentythree_SrRNA_fold_bpp,
    ppv_Viral_Phage_fold_bpp,
    ppv_snRNA_fold_bpp
    )

    fold_ppv_drtrs = (
    ppv_sixteen_SrRNA_drtr,
    ppv_GIIIntron_drtr,
    ppv_RNaseMRPRNA_drtr,
    ppv_hdvr_ribozyme_drtr,
    ppv_tmRNA_drtr,
    ppv_GIIntron_drtr,
    ppv_other_ribozyme_drtr,
    ppv_cili_telo_RNA_drtr,
    ppv_RNasePRNA_drtr,
    ppv_Vert_Telo_RNA_drtr,
    ppv_ham_ribozyme_drtr,
    ppv_SRPRNA_drtr,
    ppv_YRNA_drtr,
    ppv_cis_regulatory_element_drtr,
    ppv_tRNA_drtr,
    ppv_five_SrRNA_drtr,
    ppv_other_RNA_drtr,
    ppv_RNAIII_drtr,
    ppv_RNaseE5UTR_drtr,
    ppv_other_rRNA_drtr,
    ppv_SyntheticRNA_drtr,
    ppv_ires_drtr,
    ppv_twentythree_SrRNA_drtr,
    ppv_Viral_Phage_drtr,
    ppv_snRNA_drtr
    )

    fold_ppv_drtr_mfes = (
    ppv_sixteen_SrRNA_drtr_mfe,
    ppv_GIIIntron_drtr_mfe,
    ppv_RNaseMRPRNA_drtr_mfe,
    ppv_hdvr_ribozyme_drtr_mfe,
    ppv_tmRNA_drtr_mfe,
    ppv_GIIntron_drtr_mfe,
    ppv_other_ribozyme_drtr_mfe,
    ppv_cili_telo_RNA_drtr_mfe,
    ppv_RNasePRNA_drtr_mfe,
    ppv_Vert_Telo_RNA_drtr_mfe,
    ppv_ham_ribozyme_drtr_mfe,
    ppv_SRPRNA_drtr_mfe,
    ppv_YRNA_drtr_mfe,
    ppv_cis_regulatory_element_drtr_mfe,
    ppv_tRNA_drtr_mfe,
    ppv_five_SrRNA_drtr_mfe,
    ppv_other_RNA_drtr_mfe,
    ppv_RNAIII_drtr_mfe,
    ppv_RNaseE5UTR_drtr_mfe,
    ppv_other_rRNA_drtr_mfe,
    ppv_SyntheticRNA_drtr_mfe,
    ppv_ires_drtr_mfe,
    ppv_twentythree_SrRNA_drtr_mfe,
    ppv_Viral_Phage_drtr_mfe,
    ppv_snRNA_drtr_mfe
    )

    fold_ppv_drtr_probs = (
    ppv_sixteen_SrRNA_drtr_prob,
    ppv_GIIIntron_drtr_prob,
    ppv_RNaseMRPRNA_drtr_prob,
    ppv_hdvr_ribozyme_drtr_prob,
    ppv_tmRNA_drtr_prob,
    ppv_GIIntron_drtr_prob,
    ppv_other_ribozyme_drtr_prob,
    ppv_cili_telo_RNA_drtr_prob,
    ppv_RNasePRNA_drtr_prob,
    ppv_Vert_Telo_RNA_drtr_prob,
    ppv_ham_ribozyme_drtr_prob,
    ppv_SRPRNA_drtr_prob,
    ppv_YRNA_drtr_prob,
    ppv_cis_regulatory_element_drtr_prob,
    ppv_tRNA_drtr_prob,
    ppv_five_SrRNA_drtr_prob,
    ppv_other_RNA_drtr_prob,
    ppv_RNAIII_drtr_prob,
    ppv_RNaseE5UTR_drtr_prob,
    ppv_other_rRNA_drtr_prob,
    ppv_SyntheticRNA_drtr_prob,
    ppv_ires_drtr_prob,
    ppv_twentythree_SrRNA_drtr_prob,
    ppv_Viral_Phage_drtr_prob,
    ppv_snRNA_drtr_prob
    )




    fold_tpr_mfes = [
    tpr_sixteen_SrRNA_fold_mfe,                 #0
    tpr_GIIIntron_fold_mfe,                     #5
    tpr_RNaseMRPRNA_fold_mfe,                   #15
    tpr_hdvr_ribozyme_fold_mfe,                 #8
    tpr_tmRNA_fold_mfe,                         #20
    tpr_GIIntron_fold_mfe,                      #6
    tpr_other_ribozyme_fold_mfe,                #10
    tpr_cili_telo_RNA_fold_mfe,                 #3
    tpr_RNasePRNA_fold_mfe,                     #16
    tpr_Vert_Telo_RNA_fold_mfe,                 #22
    tpr_ham_ribozyme_fold_mfe,                  #7
    tpr_SRPRNA_fold_mfe,                        #18
    tpr_YRNA_fold_mfe,                          #24
    tpr_cis_regulatory_element_fold_mfe,        #4
    tpr_tRNA_fold_mfe,                          #21
    tpr_five_SrRNA_fold_mfe,                    #2
    tpr_other_RNA_fold_mfe,                     #11
    tpr_RNAIII_fold_mfe,                        #13
    tpr_RNaseE5UTR_fold_mfe,                    #14
    tpr_other_rRNA_fold_mfe,                    #12
    tpr_SyntheticRNA_fold_mfe,                  #19
    tpr_ires_fold_mfe,                          #9
    tpr_twentythree_SrRNA_fold_mfe,             #1
    tpr_Viral_Phage_fold_mfe,                   #23
    tpr_snRNA_fold_mfe                          #17
    ]

    fold_tpr_bpps = (
    tpr_sixteen_SrRNA_fold_bpp,
    tpr_GIIIntron_fold_bpp,
    tpr_RNaseMRPRNA_fold_bpp,
    tpr_hdvr_ribozyme_fold_bpp,
    tpr_tmRNA_fold_bpp,
    tpr_GIIntron_fold_bpp,
    tpr_other_ribozyme_fold_bpp,
    tpr_cili_telo_RNA_fold_bpp,
    tpr_RNasePRNA_fold_bpp,
    tpr_Vert_Telo_RNA_fold_bpp,
    tpr_ham_ribozyme_fold_bpp,
    tpr_SRPRNA_fold_bpp,
    tpr_YRNA_fold_bpp,
    tpr_cis_regulatory_element_fold_bpp,
    tpr_tRNA_fold_bpp,
    tpr_five_SrRNA_fold_bpp,
    tpr_other_RNA_fold_bpp,
    tpr_RNAIII_fold_bpp,
    tpr_RNaseE5UTR_fold_bpp,
    tpr_other_rRNA_fold_bpp,
    tpr_SyntheticRNA_fold_bpp,
    tpr_ires_fold_bpp,
    tpr_twentythree_SrRNA_fold_bpp,
    tpr_Viral_Phage_fold_bpp,
    tpr_snRNA_fold_bpp
    )

    fold_tpr_drtrs = (
    tpr_sixteen_SrRNA_drtr,
    tpr_GIIIntron_drtr,
    tpr_RNaseMRPRNA_drtr,
    tpr_hdvr_ribozyme_drtr,
    tpr_tmRNA_drtr,
    tpr_GIIntron_drtr,
    tpr_other_ribozyme_drtr,
    tpr_cili_telo_RNA_drtr,
    tpr_RNasePRNA_drtr,
    tpr_Vert_Telo_RNA_drtr,
    tpr_ham_ribozyme_drtr,
    tpr_SRPRNA_drtr,
    tpr_YRNA_drtr,
    tpr_cis_regulatory_element_drtr,
    tpr_tRNA_drtr,
    tpr_five_SrRNA_drtr,
    tpr_other_RNA_drtr,
    tpr_RNAIII_drtr,
    tpr_RNaseE5UTR_drtr,
    tpr_other_rRNA_drtr,
    tpr_SyntheticRNA_drtr,
    tpr_ires_drtr,
    tpr_twentythree_SrRNA_drtr,
    tpr_Viral_Phage_drtr,
    tpr_snRNA_drtr
    )

    fold_tpr_drtr_mfes = (
    tpr_sixteen_SrRNA_drtr_mfe,
    tpr_GIIIntron_drtr_mfe,
    tpr_RNaseMRPRNA_drtr_mfe,
    tpr_hdvr_ribozyme_drtr_mfe,
    tpr_tmRNA_drtr_mfe,
    tpr_GIIntron_drtr_mfe,
    tpr_other_ribozyme_drtr_mfe,
    tpr_cili_telo_RNA_drtr_mfe,
    tpr_RNasePRNA_drtr_mfe,
    tpr_Vert_Telo_RNA_drtr_mfe,
    tpr_ham_ribozyme_drtr_mfe,
    tpr_SRPRNA_drtr_mfe,
    tpr_YRNA_drtr_mfe,
    tpr_cis_regulatory_element_drtr_mfe,
    tpr_tRNA_drtr_mfe,
    tpr_five_SrRNA_drtr_mfe,
    tpr_other_RNA_drtr_mfe,
    tpr_RNAIII_drtr_mfe,
    tpr_RNaseE5UTR_drtr_mfe,
    tpr_other_rRNA_drtr_mfe,
    tpr_SyntheticRNA_drtr_mfe,
    tpr_ires_drtr_mfe,
    tpr_twentythree_SrRNA_drtr_mfe,
    tpr_Viral_Phage_drtr_mfe,
    tpr_snRNA_drtr_mfe
    )

    fold_tpr_drtr_probs = (
    tpr_sixteen_SrRNA_drtr_prob,
    tpr_GIIIntron_drtr_prob,
    tpr_RNaseMRPRNA_drtr_prob,
    tpr_hdvr_ribozyme_drtr_prob,
    tpr_tmRNA_drtr_prob,
    tpr_GIIntron_drtr_prob,
    tpr_other_ribozyme_drtr_prob,
    tpr_cili_telo_RNA_drtr_prob,
    tpr_RNasePRNA_drtr_prob,
    tpr_Vert_Telo_RNA_drtr_prob,
    tpr_ham_ribozyme_drtr_prob,
    tpr_SRPRNA_drtr_prob,
    tpr_YRNA_drtr_prob,
    tpr_cis_regulatory_element_drtr_prob,
    tpr_tRNA_drtr_prob,
    tpr_five_SrRNA_drtr_prob,
    tpr_other_RNA_drtr_prob,
    tpr_RNAIII_drtr_prob,
    tpr_RNaseE5UTR_drtr_prob,
    tpr_other_rRNA_drtr_prob,
    tpr_SyntheticRNA_drtr_prob,
    tpr_ires_drtr_prob,
    tpr_twentythree_SrRNA_drtr_prob,
    tpr_Viral_Phage_drtr_prob,
    tpr_snRNA_drtr_prob
    )




    fold_f1_mfes = [
    f1_sixteen_SrRNA_fold_mfe,                 #0
    f1_GIIIntron_fold_mfe,                     #5
    f1_RNaseMRPRNA_fold_mfe,                   #15
    f1_hdvr_ribozyme_fold_mfe,                 #8
    f1_tmRNA_fold_mfe,                         #20
    f1_GIIntron_fold_mfe,                      #6
    f1_other_ribozyme_fold_mfe,                #10
    f1_cili_telo_RNA_fold_mfe,                 #3
    f1_RNasePRNA_fold_mfe,                     #16
    f1_Vert_Telo_RNA_fold_mfe,                 #22
    f1_ham_ribozyme_fold_mfe,                  #7
    f1_SRPRNA_fold_mfe,                        #18
    f1_YRNA_fold_mfe,                          #24
    f1_cis_regulatory_element_fold_mfe,        #4
    f1_tRNA_fold_mfe,                          #21
    f1_five_SrRNA_fold_mfe,                    #2
    f1_other_RNA_fold_mfe,                     #11
    f1_RNAIII_fold_mfe,                        #13
    f1_RNaseE5UTR_fold_mfe,                    #14
    f1_other_rRNA_fold_mfe,                    #12
    f1_SyntheticRNA_fold_mfe,                  #19
    f1_ires_fold_mfe,                          #9
    f1_twentythree_SrRNA_fold_mfe,             #1
    f1_Viral_Phage_fold_mfe,                   #23
    f1_snRNA_fold_mfe                          #17
    ]

    fold_f1_bpps = (
    f1_sixteen_SrRNA_fold_bpp,
    f1_GIIIntron_fold_bpp,
    f1_RNaseMRPRNA_fold_bpp,
    f1_hdvr_ribozyme_fold_bpp,
    f1_tmRNA_fold_bpp,
    f1_GIIntron_fold_bpp,
    f1_other_ribozyme_fold_bpp,
    f1_cili_telo_RNA_fold_bpp,
    f1_RNasePRNA_fold_bpp,
    f1_Vert_Telo_RNA_fold_bpp,
    f1_ham_ribozyme_fold_bpp,
    f1_SRPRNA_fold_bpp,
    f1_YRNA_fold_bpp,
    f1_cis_regulatory_element_fold_bpp,
    f1_tRNA_fold_bpp,
    f1_five_SrRNA_fold_bpp,
    f1_other_RNA_fold_bpp,
    f1_RNAIII_fold_bpp,
    f1_RNaseE5UTR_fold_bpp,
    f1_other_rRNA_fold_bpp,
    f1_SyntheticRNA_fold_bpp,
    f1_ires_fold_bpp,
    f1_twentythree_SrRNA_fold_bpp,
    f1_Viral_Phage_fold_bpp,
    f1_snRNA_fold_bpp
    )

    fold_f1_drtrs = (
    f1_sixteen_SrRNA_drtr,
    f1_GIIIntron_drtr,
    f1_RNaseMRPRNA_drtr,
    f1_hdvr_ribozyme_drtr,
    f1_tmRNA_drtr,
    f1_GIIntron_drtr,
    f1_other_ribozyme_drtr,
    f1_cili_telo_RNA_drtr,
    f1_RNasePRNA_drtr,
    f1_Vert_Telo_RNA_drtr,
    f1_ham_ribozyme_drtr,
    f1_SRPRNA_drtr,
    f1_YRNA_drtr,
    f1_cis_regulatory_element_drtr,
    f1_tRNA_drtr,
    f1_five_SrRNA_drtr,
    f1_other_RNA_drtr,
    f1_RNAIII_drtr,
    f1_RNaseE5UTR_drtr,
    f1_other_rRNA_drtr,
    f1_SyntheticRNA_drtr,
    f1_ires_drtr,
    f1_twentythree_SrRNA_drtr,
    f1_Viral_Phage_drtr,
    f1_snRNA_drtr
    )

    fold_f1_drtr_mfes = (
    f1_sixteen_SrRNA_drtr_mfe,
    f1_GIIIntron_drtr_mfe,
    f1_RNaseMRPRNA_drtr_mfe,
    f1_hdvr_ribozyme_drtr_mfe,
    f1_tmRNA_drtr_mfe,
    f1_GIIntron_drtr_mfe,
    f1_other_ribozyme_drtr_mfe,
    f1_cili_telo_RNA_drtr_mfe,
    f1_RNasePRNA_drtr_mfe,
    f1_Vert_Telo_RNA_drtr_mfe,
    f1_ham_ribozyme_drtr_mfe,
    f1_SRPRNA_drtr_mfe,
    f1_YRNA_drtr_mfe,
    f1_cis_regulatory_element_drtr_mfe,
    f1_tRNA_drtr_mfe,
    f1_five_SrRNA_drtr_mfe,
    f1_other_RNA_drtr_mfe,
    f1_RNAIII_drtr_mfe,
    f1_RNaseE5UTR_drtr_mfe,
    f1_other_rRNA_drtr_mfe,
    f1_SyntheticRNA_drtr_mfe,
    f1_ires_drtr_mfe,
    f1_twentythree_SrRNA_drtr_mfe,
    f1_Viral_Phage_drtr_mfe,
    f1_snRNA_drtr_mfe
    )

    fold_f1_drtr_probs = (
    f1_sixteen_SrRNA_drtr_prob,
    f1_GIIIntron_drtr_prob,
    f1_RNaseMRPRNA_drtr_prob,
    f1_hdvr_ribozyme_drtr_prob,
    f1_tmRNA_drtr_prob,
    f1_GIIntron_drtr_prob,
    f1_other_ribozyme_drtr_prob,
    f1_cili_telo_RNA_drtr_prob,
    f1_RNasePRNA_drtr_prob,
    f1_Vert_Telo_RNA_drtr_prob,
    f1_ham_ribozyme_drtr_prob,
    f1_SRPRNA_drtr_prob,
    f1_YRNA_drtr_prob,
    f1_cis_regulatory_element_drtr_prob,
    f1_tRNA_drtr_prob,
    f1_five_SrRNA_drtr_prob,
    f1_other_RNA_drtr_prob,
    f1_RNAIII_drtr_prob,
    f1_RNaseE5UTR_drtr_prob,
    f1_other_rRNA_drtr_prob,
    f1_SyntheticRNA_drtr_prob,
    f1_ires_drtr_prob,
    f1_twentythree_SrRNA_drtr_prob,
    f1_Viral_Phage_drtr_prob,
    f1_snRNA_drtr_prob
    )




    titel = args.label
    barwidth = 0.19

    r_1 = np.arange(len(fold_mcc_mfes))
    r_2 = [x + barwidth for x in r_1]
    r_3 = [x + barwidth for x in r_2]
    r_4 = [x + barwidth for x in r_3]
    r_5 = [x + barwidth for x in r_4]


    plt.bar(r_1, fold_mcc_mfes, width=barwidth, color='#247E85', label='RNAfold: mfe')
    plt.bar(r_2, fold_mcc_drtr_mfes, width=barwidth, color='#87f542', label='DrTransformer: mfe')
    plt.bar(r_3, fold_mcc_bpps, width=barwidth, color='#FA5882', label='RNAfold: bpp')
    plt.bar(r_4, fold_mcc_drtrs, width=barwidth, color='#0E3151', label='DrTransformer: bpp')
    plt.bar(r_5, fold_mcc_drtr_probs, width=barwidth, color='#583151', label='DrTransformer: prob')


    plt.title(titel, fontsize=20, y=1.03)
    plt.xlabel('Family', fontsize=16)
    plt.ylabel('MCC',  fontsize=16)
    plt.axis([-0.4, 25, 0, 1])
    plt.xticks([r + barwidth for r in range (len(fold_mcc_mfes))],
                                                   [family_dict[0] + "(n=" + str(len(tp_sixteen_SrRNA_drtr)) + ")",
                                                   family_dict[5] + "(n=" + str(len(tp_GIIIntron_drtr)) + ")",
                                                   family_dict[15] + "(n=" + str(len(tp_RNaseMRPRNA_drtr)) + ")",
                                                   family_dict[8] + "(n=" + str(len(tp_hdvr_ribozyme_drtr)) + ")",
                                                   family_dict[20] + "(n=" + str(len(tp_tmRNA_drtr)) + ")",
                                                   family_dict[6] + "(n=" + str(len(tp_GIIntron_drtr)) + ")",
                                                   family_dict[10] + "(n=" + str(len(tp_other_ribozyme_drtr)) + ")",
                                                   family_dict[3] + "(n=" + str(len(tp_cili_telo_RNA_drtr)) + ")",
                                                   family_dict[16] + "(n=" + str(len(tp_RNasePRNA_drtr)) + ")",
                                                   family_dict[22] + "(n=" + str(len(tp_Vert_Telo_RNA_drtr)) + ")",
                                                   family_dict[7] + "(n=" + str(len(tp_ham_ribozyme_drtr)) + ")",
                                                   family_dict[18] + "(n=" + str(len(tp_SRPRNA_drtr)) + ")",
                                                   family_dict[24] + "(n=" + str(len(tp_YRNA_drtr)) + ")",
                                                   family_dict[4] + "(n=" + str(len(tp_cis_regulatory_element_drtr)) + ")",
                                                   family_dict[21] + "(n=" + str(len(tp_tRNA_drtr)) + ")",
                                                   family_dict[2] + "(n=" + str(len(tp_five_SrRNA_drtr)) + ")",
                                                   family_dict[11] + "(n=" + str(len(tp_other_RNA_drtr)) + ")",
                                                   family_dict[13] + "(n=" + str(len(tp_RNAIII_drtr)) + ")",
                                                   family_dict[14] + "(n=" + str(len(tp_RNaseE5UTR_drtr)) + ")",
                                                   family_dict[12] + "(n=" + str(len(tp_other_rRNA_drtr)) + ")",
                                                   family_dict[19] + "(n=" + str(len(tp_SyntheticRNA_drtr)) + ")",
                                                   family_dict[9] + "(n=" + str(len(tp_ires_drtr)) + ")",
                                                   family_dict[1] + "(n=" + str(len(tp_twentythree_SrRNA_drtr)) + ")",
                                                   family_dict[23] + "(n=" + str(len(tp_Viral_Phage_drtr)) + ")",
                                                   family_dict[17] + "(n=" + str(len(tp_snRNA_drtr)) + ")"],
                                                   rotation=56,
                                                   ha="right",
                                                   rotation_mode="anchor",
                                                   fontsize=12)

    plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    plt.gcf().subplots_adjust(bottom=0.31, left=0.06, right=0.94, top=0.91)
    plt.legend()
    plt.show()

if args.output==True:

    fold_mcc_mfes = [
    mcc_sixteen_SrRNA_fold_mfe,                 #0
    mcc_GIIIntron_fold_mfe,                     #5
    mcc_RNaseMRPRNA_fold_mfe,                   #15
    mcc_hdvr_ribozyme_fold_mfe,                 #8
    mcc_tmRNA_fold_mfe,                         #20
    mcc_GIIntron_fold_mfe,                      #6
    mcc_other_ribozyme_fold_mfe,                #10
    mcc_cili_telo_RNA_fold_mfe,                 #3
    mcc_RNasePRNA_fold_mfe,                     #16
    mcc_Vert_Telo_RNA_fold_mfe,                 #22
    mcc_ham_ribozyme_fold_mfe,                  #7
    mcc_SRPRNA_fold_mfe,                        #18
    mcc_YRNA_fold_mfe,                          #24
    mcc_cis_regulatory_element_fold_mfe,        #4
    mcc_tRNA_fold_mfe,                          #21
    mcc_five_SrRNA_fold_mfe,                    #2
    mcc_other_RNA_fold_mfe,                     #11
    mcc_RNAIII_fold_mfe,                        #13
    mcc_RNaseE5UTR_fold_mfe,                    #14
    mcc_other_rRNA_fold_mfe,                    #12
    mcc_SyntheticRNA_fold_mfe,                  #19
    mcc_ires_fold_mfe,                          #9
    mcc_twentythree_SrRNA_fold_mfe,             #1
    mcc_Viral_Phage_fold_mfe,                   #23
    mcc_snRNA_fold_mfe                          #17
    ]

    fold_mcc_bpps = (
    mcc_sixteen_SrRNA_fold_bpp,
    mcc_GIIIntron_fold_bpp,
    mcc_RNaseMRPRNA_fold_bpp,
    mcc_hdvr_ribozyme_fold_bpp,
    mcc_tmRNA_fold_bpp,
    mcc_GIIntron_fold_bpp,
    mcc_other_ribozyme_fold_bpp,
    mcc_cili_telo_RNA_fold_bpp,
    mcc_RNasePRNA_fold_bpp,
    mcc_Vert_Telo_RNA_fold_bpp,
    mcc_ham_ribozyme_fold_bpp,
    mcc_SRPRNA_fold_bpp,
    mcc_YRNA_fold_bpp,
    mcc_cis_regulatory_element_fold_bpp,
    mcc_tRNA_fold_bpp,
    mcc_five_SrRNA_fold_bpp,
    mcc_other_RNA_fold_bpp,
    mcc_RNAIII_fold_bpp,
    mcc_RNaseE5UTR_fold_bpp,
    mcc_other_rRNA_fold_bpp,
    mcc_SyntheticRNA_fold_bpp,
    mcc_ires_fold_bpp,
    mcc_twentythree_SrRNA_fold_bpp,
    mcc_Viral_Phage_fold_bpp,
    mcc_snRNA_fold_bpp
    )

    fold_mcc_drtrs = (
    mcc_sixteen_SrRNA_drtr,
    mcc_GIIIntron_drtr,
    mcc_RNaseMRPRNA_drtr,
    mcc_hdvr_ribozyme_drtr,
    mcc_tmRNA_drtr,
    mcc_GIIntron_drtr,
    mcc_other_ribozyme_drtr,
    mcc_cili_telo_RNA_drtr,
    mcc_RNasePRNA_drtr,
    mcc_Vert_Telo_RNA_drtr,
    mcc_ham_ribozyme_drtr,
    mcc_SRPRNA_drtr,
    mcc_YRNA_drtr,
    mcc_cis_regulatory_element_drtr,
    mcc_tRNA_drtr,
    mcc_five_SrRNA_drtr,
    mcc_other_RNA_drtr,
    mcc_RNAIII_drtr,
    mcc_RNaseE5UTR_drtr,
    mcc_other_rRNA_drtr,
    mcc_SyntheticRNA_drtr,
    mcc_ires_drtr,
    mcc_twentythree_SrRNA_drtr,
    mcc_Viral_Phage_drtr,
    mcc_snRNA_drtr
    )

    fold_mcc_drtr_mfes = (
    mcc_sixteen_SrRNA_drtr_mfe,
    mcc_GIIIntron_drtr_mfe,
    mcc_RNaseMRPRNA_drtr_mfe,
    mcc_hdvr_ribozyme_drtr_mfe,
    mcc_tmRNA_drtr_mfe,
    mcc_GIIntron_drtr_mfe,
    mcc_other_ribozyme_drtr_mfe,
    mcc_cili_telo_RNA_drtr_mfe,
    mcc_RNasePRNA_drtr_mfe,
    mcc_Vert_Telo_RNA_drtr_mfe,
    mcc_ham_ribozyme_drtr_mfe,
    mcc_SRPRNA_drtr_mfe,
    mcc_YRNA_drtr_mfe,
    mcc_cis_regulatory_element_drtr_mfe,
    mcc_tRNA_drtr_mfe,
    mcc_five_SrRNA_drtr_mfe,
    mcc_other_RNA_drtr_mfe,
    mcc_RNAIII_drtr_mfe,
    mcc_RNaseE5UTR_drtr_mfe,
    mcc_other_rRNA_drtr_mfe,
    mcc_SyntheticRNA_drtr_mfe,
    mcc_ires_drtr_mfe,
    mcc_twentythree_SrRNA_drtr_mfe,
    mcc_Viral_Phage_drtr_mfe,
    mcc_snRNA_drtr_mfe
    )

    fold_mcc_drtr_probs = (
    mcc_sixteen_SrRNA_drtr_prob,
    mcc_GIIIntron_drtr_prob,
    mcc_RNaseMRPRNA_drtr_prob,
    mcc_hdvr_ribozyme_drtr_prob,
    mcc_tmRNA_drtr_prob,
    mcc_GIIntron_drtr_prob,
    mcc_other_ribozyme_drtr_prob,
    mcc_cili_telo_RNA_drtr_prob,
    mcc_RNasePRNA_drtr_prob,
    mcc_Vert_Telo_RNA_drtr_prob,
    mcc_ham_ribozyme_drtr_prob,
    mcc_SRPRNA_drtr_prob,
    mcc_YRNA_drtr_prob,
    mcc_cis_regulatory_element_drtr_prob,
    mcc_tRNA_drtr_prob,
    mcc_five_SrRNA_drtr_prob,
    mcc_other_RNA_drtr_prob,
    mcc_RNAIII_drtr_prob,
    mcc_RNaseE5UTR_drtr_prob,
    mcc_other_rRNA_drtr_prob,
    mcc_SyntheticRNA_drtr_prob,
    mcc_ires_drtr_prob,
    mcc_twentythree_SrRNA_drtr_prob,
    mcc_Viral_Phage_drtr_prob,
    mcc_snRNA_drtr_prob
    )

    fold_ppv_mfes = [
    ppv_sixteen_SrRNA_fold_mfe,                 #0
    ppv_GIIIntron_fold_mfe,                     #5
    ppv_RNaseMRPRNA_fold_mfe,                   #15
    ppv_hdvr_ribozyme_fold_mfe,                 #8
    ppv_tmRNA_fold_mfe,                         #20
    ppv_GIIntron_fold_mfe,                      #6
    ppv_other_ribozyme_fold_mfe,                #10
    ppv_cili_telo_RNA_fold_mfe,                 #3
    ppv_RNasePRNA_fold_mfe,                     #16
    ppv_Vert_Telo_RNA_fold_mfe,                 #22
    ppv_ham_ribozyme_fold_mfe,                  #7
    ppv_SRPRNA_fold_mfe,                        #18
    ppv_YRNA_fold_mfe,                          #24
    ppv_cis_regulatory_element_fold_mfe,        #4
    ppv_tRNA_fold_mfe,                          #21
    ppv_five_SrRNA_fold_mfe,                    #2
    ppv_other_RNA_fold_mfe,                     #11
    ppv_RNAIII_fold_mfe,                        #13
    ppv_RNaseE5UTR_fold_mfe,                    #14
    ppv_other_rRNA_fold_mfe,                    #12
    ppv_SyntheticRNA_fold_mfe,                  #19
    ppv_ires_fold_mfe,                          #9
    ppv_twentythree_SrRNA_fold_mfe,             #1
    ppv_Viral_Phage_fold_mfe,                   #23
    ppv_snRNA_fold_mfe                          #17
    ]

    fold_ppv_bpps = (
    ppv_sixteen_SrRNA_fold_bpp,
    ppv_GIIIntron_fold_bpp,
    ppv_RNaseMRPRNA_fold_bpp,
    ppv_hdvr_ribozyme_fold_bpp,
    ppv_tmRNA_fold_bpp,
    ppv_GIIntron_fold_bpp,
    ppv_other_ribozyme_fold_bpp,
    ppv_cili_telo_RNA_fold_bpp,
    ppv_RNasePRNA_fold_bpp,
    ppv_Vert_Telo_RNA_fold_bpp,
    ppv_ham_ribozyme_fold_bpp,
    ppv_SRPRNA_fold_bpp,
    ppv_YRNA_fold_bpp,
    ppv_cis_regulatory_element_fold_bpp,
    ppv_tRNA_fold_bpp,
    ppv_five_SrRNA_fold_bpp,
    ppv_other_RNA_fold_bpp,
    ppv_RNAIII_fold_bpp,
    ppv_RNaseE5UTR_fold_bpp,
    ppv_other_rRNA_fold_bpp,
    ppv_SyntheticRNA_fold_bpp,
    ppv_ires_fold_bpp,
    ppv_twentythree_SrRNA_fold_bpp,
    ppv_Viral_Phage_fold_bpp,
    ppv_snRNA_fold_bpp
    )

    fold_ppv_drtrs = (
    ppv_sixteen_SrRNA_drtr,
    ppv_GIIIntron_drtr,
    ppv_RNaseMRPRNA_drtr,
    ppv_hdvr_ribozyme_drtr,
    ppv_tmRNA_drtr,
    ppv_GIIntron_drtr,
    ppv_other_ribozyme_drtr,
    ppv_cili_telo_RNA_drtr,
    ppv_RNasePRNA_drtr,
    ppv_Vert_Telo_RNA_drtr,
    ppv_ham_ribozyme_drtr,
    ppv_SRPRNA_drtr,
    ppv_YRNA_drtr,
    ppv_cis_regulatory_element_drtr,
    ppv_tRNA_drtr,
    ppv_five_SrRNA_drtr,
    ppv_other_RNA_drtr,
    ppv_RNAIII_drtr,
    ppv_RNaseE5UTR_drtr,
    ppv_other_rRNA_drtr,
    ppv_SyntheticRNA_drtr,
    ppv_ires_drtr,
    ppv_twentythree_SrRNA_drtr,
    ppv_Viral_Phage_drtr,
    ppv_snRNA_drtr
    )

    fold_ppv_drtr_mfes = (
    ppv_sixteen_SrRNA_drtr_mfe,
    ppv_GIIIntron_drtr_mfe,
    ppv_RNaseMRPRNA_drtr_mfe,
    ppv_hdvr_ribozyme_drtr_mfe,
    ppv_tmRNA_drtr_mfe,
    ppv_GIIntron_drtr_mfe,
    ppv_other_ribozyme_drtr_mfe,
    ppv_cili_telo_RNA_drtr_mfe,
    ppv_RNasePRNA_drtr_mfe,
    ppv_Vert_Telo_RNA_drtr_mfe,
    ppv_ham_ribozyme_drtr_mfe,
    ppv_SRPRNA_drtr_mfe,
    ppv_YRNA_drtr_mfe,
    ppv_cis_regulatory_element_drtr_mfe,
    ppv_tRNA_drtr_mfe,
    ppv_five_SrRNA_drtr_mfe,
    ppv_other_RNA_drtr_mfe,
    ppv_RNAIII_drtr_mfe,
    ppv_RNaseE5UTR_drtr_mfe,
    ppv_other_rRNA_drtr_mfe,
    ppv_SyntheticRNA_drtr_mfe,
    ppv_ires_drtr_mfe,
    ppv_twentythree_SrRNA_drtr_mfe,
    ppv_Viral_Phage_drtr_mfe,
    ppv_snRNA_drtr_mfe
    )

    fold_ppv_drtr_probs = (
    ppv_sixteen_SrRNA_drtr_prob,
    ppv_GIIIntron_drtr_prob,
    ppv_RNaseMRPRNA_drtr_prob,
    ppv_hdvr_ribozyme_drtr_prob,
    ppv_tmRNA_drtr_prob,
    ppv_GIIntron_drtr_prob,
    ppv_other_ribozyme_drtr_prob,
    ppv_cili_telo_RNA_drtr_prob,
    ppv_RNasePRNA_drtr_prob,
    ppv_Vert_Telo_RNA_drtr_prob,
    ppv_ham_ribozyme_drtr_prob,
    ppv_SRPRNA_drtr_prob,
    ppv_YRNA_drtr_prob,
    ppv_cis_regulatory_element_drtr_prob,
    ppv_tRNA_drtr_prob,
    ppv_five_SrRNA_drtr_prob,
    ppv_other_RNA_drtr_prob,
    ppv_RNAIII_drtr_prob,
    ppv_RNaseE5UTR_drtr_prob,
    ppv_other_rRNA_drtr_prob,
    ppv_SyntheticRNA_drtr_prob,
    ppv_ires_drtr_prob,
    ppv_twentythree_SrRNA_drtr_prob,
    ppv_Viral_Phage_drtr_prob,
    ppv_snRNA_drtr_prob
    )




    fold_tpr_mfes = [
    tpr_sixteen_SrRNA_fold_mfe,                 #0
    tpr_GIIIntron_fold_mfe,                     #5
    tpr_RNaseMRPRNA_fold_mfe,                   #15
    tpr_hdvr_ribozyme_fold_mfe,                 #8
    tpr_tmRNA_fold_mfe,                         #20
    tpr_GIIntron_fold_mfe,                      #6
    tpr_other_ribozyme_fold_mfe,                #10
    tpr_cili_telo_RNA_fold_mfe,                 #3
    tpr_RNasePRNA_fold_mfe,                     #16
    tpr_Vert_Telo_RNA_fold_mfe,                 #22
    tpr_ham_ribozyme_fold_mfe,                  #7
    tpr_SRPRNA_fold_mfe,                        #18
    tpr_YRNA_fold_mfe,                          #24
    tpr_cis_regulatory_element_fold_mfe,        #4
    tpr_tRNA_fold_mfe,                          #21
    tpr_five_SrRNA_fold_mfe,                    #2
    tpr_other_RNA_fold_mfe,                     #11
    tpr_RNAIII_fold_mfe,                        #13
    tpr_RNaseE5UTR_fold_mfe,                    #14
    tpr_other_rRNA_fold_mfe,                    #12
    tpr_SyntheticRNA_fold_mfe,                  #19
    tpr_ires_fold_mfe,                          #9
    tpr_twentythree_SrRNA_fold_mfe,             #1
    tpr_Viral_Phage_fold_mfe,                   #23
    tpr_snRNA_fold_mfe                          #17
    ]

    fold_tpr_bpps = (
    tpr_sixteen_SrRNA_fold_bpp,
    tpr_GIIIntron_fold_bpp,
    tpr_RNaseMRPRNA_fold_bpp,
    tpr_hdvr_ribozyme_fold_bpp,
    tpr_tmRNA_fold_bpp,
    tpr_GIIntron_fold_bpp,
    tpr_other_ribozyme_fold_bpp,
    tpr_cili_telo_RNA_fold_bpp,
    tpr_RNasePRNA_fold_bpp,
    tpr_Vert_Telo_RNA_fold_bpp,
    tpr_ham_ribozyme_fold_bpp,
    tpr_SRPRNA_fold_bpp,
    tpr_YRNA_fold_bpp,
    tpr_cis_regulatory_element_fold_bpp,
    tpr_tRNA_fold_bpp,
    tpr_five_SrRNA_fold_bpp,
    tpr_other_RNA_fold_bpp,
    tpr_RNAIII_fold_bpp,
    tpr_RNaseE5UTR_fold_bpp,
    tpr_other_rRNA_fold_bpp,
    tpr_SyntheticRNA_fold_bpp,
    tpr_ires_fold_bpp,
    tpr_twentythree_SrRNA_fold_bpp,
    tpr_Viral_Phage_fold_bpp,
    tpr_snRNA_fold_bpp
    )

    fold_tpr_drtrs = (
    tpr_sixteen_SrRNA_drtr,
    tpr_GIIIntron_drtr,
    tpr_RNaseMRPRNA_drtr,
    tpr_hdvr_ribozyme_drtr,
    tpr_tmRNA_drtr,
    tpr_GIIntron_drtr,
    tpr_other_ribozyme_drtr,
    tpr_cili_telo_RNA_drtr,
    tpr_RNasePRNA_drtr,
    tpr_Vert_Telo_RNA_drtr,
    tpr_ham_ribozyme_drtr,
    tpr_SRPRNA_drtr,
    tpr_YRNA_drtr,
    tpr_cis_regulatory_element_drtr,
    tpr_tRNA_drtr,
    tpr_five_SrRNA_drtr,
    tpr_other_RNA_drtr,
    tpr_RNAIII_drtr,
    tpr_RNaseE5UTR_drtr,
    tpr_other_rRNA_drtr,
    tpr_SyntheticRNA_drtr,
    tpr_ires_drtr,
    tpr_twentythree_SrRNA_drtr,
    tpr_Viral_Phage_drtr,
    tpr_snRNA_drtr
    )

    fold_tpr_drtr_mfes = (
    tpr_sixteen_SrRNA_drtr_mfe,
    tpr_GIIIntron_drtr_mfe,
    tpr_RNaseMRPRNA_drtr_mfe,
    tpr_hdvr_ribozyme_drtr_mfe,
    tpr_tmRNA_drtr_mfe,
    tpr_GIIntron_drtr_mfe,
    tpr_other_ribozyme_drtr_mfe,
    tpr_cili_telo_RNA_drtr_mfe,
    tpr_RNasePRNA_drtr_mfe,
    tpr_Vert_Telo_RNA_drtr_mfe,
    tpr_ham_ribozyme_drtr_mfe,
    tpr_SRPRNA_drtr_mfe,
    tpr_YRNA_drtr_mfe,
    tpr_cis_regulatory_element_drtr_mfe,
    tpr_tRNA_drtr_mfe,
    tpr_five_SrRNA_drtr_mfe,
    tpr_other_RNA_drtr_mfe,
    tpr_RNAIII_drtr_mfe,
    tpr_RNaseE5UTR_drtr_mfe,
    tpr_other_rRNA_drtr_mfe,
    tpr_SyntheticRNA_drtr_mfe,
    tpr_ires_drtr_mfe,
    tpr_twentythree_SrRNA_drtr_mfe,
    tpr_Viral_Phage_drtr_mfe,
    tpr_snRNA_drtr_mfe
    )

    fold_tpr_drtr_probs = (
    tpr_sixteen_SrRNA_drtr_prob,
    tpr_GIIIntron_drtr_prob,
    tpr_RNaseMRPRNA_drtr_prob,
    tpr_hdvr_ribozyme_drtr_prob,
    tpr_tmRNA_drtr_prob,
    tpr_GIIntron_drtr_prob,
    tpr_other_ribozyme_drtr_prob,
    tpr_cili_telo_RNA_drtr_prob,
    tpr_RNasePRNA_drtr_prob,
    tpr_Vert_Telo_RNA_drtr_prob,
    tpr_ham_ribozyme_drtr_prob,
    tpr_SRPRNA_drtr_prob,
    tpr_YRNA_drtr_prob,
    tpr_cis_regulatory_element_drtr_prob,
    tpr_tRNA_drtr_prob,
    tpr_five_SrRNA_drtr_prob,
    tpr_other_RNA_drtr_prob,
    tpr_RNAIII_drtr_prob,
    tpr_RNaseE5UTR_drtr_prob,
    tpr_other_rRNA_drtr_prob,
    tpr_SyntheticRNA_drtr_prob,
    tpr_ires_drtr_prob,
    tpr_twentythree_SrRNA_drtr_prob,
    tpr_Viral_Phage_drtr_prob,
    tpr_snRNA_drtr_prob
    )




    fold_f1_mfes = [
    f1_sixteen_SrRNA_fold_mfe,                 #0
    f1_GIIIntron_fold_mfe,                     #5
    f1_RNaseMRPRNA_fold_mfe,                   #15
    f1_hdvr_ribozyme_fold_mfe,                 #8
    f1_tmRNA_fold_mfe,                         #20
    f1_GIIntron_fold_mfe,                      #6
    f1_other_ribozyme_fold_mfe,                #10
    f1_cili_telo_RNA_fold_mfe,                 #3
    f1_RNasePRNA_fold_mfe,                     #16
    f1_Vert_Telo_RNA_fold_mfe,                 #22
    f1_ham_ribozyme_fold_mfe,                  #7
    f1_SRPRNA_fold_mfe,                        #18
    f1_YRNA_fold_mfe,                          #24
    f1_cis_regulatory_element_fold_mfe,        #4
    f1_tRNA_fold_mfe,                          #21
    f1_five_SrRNA_fold_mfe,                    #2
    f1_other_RNA_fold_mfe,                     #11
    f1_RNAIII_fold_mfe,                        #13
    f1_RNaseE5UTR_fold_mfe,                    #14
    f1_other_rRNA_fold_mfe,                    #12
    f1_SyntheticRNA_fold_mfe,                  #19
    f1_ires_fold_mfe,                          #9
    f1_twentythree_SrRNA_fold_mfe,             #1
    f1_Viral_Phage_fold_mfe,                   #23
    f1_snRNA_fold_mfe                          #17
    ]

    fold_f1_bpps = (
    f1_sixteen_SrRNA_fold_bpp,
    f1_GIIIntron_fold_bpp,
    f1_RNaseMRPRNA_fold_bpp,
    f1_hdvr_ribozyme_fold_bpp,
    f1_tmRNA_fold_bpp,
    f1_GIIntron_fold_bpp,
    f1_other_ribozyme_fold_bpp,
    f1_cili_telo_RNA_fold_bpp,
    f1_RNasePRNA_fold_bpp,
    f1_Vert_Telo_RNA_fold_bpp,
    f1_ham_ribozyme_fold_bpp,
    f1_SRPRNA_fold_bpp,
    f1_YRNA_fold_bpp,
    f1_cis_regulatory_element_fold_bpp,
    f1_tRNA_fold_bpp,
    f1_five_SrRNA_fold_bpp,
    f1_other_RNA_fold_bpp,
    f1_RNAIII_fold_bpp,
    f1_RNaseE5UTR_fold_bpp,
    f1_other_rRNA_fold_bpp,
    f1_SyntheticRNA_fold_bpp,
    f1_ires_fold_bpp,
    f1_twentythree_SrRNA_fold_bpp,
    f1_Viral_Phage_fold_bpp,
    f1_snRNA_fold_bpp
    )

    fold_f1_drtrs = (
    f1_sixteen_SrRNA_drtr,
    f1_GIIIntron_drtr,
    f1_RNaseMRPRNA_drtr,
    f1_hdvr_ribozyme_drtr,
    f1_tmRNA_drtr,
    f1_GIIntron_drtr,
    f1_other_ribozyme_drtr,
    f1_cili_telo_RNA_drtr,
    f1_RNasePRNA_drtr,
    f1_Vert_Telo_RNA_drtr,
    f1_ham_ribozyme_drtr,
    f1_SRPRNA_drtr,
    f1_YRNA_drtr,
    f1_cis_regulatory_element_drtr,
    f1_tRNA_drtr,
    f1_five_SrRNA_drtr,
    f1_other_RNA_drtr,
    f1_RNAIII_drtr,
    f1_RNaseE5UTR_drtr,
    f1_other_rRNA_drtr,
    f1_SyntheticRNA_drtr,
    f1_ires_drtr,
    f1_twentythree_SrRNA_drtr,
    f1_Viral_Phage_drtr,
    f1_snRNA_drtr
    )

    fold_f1_drtr_mfes = (
    f1_sixteen_SrRNA_drtr_mfe,
    f1_GIIIntron_drtr_mfe,
    f1_RNaseMRPRNA_drtr_mfe,
    f1_hdvr_ribozyme_drtr_mfe,
    f1_tmRNA_drtr_mfe,
    f1_GIIntron_drtr_mfe,
    f1_other_ribozyme_drtr_mfe,
    f1_cili_telo_RNA_drtr_mfe,
    f1_RNasePRNA_drtr_mfe,
    f1_Vert_Telo_RNA_drtr_mfe,
    f1_ham_ribozyme_drtr_mfe,
    f1_SRPRNA_drtr_mfe,
    f1_YRNA_drtr_mfe,
    f1_cis_regulatory_element_drtr_mfe,
    f1_tRNA_drtr_mfe,
    f1_five_SrRNA_drtr_mfe,
    f1_other_RNA_drtr_mfe,
    f1_RNAIII_drtr_mfe,
    f1_RNaseE5UTR_drtr_mfe,
    f1_other_rRNA_drtr_mfe,
    f1_SyntheticRNA_drtr_mfe,
    f1_ires_drtr_mfe,
    f1_twentythree_SrRNA_drtr_mfe,
    f1_Viral_Phage_drtr_mfe,
    f1_snRNA_drtr_mfe
    )

    fold_f1_drtr_probs = (
    f1_sixteen_SrRNA_drtr_prob,
    f1_GIIIntron_drtr_prob,
    f1_RNaseMRPRNA_drtr_prob,
    f1_hdvr_ribozyme_drtr_prob,
    f1_tmRNA_drtr_prob,
    f1_GIIntron_drtr_prob,
    f1_other_ribozyme_drtr_prob,
    f1_cili_telo_RNA_drtr_prob,
    f1_RNasePRNA_drtr_prob,
    f1_Vert_Telo_RNA_drtr_prob,
    f1_ham_ribozyme_drtr_prob,
    f1_SRPRNA_drtr_prob,
    f1_YRNA_drtr_prob,
    f1_cis_regulatory_element_drtr_prob,
    f1_tRNA_drtr_prob,
    f1_five_SrRNA_drtr_prob,
    f1_other_RNA_drtr_prob,
    f1_RNAIII_drtr_prob,
    f1_RNaseE5UTR_drtr_prob,
    f1_other_rRNA_drtr_prob,
    f1_SyntheticRNA_drtr_prob,
    f1_ires_drtr_prob,
    f1_twentythree_SrRNA_drtr_prob,
    f1_Viral_Phage_drtr_prob,
    f1_snRNA_drtr_prob
    )

    print(output_name)
    data_output = [
    fold_ppv_mfes,
    fold_ppv_bpps,
    fold_ppv_drtr_mfes,
    fold_ppv_drtr_probs,
    fold_ppv_drtrs,
    fold_tpr_mfes,
    fold_tpr_bpps,
    fold_tpr_drtr_mfes,
    fold_tpr_drtr_probs,
    fold_tpr_drtrs,
    fold_f1_mfes,
    fold_f1_bpps,
    fold_f1_drtr_mfes,
    fold_f1_drtr_probs,
    fold_f1_drtrs,
    fold_mcc_mfes,
    fold_mcc_bpps,
    fold_mcc_drtr_mfes,
    fold_mcc_drtr_probs,
    fold_mcc_drtrs
    ]


fam_names =[
    family_dict[0],
    family_dict[5],
    family_dict[15],
    family_dict[8],
    family_dict[20],
    family_dict[6],
    family_dict[10],
    family_dict[3],
    family_dict[16],
    family_dict[22],
    family_dict[7],
    family_dict[18],
    family_dict[24],
    family_dict[4],
    family_dict[21],
    family_dict[2],
    family_dict[11],
    family_dict[13],
    family_dict[14],
    family_dict[12],
    family_dict[19],
    family_dict[9],
    family_dict[1],
    family_dict[23],
    family_dict[17]]
    # row_names = pd.Index()

if args.output==True:
    row_names = ["fold_mfe_ppv", "fold_bpp_ppv", "drtr_mfe_ppv", "drtr_prob_ppv", "drtr_bpp_ppv",
                 "fold_mfe_tpr", "fold_bpp_tpr", "drtr_mfe_tpr", "drtr_prob_tpr", "drtr_bpp_tpr",
                 "fold_mfe_f1", "fold_bpp_f1", "drtr_mfe_f1", "drtr_prob_f1", "drtr_bpp_f1",
                 "fold_mfe_mcc", "fold_bpp_mcc", "drtr_mfe_mcc", "drtr_prob_mcc", "drtr_bpp_mcc"]
    data_output = pd.DataFrame(data_output, columns=fam_names, index=row_names)
    print(data_output)
    export_csv = data_output.to_csv(output_name)






if args.printout==True:
    count_of_fam = len(tp_sixteen_SrRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[0],
                tpr_sixteen_SrRNA_fold_mfe,
                tpr_sixteen_SrRNA_fold_bpp,
                tpr_sixteen_SrRNA_drtr_mfe,
                tpr_sixteen_SrRNA_drtr,
                ppv_sixteen_SrRNA_fold_mfe,
                ppv_sixteen_SrRNA_fold_bpp,
                ppv_sixteen_SrRNA_drtr_mfe,
                ppv_sixteen_SrRNA_drtr,
                f1_sixteen_SrRNA_fold_mfe,
                f1_sixteen_SrRNA_fold_bpp,
                f1_sixteen_SrRNA_drtr_mfe,
                f1_sixteen_SrRNA_drtr,
                mcc_sixteen_SrRNA_fold_mfe,
                mcc_sixteen_SrRNA_fold_bpp,
                mcc_sixteen_SrRNA_drtr_mfe,
                mcc_sixteen_SrRNA_drtr,
                count_of_fam
                ))
    count_of_fam = len(tp_twentythree_SrRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[1],
                tpr_twentythree_SrRNA_fold_mfe,
                tpr_twentythree_SrRNA_fold_bpp,
                tpr_twentythree_SrRNA_drtr_mfe,
                tpr_twentythree_SrRNA_drtr,
                ppv_twentythree_SrRNA_fold_mfe,
                ppv_twentythree_SrRNA_fold_bpp,
                ppv_twentythree_SrRNA_drtr_mfe,
                ppv_twentythree_SrRNA_drtr,
                f1_twentythree_SrRNA_fold_mfe,
                f1_twentythree_SrRNA_fold_bpp,
                f1_twentythree_SrRNA_drtr_mfe,
                f1_twentythree_SrRNA_drtr,
                mcc_twentythree_SrRNA_fold_mfe,
                mcc_twentythree_SrRNA_fold_bpp,
                mcc_twentythree_SrRNA_drtr_mfe,
                mcc_twentythree_SrRNA_drtr,
                count_of_fam
                ))
    count_of_fam = len(tp_five_SrRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[2],
                tpr_five_SrRNA_fold_mfe,
                tpr_five_SrRNA_fold_bpp,
                tpr_five_SrRNA_drtr_mfe,
                tpr_five_SrRNA_drtr,
                ppv_five_SrRNA_fold_mfe,
                ppv_five_SrRNA_fold_bpp,
                ppv_five_SrRNA_drtr_mfe,
                ppv_five_SrRNA_drtr,
                f1_five_SrRNA_fold_mfe,
                f1_five_SrRNA_fold_bpp,
                f1_five_SrRNA_drtr_mfe,
                f1_five_SrRNA_drtr,
                mcc_five_SrRNA_fold_mfe,
                mcc_five_SrRNA_fold_bpp,
                mcc_five_SrRNA_drtr_mfe,
                mcc_five_SrRNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_cili_telo_RNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[3],
                tpr_cili_telo_RNA_fold_mfe,
                tpr_cili_telo_RNA_fold_bpp,
                tpr_cili_telo_RNA_drtr_mfe,
                tpr_cili_telo_RNA_drtr,
                ppv_cili_telo_RNA_fold_mfe,
                ppv_cili_telo_RNA_fold_bpp,
                ppv_cili_telo_RNA_drtr_mfe,
                ppv_cili_telo_RNA_drtr,
                f1_cili_telo_RNA_fold_mfe,
                f1_cili_telo_RNA_fold_bpp,
                f1_cili_telo_RNA_drtr_mfe,
                f1_cili_telo_RNA_drtr,
                mcc_cili_telo_RNA_fold_mfe,
                mcc_cili_telo_RNA_fold_bpp,
                mcc_cili_telo_RNA_drtr_mfe,
                mcc_cili_telo_RNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_cis_regulatory_element_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[4],
                tpr_cis_regulatory_element_fold_mfe,
                tpr_cis_regulatory_element_fold_bpp,
                tpr_cis_regulatory_element_drtr_mfe,
                tpr_cis_regulatory_element_drtr,
                ppv_cis_regulatory_element_fold_mfe,
                ppv_cis_regulatory_element_fold_bpp,
                ppv_cis_regulatory_element_drtr_mfe,
                ppv_cis_regulatory_element_drtr,
                f1_cis_regulatory_element_fold_mfe,
                f1_cis_regulatory_element_fold_bpp,
                f1_cis_regulatory_element_drtr_mfe,
                f1_cis_regulatory_element_drtr,
                mcc_cis_regulatory_element_fold_mfe,
                mcc_cis_regulatory_element_fold_bpp,
                mcc_cis_regulatory_element_drtr_mfe,
                mcc_cis_regulatory_element_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_GIIIntron_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[5],
                tpr_GIIIntron_fold_mfe,
                tpr_GIIIntron_fold_bpp,
                tpr_GIIIntron_drtr_mfe,
                tpr_GIIIntron_drtr,
                ppv_GIIIntron_fold_mfe,
                ppv_GIIIntron_fold_bpp,
                ppv_GIIIntron_drtr_mfe,
                ppv_GIIIntron_drtr,
                f1_GIIIntron_fold_mfe,
                f1_GIIIntron_fold_bpp,
                f1_GIIIntron_drtr_mfe,
                f1_GIIIntron_drtr,
                mcc_GIIIntron_fold_mfe,
                mcc_GIIIntron_fold_bpp,
                mcc_GIIIntron_drtr_mfe,
                mcc_GIIIntron_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_GIIntron_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[6],
                tpr_GIIntron_fold_mfe,
                tpr_GIIntron_fold_bpp,
                tpr_GIIntron_drtr_mfe,
                tpr_GIIntron_drtr,
                ppv_GIIntron_fold_mfe,
                ppv_GIIntron_fold_bpp,
                ppv_GIIntron_drtr_mfe,
                ppv_GIIntron_drtr,
                f1_GIIntron_fold_mfe,
                f1_GIIntron_fold_bpp,
                f1_GIIntron_drtr_mfe,
                f1_GIIntron_drtr,
                mcc_GIIntron_fold_mfe,
                mcc_GIIntron_fold_bpp,
                mcc_GIIntron_drtr_mfe,
                mcc_GIIntron_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_ham_ribozyme_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[7],
                tpr_ham_ribozyme_fold_mfe,
                tpr_ham_ribozyme_fold_bpp,
                tpr_ham_ribozyme_drtr_mfe,
                tpr_ham_ribozyme_drtr,
                ppv_ham_ribozyme_fold_mfe,
                ppv_ham_ribozyme_fold_bpp,
                ppv_ham_ribozyme_drtr_mfe,
                ppv_ham_ribozyme_drtr,
                f1_ham_ribozyme_fold_mfe,
                f1_ham_ribozyme_fold_bpp,
                f1_ham_ribozyme_drtr_mfe,
                f1_ham_ribozyme_drtr,
                mcc_ham_ribozyme_fold_mfe,
                mcc_ham_ribozyme_fold_bpp,
                mcc_ham_ribozyme_drtr_mfe,
                mcc_ham_ribozyme_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_hdvr_ribozyme_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[8],
                tpr_hdvr_ribozyme_fold_mfe,
                tpr_hdvr_ribozyme_fold_bpp,
                tpr_hdvr_ribozyme_drtr_mfe,
                tpr_hdvr_ribozyme_drtr,
                ppv_hdvr_ribozyme_fold_mfe,
                ppv_hdvr_ribozyme_fold_bpp,
                ppv_hdvr_ribozyme_drtr_mfe,
                ppv_hdvr_ribozyme_drtr,
                f1_hdvr_ribozyme_fold_mfe,
                f1_hdvr_ribozyme_fold_bpp,
                f1_hdvr_ribozyme_drtr_mfe,
                f1_hdvr_ribozyme_drtr,
                mcc_hdvr_ribozyme_fold_mfe,
                mcc_hdvr_ribozyme_fold_bpp,
                mcc_hdvr_ribozyme_drtr_mfe,
                mcc_hdvr_ribozyme_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_ires_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[9],
                tpr_ires_fold_mfe,
                tpr_ires_fold_bpp,
                tpr_ires_drtr_mfe,
                tpr_ires_drtr,
                ppv_ires_fold_mfe,
                ppv_ires_fold_bpp,
                ppv_ires_drtr_mfe,
                ppv_ires_drtr,
                f1_ires_fold_mfe,
                f1_ires_fold_bpp,
                f1_ires_drtr_mfe,
                f1_ires_drtr,
                mcc_ires_fold_mfe,
                mcc_ires_fold_bpp,
                mcc_ires_drtr_mfe,
                mcc_ires_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_other_ribozyme_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[10],
                tpr_other_ribozyme_fold_mfe,
                tpr_other_ribozyme_fold_bpp,
                tpr_other_ribozyme_drtr_mfe,
                tpr_other_ribozyme_drtr,
                ppv_other_ribozyme_fold_mfe,
                ppv_other_ribozyme_fold_bpp,
                ppv_other_ribozyme_drtr_mfe,
                ppv_other_ribozyme_drtr,
                f1_other_ribozyme_fold_mfe,
                f1_other_ribozyme_fold_bpp,
                f1_other_ribozyme_drtr_mfe,
                f1_other_ribozyme_drtr,
                mcc_other_ribozyme_fold_mfe,
                mcc_other_ribozyme_fold_bpp,
                mcc_other_ribozyme_drtr_mfe,
                mcc_other_ribozyme_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_other_RNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[11],
                tpr_other_RNA_fold_mfe,
                tpr_other_RNA_fold_bpp,
                tpr_other_RNA_drtr_mfe,
                tpr_other_RNA_drtr,
                ppv_other_RNA_fold_mfe,
                ppv_other_RNA_fold_bpp,
                ppv_other_RNA_drtr_mfe,
                ppv_other_RNA_drtr,
                f1_other_RNA_fold_mfe,
                f1_other_RNA_fold_bpp,
                f1_other_RNA_drtr_mfe,
                f1_other_RNA_drtr,
                mcc_other_RNA_fold_mfe,
                mcc_other_RNA_fold_bpp,
                mcc_other_RNA_drtr_mfe,
                mcc_other_RNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_other_rRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[12],
                tpr_other_rRNA_fold_mfe,
                tpr_other_rRNA_fold_bpp,
                tpr_other_rRNA_drtr_mfe,
                tpr_other_rRNA_drtr,
                ppv_other_rRNA_fold_mfe,
                ppv_other_rRNA_fold_bpp,
                ppv_other_rRNA_drtr_mfe,
                ppv_other_rRNA_drtr,
                f1_other_rRNA_fold_mfe,
                f1_other_rRNA_fold_bpp,
                f1_other_rRNA_drtr_mfe,
                f1_other_rRNA_drtr,
                mcc_other_rRNA_fold_mfe,
                mcc_other_rRNA_fold_bpp,
                mcc_other_rRNA_drtr_mfe,
                mcc_other_rRNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_RNAIII_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[13],
                tpr_RNAIII_fold_mfe,
                tpr_RNAIII_fold_bpp,
                tpr_RNAIII_drtr_mfe,
                tpr_RNAIII_drtr,
                ppv_RNAIII_fold_mfe,
                ppv_RNAIII_fold_bpp,
                ppv_RNAIII_drtr_mfe,
                ppv_RNAIII_drtr,
                f1_RNAIII_fold_mfe,
                f1_RNAIII_fold_bpp,
                f1_RNAIII_drtr_mfe,
                f1_RNAIII_drtr,
                mcc_RNAIII_fold_mfe,
                mcc_RNAIII_fold_bpp,
                mcc_RNAIII_drtr_mfe,
                mcc_RNAIII_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_RNaseE5UTR_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[14],
                tpr_RNaseE5UTR_fold_mfe,
                tpr_RNaseE5UTR_fold_bpp,
                tpr_RNaseE5UTR_drtr_mfe,
                tpr_RNaseE5UTR_drtr,
                ppv_RNaseE5UTR_fold_mfe,
                ppv_RNaseE5UTR_fold_bpp,
                ppv_RNaseE5UTR_drtr_mfe,
                ppv_RNaseE5UTR_drtr,
                f1_RNaseE5UTR_fold_mfe,
                f1_RNaseE5UTR_fold_bpp,
                f1_RNaseE5UTR_drtr_mfe,
                f1_RNaseE5UTR_drtr,
                mcc_RNaseE5UTR_fold_mfe,
                mcc_RNaseE5UTR_fold_bpp,
                mcc_RNaseE5UTR_drtr_mfe,
                mcc_RNaseE5UTR_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_RNaseMRPRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[15],
                tpr_RNaseMRPRNA_fold_mfe,
                tpr_RNaseMRPRNA_fold_bpp,
                tpr_RNaseMRPRNA_drtr_mfe,
                tpr_RNaseMRPRNA_drtr,
                ppv_RNaseMRPRNA_fold_mfe,
                ppv_RNaseMRPRNA_fold_bpp,
                ppv_RNaseMRPRNA_drtr_mfe,
                ppv_RNaseMRPRNA_drtr,
                f1_RNaseMRPRNA_fold_mfe,
                f1_RNaseMRPRNA_fold_bpp,
                f1_RNaseMRPRNA_drtr_mfe,
                f1_RNaseMRPRNA_drtr,
                mcc_RNaseMRPRNA_fold_mfe,
                mcc_RNaseMRPRNA_fold_bpp,
                mcc_RNaseMRPRNA_drtr_mfe,
                mcc_RNaseMRPRNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_RNasePRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[16],
                tpr_RNasePRNA_fold_mfe,
                tpr_RNasePRNA_fold_bpp,
                tpr_RNasePRNA_drtr_mfe,
                tpr_RNasePRNA_drtr,
                ppv_RNasePRNA_fold_mfe,
                ppv_RNasePRNA_fold_bpp,
                ppv_RNasePRNA_drtr,
                ppv_RNasePRNA_drtr_mfe,
                f1_RNasePRNA_fold_mfe,
                f1_RNasePRNA_fold_bpp,
                f1_RNasePRNA_drtr_mfe,
                f1_RNasePRNA_drtr,
                mcc_RNasePRNA_fold_mfe,
                mcc_RNasePRNA_fold_bpp,
                mcc_RNasePRNA_drtr_mfe,
                mcc_RNasePRNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_snRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[17],
                tpr_snRNA_fold_mfe,
                tpr_snRNA_fold_bpp,
                tpr_snRNA_drtr_mfe,
                tpr_snRNA_drtr,
                ppv_snRNA_fold_mfe,
                ppv_snRNA_fold_bpp,
                ppv_snRNA_drtr_mfe,
                ppv_snRNA_drtr,
                f1_snRNA_fold_mfe,
                f1_snRNA_fold_bpp,
                f1_snRNA_drtr_mfe,
                f1_snRNA_drtr,
                mcc_snRNA_fold_mfe,
                mcc_snRNA_fold_bpp,
                mcc_snRNA_drtr_mfe,
                mcc_snRNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_SRPRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[18],
                tpr_SRPRNA_fold_mfe,
                tpr_SRPRNA_fold_bpp,
                tpr_SRPRNA_drtr_mfe,
                tpr_SRPRNA_drtr,
                ppv_SRPRNA_fold_mfe,
                ppv_SRPRNA_fold_bpp,
                ppv_SRPRNA_drtr_mfe,
                ppv_SRPRNA_drtr,
                f1_SRPRNA_fold_mfe,
                f1_SRPRNA_fold_bpp,
                f1_SRPRNA_drtr_mfe,
                f1_SRPRNA_drtr,
                mcc_SRPRNA_fold_mfe,
                mcc_SRPRNA_fold_bpp,
                mcc_SRPRNA_drtr_mfe,
                mcc_SRPRNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_SyntheticRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[19],
                tpr_SyntheticRNA_fold_mfe,
                tpr_SyntheticRNA_fold_bpp,
                tpr_SyntheticRNA_drtr_mfe,
                tpr_SyntheticRNA_drtr,
                ppv_SyntheticRNA_fold_mfe,
                ppv_SyntheticRNA_fold_bpp,
                ppv_SyntheticRNA_drtr_mfe,
                ppv_SyntheticRNA_drtr,
                f1_SyntheticRNA_fold_mfe,
                f1_SyntheticRNA_fold_bpp,
                f1_SyntheticRNA_drtr_mfe,
                f1_SyntheticRNA_drtr,
                mcc_SyntheticRNA_fold_mfe,
                mcc_SyntheticRNA_fold_bpp,
                mcc_SyntheticRNA_drtr_mfe,
                mcc_SyntheticRNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_tmRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[20],
                tpr_tmRNA_fold_mfe,
                tpr_tmRNA_fold_bpp,
                tpr_tmRNA_drtr_mfe,
                tpr_tmRNA_drtr,
                ppv_tmRNA_fold_mfe,
                ppv_tmRNA_fold_bpp,
                ppv_tmRNA_drtr_mfe,
                ppv_tmRNA_drtr,
                f1_tmRNA_fold_mfe,
                f1_tmRNA_fold_bpp,
                f1_tmRNA_drtr_mfe,
                f1_tmRNA_drtr,
                mcc_tmRNA_fold_mfe,
                mcc_tmRNA_fold_bpp,
                mcc_tmRNA_drtr_mfe,
                mcc_tmRNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_tRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[21],
                tpr_tRNA_fold_mfe,
                tpr_tRNA_fold_bpp,
                tpr_tRNA_drtr_mfe,
                tpr_tRNA_drtr,
                ppv_tRNA_fold_mfe,
                ppv_tRNA_fold_bpp,
                ppv_tRNA_drtr_mfe,
                ppv_tRNA_drtr,
                f1_tRNA_fold_mfe,
                f1_tRNA_fold_bpp,
                f1_tRNA_drtr_mfe,
                f1_tRNA_drtr,
                mcc_tRNA_fold_mfe,
                mcc_tRNA_fold_bpp,
                mcc_tRNA_drtr_mfe,
                mcc_tRNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_Vert_Telo_RNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[22],
                tpr_Vert_Telo_RNA_fold_mfe,
                tpr_Vert_Telo_RNA_fold_bpp,
                tpr_Vert_Telo_RNA_drtr_mfe,
                tpr_Vert_Telo_RNA_drtr,
                ppv_Vert_Telo_RNA_fold_mfe,
                ppv_Vert_Telo_RNA_fold_bpp,
                ppv_Vert_Telo_RNA_drtr_mfe,
                ppv_Vert_Telo_RNA_drtr,
                f1_Vert_Telo_RNA_fold_mfe,
                f1_Vert_Telo_RNA_fold_bpp,
                f1_Vert_Telo_RNA_drtr_mfe,
                f1_Vert_Telo_RNA_drtr,
                mcc_Vert_Telo_RNA_fold_mfe,
                mcc_Vert_Telo_RNA_fold_bpp,
                mcc_Vert_Telo_RNA_drtr_mfe,
                mcc_Vert_Telo_RNA_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_Viral_Phage_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[23],
                tpr_Viral_Phage_fold_mfe,
                tpr_Viral_Phage_fold_bpp,
                tpr_Viral_Phage_drtr_mfe,
                tpr_Viral_Phage_drtr,
                ppv_Viral_Phage_fold_mfe,
                ppv_Viral_Phage_fold_bpp,
                ppv_Viral_Phage_drtr_mfe,
                ppv_Viral_Phage_drtr,
                f1_Viral_Phage_fold_mfe,
                f1_Viral_Phage_fold_bpp,
                f1_Viral_Phage_drtr_mfe,
                f1_Viral_Phage_drtr,
                mcc_Viral_Phage_fold_mfe,
                mcc_Viral_Phage_fold_bpp,
                mcc_Viral_Phage_drtr_mfe,
                mcc_Viral_Phage_drtr,
                count_of_fam
                ))

    count_of_fam = len(tp_YRNA_drtr)
    print(
          "%s\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%8.5f\t"
          "%10i\t"
           % (
                family_dict[24],
                tpr_YRNA_fold_mfe,
                tpr_YRNA_fold_bpp,
                tpr_YRNA_drtr_mfe,
                tpr_YRNA_drtr,
                ppv_YRNA_fold_mfe,
                ppv_YRNA_fold_bpp,
                ppv_YRNA_drtr_mfe,
                ppv_YRNA_drtr,
                f1_YRNA_fold_mfe,
                f1_YRNA_fold_bpp,
                f1_YRNA_drtr_mfe,
                f1_YRNA_drtr,
                mcc_YRNA_fold_mfe,
                mcc_YRNA_fold_bpp,
                mcc_YRNA_drtr_mfe,
                mcc_YRNA_drtr,
                count_of_fam
                ))
