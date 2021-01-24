#!python


### ======================================= IMPORT ====================================== ###
import numpy as np
import pandas as pd
import math
import sys
import matplotlib.pyplot as plt
import argparse
from matplotlib.pyplot import figure
from random import seed
from random import randint

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
                    default="comp_m_n_fin.csv",
                    help='parse directory')

parser.add_argument('--output','-o', action='store_true',
                    default=False,
                    help='recalculate boot/csvs')

args=parser.parse_args()


if args.test == False:
    input_name = args.input

if args.test == True:
    input_name = "comp_m_cp_mfe_prob.csv"

output_name = input_name[:-4] + "_f_boot.csv"


### ================================ FUNCTIONS ==================================== ###
def tpr(tp_data,fn_data):
    result = 0
    list_r = list()
    for i in range(0, len(tp_data)):
        result = float(tp_data[i]) / (float(tp_data[i]) + float(fn_data[i]))
        list_r.append(result)
    return list_r


def ppv(tp_data,fp_data):
    result = 0
    list_r = list()
    for i in range(0, len(tp_data)):
        result = float(tp_data[i]) / (float(tp_data[i]) + float(fp_data[i]))
        list_r.append(result)
    return list_r


def f_1(tp_data, fp_data, fn_data):
    result = 0
    list_r = list()
    f = 0
    for i in range(0, len(tp_data)):
        ppv = float(tp_data[i]) / (float(tp_data[i]) + float(fp_data[i]))
        tpr = float(tp_data[i]) / (float(tp_data[i]) + float(fn_data[i]))
        if (ppv + tpr) > 0:
            f = 2 * (ppv * tpr)/(ppv + tpr)
        else: f = 0
        result = f
        list_r.append(result)
    return list_r


def mcc(tp_data, fp_data, tn_data, fn_data):
    result = 0
    list_r = list()
    for i in range (0, len(tp_data)):
        x = ((float(tp_data[i]) + float(fp_data[i])) * (float(tp_data[i]) + float(fn_data[i])) * (float(tn_data[i]) + float(fp_data[i])) * (float(tn_data[i]) + float(fn_data[i])))
        # x = ((np.float64(tp_data[i]) + fp_data[i]) * (tp_data[i] + fn_data[i]) * (tn_data[i] + fp_data[i]) * (tn_data[i] + fn_data[i]))
        # # print(x)
        root = math.sqrt(x)
        mcc = (float(tp_data[i]) * float(tn_data[i]) - float(fp_data[i]) * float(fn_data[i])) / root
        list_r.append(mcc)
        # # print(tp_data[i])
    return list_r

if args.test == True:
    boots = 10
if args.test == False:
    boots = 1000
def boot(data):
    n = len(data)
    value = 0
    sum = 0
    list_v = list()
    list_s = list()
    for i in range(0, boots):
        for j in range(0, n - 1):
            rng = randint(0, n - 1)
            # print(rng)
            value = (data.iloc[rng])
            list_v.append(value)
        # print(list_v)
        for k in range(0, len(list_v)):
            sum += list_v[k]
        list_s.append(sum)
        list_v = list()
        sum = 0
    return list_s
    # return sum


### ====================================== IMPORT =================================== ###
raw_data = pd.read_csv( open(input_name),float_precision = "high")
data_length = len(raw_data)

# print(raw_data)

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
RNAfold_mfe_tp = raw_data['RNAfold_mfe_tp']
RNAfold_mfe_fp = raw_data['RNAfold_mfe_fp']
RNAfold_mfe_tn = raw_data['RNAfold_mfe_tn']
RNAfold_mfe_fn = raw_data['RNAfold_mfe_fn']
RNAfold_bpp_tp = raw_data['RNAfold_bpp_tp']
RNAfold_bpp_fp = raw_data['RNAfold_bpp_fp']
RNAfold_bpp_tn = raw_data['RNAfold_bpp_tn']
RNAfold_bpp_fn = raw_data['RNAfold_bpp_fn']
lendata = len(pps_fold)


### ==================================== split by family ============================ ###
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

# print(family_dict[1])

# for i in range(0, len(family_dict)):


#16SrRNA _fold, _drtr
sixteen_SrRNA=family_data.loc[family_data.family=='16SrRNA']

tp_sixteen_SrRNA_fold_mfe = (sixteen_SrRNA[['RNAfold_mfe_tp']])
tp_sixteen_SrRNA_fold_mfe_sum = boot(tp_sixteen_SrRNA_fold_mfe)

fp_sixteen_SrRNA_fold_mfe = (sixteen_SrRNA[['RNAfold_mfe_fp']])
fp_sixteen_SrRNA_fold_mfe_sum = boot(fp_sixteen_SrRNA_fold_mfe)

tn_sixteen_SrRNA_fold_mfe = (sixteen_SrRNA[['RNAfold_mfe_tn']])
tn_sixteen_SrRNA_fold_mfe_sum = boot(tn_sixteen_SrRNA_fold_mfe)

fn_sixteen_SrRNA_fold_mfe = (sixteen_SrRNA[['RNAfold_mfe_fn']])
fn_sixteen_SrRNA_fold_mfe_sum = boot(fn_sixteen_SrRNA_fold_mfe)


tp_sixteen_SrRNA_fold_bpp = (sixteen_SrRNA[['RNAfold_bpp_tp']])
tp_sixteen_SrRNA_fold_bpp_sum = boot(tp_sixteen_SrRNA_fold_bpp)

fp_sixteen_SrRNA_fold_bpp = (sixteen_SrRNA[['RNAfold_bpp_fp']])
fp_sixteen_SrRNA_fold_bpp_sum = boot(fp_sixteen_SrRNA_fold_bpp)

tn_sixteen_SrRNA_fold_bpp = (sixteen_SrRNA[['RNAfold_bpp_tn']])
tn_sixteen_SrRNA_fold_bpp_sum = boot(tn_sixteen_SrRNA_fold_bpp)

fn_sixteen_SrRNA_fold_bpp = (sixteen_SrRNA[['RNAfold_bpp_fn']])
fn_sixteen_SrRNA_fold_bpp_sum = boot(fn_sixteen_SrRNA_fold_bpp)


tp_sixteen_SrRNA_drtr_mfe = (sixteen_SrRNA[['drtr_tp_mfe']])
tp_sixteen_SrRNA_drtr_mfe_sum = boot(tp_sixteen_SrRNA_drtr_mfe)

fp_sixteen_SrRNA_drtr_mfe = (sixteen_SrRNA[['drtr_fp_mfe']])
fp_sixteen_SrRNA_drtr_mfe_sum = boot(fp_sixteen_SrRNA_drtr_mfe)

tn_sixteen_SrRNA_drtr_mfe = (sixteen_SrRNA[['drtr_tn_mfe']])
tn_sixteen_SrRNA_drtr_mfe_sum = boot(tn_sixteen_SrRNA_drtr_mfe)

fn_sixteen_SrRNA_drtr_mfe = (sixteen_SrRNA[['drtr_fn_mfe']])
fn_sixteen_SrRNA_drtr_mfe_sum = boot(fn_sixteen_SrRNA_drtr_mfe)


tp_sixteen_SrRNA_drtr_prob = (sixteen_SrRNA[['drtr_tp_prob']])
tp_sixteen_SrRNA_drtr_prob_sum = boot(tp_sixteen_SrRNA_drtr_prob)

fp_sixteen_SrRNA_drtr_prob = (sixteen_SrRNA[['drtr_fp_prob']])
fp_sixteen_SrRNA_drtr_prob_sum = boot(fp_sixteen_SrRNA_drtr_prob)

tn_sixteen_SrRNA_drtr_prob = (sixteen_SrRNA[['drtr_tn_prob']])
tn_sixteen_SrRNA_drtr_prob_sum = boot(tn_sixteen_SrRNA_drtr_prob)

fn_sixteen_SrRNA_drtr_prob = (sixteen_SrRNA[['drtr_fn_prob']])
fn_sixteen_SrRNA_drtr_prob_sum = boot(fn_sixteen_SrRNA_drtr_prob)


tp_sixteen_SrRNA_drtr_bpp = (sixteen_SrRNA[['drtr_tp']])
tp_sixteen_SrRNA_drtr_bpp_sum = boot(tp_sixteen_SrRNA_drtr_bpp)

fp_sixteen_SrRNA_drtr_bpp = (sixteen_SrRNA[['drtr_fp']])
fp_sixteen_SrRNA_drtr_bpp_sum = boot(fp_sixteen_SrRNA_drtr_bpp)

tn_sixteen_SrRNA_drtr_bpp = (sixteen_SrRNA[['drtr_tn']])
tn_sixteen_SrRNA_drtr_bpp_sum = boot(tn_sixteen_SrRNA_drtr_bpp)

fn_sixteen_SrRNA_drtr_bpp = (sixteen_SrRNA[['drtr_fn']])
fn_sixteen_SrRNA_drtr_bpp_sum = boot(fn_sixteen_SrRNA_drtr_bpp)

#TPR
tpr_sixteen_SrRNA_fold_mfe = tpr(tp_sixteen_SrRNA_fold_mfe_sum, fn_sixteen_SrRNA_fold_mfe_sum)
# print(tpr_sixteen_SrRNA_fold_mfe)

tpr_sixteen_SrRNA_fold_bpp = tpr(tp_sixteen_SrRNA_fold_bpp_sum, fn_sixteen_SrRNA_fold_bpp_sum)
# print(tpr_sixteen_SrRNA_fold_bpp)

tpr_sixteen_SrRNA_drtr_mfe = tpr(tp_sixteen_SrRNA_drtr_mfe_sum, fn_sixteen_SrRNA_drtr_mfe_sum)
# print(tpr_sixteen_SrRNA_drtr_mfe)

tpr_sixteen_SrRNA_drtr_prob = tpr(tp_sixteen_SrRNA_drtr_prob_sum, fn_sixteen_SrRNA_drtr_prob_sum)
# print(tpr_sixteen_SrRNA_drtr_prob)

tpr_sixteen_SrRNA_drtr_bpp = tpr(tp_sixteen_SrRNA_drtr_bpp_sum, fn_sixteen_SrRNA_drtr_bpp_sum)
# print(tpr_sixteen_SrRNA_drtr_bpp)


#PPV
ppv_sixteen_SrRNA_fold_mfe = ppv(tp_sixteen_SrRNA_fold_mfe_sum, fn_sixteen_SrRNA_fold_mfe_sum)
# print(ppv_sixteen_SrRNA_fold_mfe)

ppv_sixteen_SrRNA_fold_bpp = ppv(tp_sixteen_SrRNA_fold_bpp_sum, fn_sixteen_SrRNA_fold_bpp_sum)
# print(ppv_sixteen_SrRNA_fold_bpp)

ppv_sixteen_SrRNA_drtr_mfe = ppv(tp_sixteen_SrRNA_drtr_mfe_sum, fn_sixteen_SrRNA_drtr_mfe_sum)
# print(ppv_sixteen_SrRNA_drtr_mfe)

ppv_sixteen_SrRNA_drtr_prob = ppv(tp_sixteen_SrRNA_drtr_prob_sum, fn_sixteen_SrRNA_drtr_prob_sum)
# print(ppv_sixteen_SrRNA_drtr_prob)

ppv_sixteen_SrRNA_drtr_bpp = ppv(tp_sixteen_SrRNA_drtr_bpp_sum, fn_sixteen_SrRNA_drtr_bpp_sum)
# print(ppv_sixteen_SrRNA_drtr_bpp)

#F1
f_1_sixteen_SrRNA_fold_mfe = f_1(tp_sixteen_SrRNA_fold_mfe_sum, fp_sixteen_SrRNA_fold_mfe_sum, fn_sixteen_SrRNA_fold_mfe_sum)
# print(f_1_sixteen_SrRNA_fold_mfe)

f_1_sixteen_SrRNA_fold_bpp = f_1(tp_sixteen_SrRNA_fold_bpp_sum, fp_sixteen_SrRNA_fold_bpp_sum, fn_sixteen_SrRNA_fold_bpp_sum)
# print(f_1_sixteen_SrRNA_fold_bpp)

f_1_sixteen_SrRNA_drtr_mfe = f_1(tp_sixteen_SrRNA_drtr_mfe_sum, fp_sixteen_SrRNA_drtr_mfe_sum, fn_sixteen_SrRNA_drtr_mfe_sum)
# print(f_1_sixteen_SrRNA_drtr_mfe)

f_1_sixteen_SrRNA_drtr_prob = f_1(tp_sixteen_SrRNA_drtr_prob_sum, fp_sixteen_SrRNA_drtr_prob_sum, fn_sixteen_SrRNA_drtr_prob_sum)
# print(f_1_sixteen_SrRNA_drtr_prob)

f_1_sixteen_SrRNA_drtr_bpp = f_1(tp_sixteen_SrRNA_drtr_bpp_sum, fp_sixteen_SrRNA_drtr_bpp_sum, fn_sixteen_SrRNA_drtr_bpp_sum)
# print(f_1_sixteen_SrRNA_drtr_bpp)

# #MCC
mcc_sixteen_SrRNA_fold_mfe = mcc(tp_sixteen_SrRNA_fold_mfe_sum, fp_sixteen_SrRNA_fold_mfe_sum, tn_sixteen_SrRNA_fold_mfe_sum, fn_sixteen_SrRNA_fold_mfe_sum)
# print(mcc_sixteen_SrRNA_fold_mfe)

mcc_sixteen_SrRNA_fold_bpp = mcc(tp_sixteen_SrRNA_fold_bpp_sum, fp_sixteen_SrRNA_fold_bpp_sum, tn_sixteen_SrRNA_fold_bpp_sum, fn_sixteen_SrRNA_fold_bpp_sum)
# print(mcc_sixteen_SrRNA_fold_bpp)

mcc_sixteen_SrRNA_drtr_mfe = mcc(tp_sixteen_SrRNA_drtr_mfe_sum, fp_sixteen_SrRNA_drtr_mfe_sum, tn_sixteen_SrRNA_drtr_mfe_sum, fn_sixteen_SrRNA_drtr_mfe_sum)
# print(mcc_sixteen_SrRNA_drtr_mfe)

mcc_sixteen_SrRNA_drtr_prob = mcc(tp_sixteen_SrRNA_drtr_prob_sum, fp_sixteen_SrRNA_drtr_prob_sum, tn_sixteen_SrRNA_drtr_prob_sum, fn_sixteen_SrRNA_drtr_prob_sum)
# print(mcc_sixteen_SrRNA_drtr_prob)

mcc_sixteen_SrRNA_drtr_bpp = mcc(tp_sixteen_SrRNA_drtr_bpp_sum, fp_sixteen_SrRNA_drtr_bpp_sum, tn_sixteen_SrRNA_drtr_bpp_sum, fn_sixteen_SrRNA_drtr_bpp_sum)
# print(mcc_sixteen_SrRNA_drtr_bpp)

### ============================================================== 23SrRNA
twentythree_SrRNA=family_data.loc[family_data.family=='23SrRNA']

#tp_twentythree_SrRNA_fold_mfe = (twentythree_SrRNA[['RNAfold_mfe_tp']])

tp_twentythree_SrRNA_fold_mfe = (twentythree_SrRNA[['RNAfold_mfe_tp']])
tp_twentythree_SrRNA_fold_mfe_sum = boot(tp_twentythree_SrRNA_fold_mfe)

fp_twentythree_SrRNA_fold_mfe = (twentythree_SrRNA[['RNAfold_mfe_fp']])
fp_twentythree_SrRNA_fold_mfe_sum = boot(fp_twentythree_SrRNA_fold_mfe)

tn_twentythree_SrRNA_fold_mfe = (twentythree_SrRNA[['RNAfold_mfe_tn']])
tn_twentythree_SrRNA_fold_mfe_sum = boot(tn_twentythree_SrRNA_fold_mfe)

fn_twentythree_SrRNA_fold_mfe = (twentythree_SrRNA[['RNAfold_mfe_fn']])
fn_twentythree_SrRNA_fold_mfe_sum = boot(fn_twentythree_SrRNA_fold_mfe)


tp_twentythree_SrRNA_fold_bpp = (twentythree_SrRNA[['RNAfold_bpp_tp']])
tp_twentythree_SrRNA_fold_bpp_sum = boot(tp_twentythree_SrRNA_fold_bpp)

fp_twentythree_SrRNA_fold_bpp = (twentythree_SrRNA[['RNAfold_bpp_fp']])
fp_twentythree_SrRNA_fold_bpp_sum = boot(fp_twentythree_SrRNA_fold_bpp)

tn_twentythree_SrRNA_fold_bpp = (twentythree_SrRNA[['RNAfold_bpp_tn']])
tn_twentythree_SrRNA_fold_bpp_sum = boot(tn_twentythree_SrRNA_fold_bpp)

fn_twentythree_SrRNA_fold_bpp = (twentythree_SrRNA[['RNAfold_bpp_fn']])
fn_twentythree_SrRNA_fold_bpp_sum = boot(fn_twentythree_SrRNA_fold_bpp)


tp_twentythree_SrRNA_drtr_mfe = (twentythree_SrRNA[['drtr_tp_mfe']])
tp_twentythree_SrRNA_drtr_mfe_sum = boot(tp_twentythree_SrRNA_drtr_mfe)

fp_twentythree_SrRNA_drtr_mfe = (twentythree_SrRNA[['drtr_fp_mfe']])
fp_twentythree_SrRNA_drtr_mfe_sum = boot(fp_twentythree_SrRNA_drtr_mfe)

tn_twentythree_SrRNA_drtr_mfe = (twentythree_SrRNA[['drtr_tn_mfe']])
tn_twentythree_SrRNA_drtr_mfe_sum = boot(tn_twentythree_SrRNA_drtr_mfe)

fn_twentythree_SrRNA_drtr_mfe = (twentythree_SrRNA[['drtr_fn_mfe']])
fn_twentythree_SrRNA_drtr_mfe_sum = boot(fn_twentythree_SrRNA_drtr_mfe)


tp_twentythree_SrRNA_drtr_prob = (twentythree_SrRNA[['drtr_tp_prob']])
tp_twentythree_SrRNA_drtr_prob_sum = boot(tp_twentythree_SrRNA_drtr_prob)

fp_twentythree_SrRNA_drtr_prob = (twentythree_SrRNA[['drtr_fp_prob']])
fp_twentythree_SrRNA_drtr_prob_sum = boot(fp_twentythree_SrRNA_drtr_prob)

tn_twentythree_SrRNA_drtr_prob = (twentythree_SrRNA[['drtr_tn_prob']])
tn_twentythree_SrRNA_drtr_prob_sum = boot(tn_twentythree_SrRNA_drtr_prob)

fn_twentythree_SrRNA_drtr_prob = (twentythree_SrRNA[['drtr_fn_prob']])
fn_twentythree_SrRNA_drtr_prob_sum = boot(fn_twentythree_SrRNA_drtr_prob)


tp_twentythree_SrRNA_drtr_bpp = (twentythree_SrRNA[['drtr_tp']])
tp_twentythree_SrRNA_drtr_bpp_sum = boot(tp_twentythree_SrRNA_drtr_bpp)

fp_twentythree_SrRNA_drtr_bpp = (twentythree_SrRNA[['drtr_fp']])
fp_twentythree_SrRNA_drtr_bpp_sum = boot(fp_twentythree_SrRNA_drtr_bpp)

tn_twentythree_SrRNA_drtr_bpp = (twentythree_SrRNA[['drtr_tn']])
tn_twentythree_SrRNA_drtr_bpp_sum = boot(tn_twentythree_SrRNA_drtr_bpp)

fn_twentythree_SrRNA_drtr_bpp = (twentythree_SrRNA[['drtr_fn']])
fn_twentythree_SrRNA_drtr_bpp_sum = boot(fn_twentythree_SrRNA_drtr_bpp)

#TPR
tpr_twentythree_SrRNA_fold_mfe = tpr(tp_twentythree_SrRNA_fold_mfe_sum, fn_twentythree_SrRNA_fold_mfe_sum)
# print(tpr_twentythree_SrRNA_fold_mfe)

tpr_twentythree_SrRNA_fold_bpp = tpr(tp_twentythree_SrRNA_fold_bpp_sum, fn_twentythree_SrRNA_fold_bpp_sum)
# print(tpr_twentythree_SrRNA_fold_bpp)

tpr_twentythree_SrRNA_drtr_mfe = tpr(tp_twentythree_SrRNA_drtr_mfe_sum, fn_twentythree_SrRNA_drtr_mfe_sum)
# print(tpr_twentythree_SrRNA_drtr_mfe)

tpr_twentythree_SrRNA_drtr_prob = tpr(tp_twentythree_SrRNA_drtr_prob_sum, fn_twentythree_SrRNA_drtr_prob_sum)
# print(tpr_twentythree_SrRNA_drtr_prob)

tpr_twentythree_SrRNA_drtr_bpp = tpr(tp_twentythree_SrRNA_drtr_bpp_sum, fn_twentythree_SrRNA_drtr_bpp_sum)
# print(tpr_twentythree_SrRNA_drtr_bpp)


#PPV
ppv_twentythree_SrRNA_fold_mfe = ppv(tp_twentythree_SrRNA_fold_mfe_sum, fn_twentythree_SrRNA_fold_mfe_sum)
# print(ppv_twentythree_SrRNA_fold_mfe)

ppv_twentythree_SrRNA_fold_bpp = ppv(tp_twentythree_SrRNA_fold_bpp_sum, fn_twentythree_SrRNA_fold_bpp_sum)
# print(ppv_twentythree_SrRNA_fold_bpp)

ppv_twentythree_SrRNA_drtr_mfe = ppv(tp_twentythree_SrRNA_drtr_mfe_sum, fn_twentythree_SrRNA_drtr_mfe_sum)
# print(ppv_twentythree_SrRNA_drtr_mfe)

ppv_twentythree_SrRNA_drtr_prob = ppv(tp_twentythree_SrRNA_drtr_prob_sum, fn_twentythree_SrRNA_drtr_prob_sum)
# print(ppv_twentythree_SrRNA_drtr_prob)

ppv_twentythree_SrRNA_drtr_bpp = ppv(tp_twentythree_SrRNA_drtr_bpp_sum, fn_twentythree_SrRNA_drtr_bpp_sum)
# print(ppv_twentythree_SrRNA_drtr_bpp)

#F1
f_1_twentythree_SrRNA_fold_mfe = f_1(tp_twentythree_SrRNA_fold_mfe_sum, fp_twentythree_SrRNA_fold_mfe_sum, fn_twentythree_SrRNA_fold_mfe_sum)
# print(f_1_twentythree_SrRNA_fold_mfe)

f_1_twentythree_SrRNA_fold_bpp = f_1(tp_twentythree_SrRNA_fold_bpp_sum, fp_twentythree_SrRNA_fold_bpp_sum, fn_twentythree_SrRNA_fold_bpp_sum)
# print(f_1_twentythree_SrRNA_fold_bpp)

f_1_twentythree_SrRNA_drtr_mfe = f_1(tp_twentythree_SrRNA_drtr_mfe_sum, fp_twentythree_SrRNA_drtr_mfe_sum, fn_twentythree_SrRNA_drtr_mfe_sum)
# print(f_1_twentythree_SrRNA_drtr_mfe)

f_1_twentythree_SrRNA_drtr_prob = f_1(tp_twentythree_SrRNA_drtr_prob_sum, fp_twentythree_SrRNA_drtr_prob_sum, fn_twentythree_SrRNA_drtr_prob_sum)
# print(f_1_twentythree_SrRNA_drtr_prob)

f_1_twentythree_SrRNA_drtr_bpp = f_1(tp_twentythree_SrRNA_drtr_bpp_sum, fp_twentythree_SrRNA_drtr_bpp_sum, fn_twentythree_SrRNA_drtr_bpp_sum)
# print(f_1_twentythree_SrRNA_drtr_bpp)

# #MCC
mcc_twentythree_SrRNA_fold_mfe = mcc(tp_twentythree_SrRNA_fold_mfe_sum, fp_twentythree_SrRNA_fold_mfe_sum, tn_twentythree_SrRNA_fold_mfe_sum, fn_twentythree_SrRNA_fold_mfe_sum)
# print(mcc_twentythree_SrRNA_fold_mfe)

mcc_twentythree_SrRNA_fold_bpp = mcc(tp_twentythree_SrRNA_fold_bpp_sum, fp_twentythree_SrRNA_fold_bpp_sum, tn_twentythree_SrRNA_fold_bpp_sum, fn_twentythree_SrRNA_fold_bpp_sum)
# print(mcc_twentythree_SrRNA_fold_bpp)

mcc_twentythree_SrRNA_drtr_mfe = mcc(tp_twentythree_SrRNA_drtr_mfe_sum, fp_twentythree_SrRNA_drtr_mfe_sum, tn_twentythree_SrRNA_drtr_mfe_sum, fn_twentythree_SrRNA_drtr_mfe_sum)
# print(mcc_twentythree_SrRNA_drtr_mfe)

mcc_twentythree_SrRNA_drtr_prob = mcc(tp_twentythree_SrRNA_drtr_prob_sum, fp_twentythree_SrRNA_drtr_prob_sum, tn_twentythree_SrRNA_drtr_prob_sum, fn_twentythree_SrRNA_drtr_prob_sum)
# print(mcc_twentythree_SrRNA_drtr_prob)

mcc_twentythree_SrRNA_drtr_bpp = mcc(tp_twentythree_SrRNA_drtr_bpp_sum, fp_twentythree_SrRNA_drtr_bpp_sum, tn_twentythree_SrRNA_drtr_bpp_sum, fn_twentythree_SrRNA_drtr_bpp_sum)
# print(mcc_twentythree_SrRNA_drtr_bpp)


# 5SrRNA
five_SrRNA=family_data.loc[family_data.family=='5SrRNA']

tp_five_SrRNA_fold_mfe = (five_SrRNA[['RNAfold_mfe_tp']])
tp_five_SrRNA_fold_mfe_sum = boot(tp_five_SrRNA_fold_mfe)

fp_five_SrRNA_fold_mfe = (five_SrRNA[['RNAfold_mfe_fp']])
fp_five_SrRNA_fold_mfe_sum = boot(fp_five_SrRNA_fold_mfe)

tn_five_SrRNA_fold_mfe = (five_SrRNA[['RNAfold_mfe_tn']])
tn_five_SrRNA_fold_mfe_sum = boot(tn_five_SrRNA_fold_mfe)

fn_five_SrRNA_fold_mfe = (five_SrRNA[['RNAfold_mfe_fn']])
fn_five_SrRNA_fold_mfe_sum = boot(fn_five_SrRNA_fold_mfe)


tp_five_SrRNA_fold_bpp = (five_SrRNA[['RNAfold_bpp_tp']])
tp_five_SrRNA_fold_bpp_sum = boot(tp_five_SrRNA_fold_bpp)

fp_five_SrRNA_fold_bpp = (five_SrRNA[['RNAfold_bpp_fp']])
fp_five_SrRNA_fold_bpp_sum = boot(fp_five_SrRNA_fold_bpp)

tn_five_SrRNA_fold_bpp = (five_SrRNA[['RNAfold_bpp_tn']])
tn_five_SrRNA_fold_bpp_sum = boot(tn_five_SrRNA_fold_bpp)

fn_five_SrRNA_fold_bpp = (five_SrRNA[['RNAfold_bpp_fn']])
fn_five_SrRNA_fold_bpp_sum = boot(fn_five_SrRNA_fold_bpp)


tp_five_SrRNA_drtr_mfe = (five_SrRNA[['drtr_tp_mfe']])
tp_five_SrRNA_drtr_mfe_sum = boot(tp_five_SrRNA_drtr_mfe)

fp_five_SrRNA_drtr_mfe = (five_SrRNA[['drtr_fp_mfe']])
fp_five_SrRNA_drtr_mfe_sum = boot(fp_five_SrRNA_drtr_mfe)

tn_five_SrRNA_drtr_mfe = (five_SrRNA[['drtr_tn_mfe']])
tn_five_SrRNA_drtr_mfe_sum = boot(tn_five_SrRNA_drtr_mfe)

fn_five_SrRNA_drtr_mfe = (five_SrRNA[['drtr_fn_mfe']])
fn_five_SrRNA_drtr_mfe_sum = boot(fn_five_SrRNA_drtr_mfe)


tp_five_SrRNA_drtr_prob = (five_SrRNA[['drtr_tp_prob']])
tp_five_SrRNA_drtr_prob_sum = boot(tp_five_SrRNA_drtr_prob)

fp_five_SrRNA_drtr_prob = (five_SrRNA[['drtr_fp_prob']])
fp_five_SrRNA_drtr_prob_sum = boot(fp_five_SrRNA_drtr_prob)

tn_five_SrRNA_drtr_prob = (five_SrRNA[['drtr_tn_prob']])
tn_five_SrRNA_drtr_prob_sum = boot(tn_five_SrRNA_drtr_prob)

fn_five_SrRNA_drtr_prob = (five_SrRNA[['drtr_fn_prob']])
fn_five_SrRNA_drtr_prob_sum = boot(fn_five_SrRNA_drtr_prob)


tp_five_SrRNA_drtr_bpp = (five_SrRNA[['drtr_tp']])
tp_five_SrRNA_drtr_bpp_sum = boot(tp_five_SrRNA_drtr_bpp)

fp_five_SrRNA_drtr_bpp = (five_SrRNA[['drtr_fp']])
fp_five_SrRNA_drtr_bpp_sum = boot(fp_five_SrRNA_drtr_bpp)

tn_five_SrRNA_drtr_bpp = (five_SrRNA[['drtr_tn']])
tn_five_SrRNA_drtr_bpp_sum = boot(tn_five_SrRNA_drtr_bpp)

fn_five_SrRNA_drtr_bpp = (five_SrRNA[['drtr_fn']])
fn_five_SrRNA_drtr_bpp_sum = boot(fn_five_SrRNA_drtr_bpp)

#TPR
tpr_five_SrRNA_fold_mfe = tpr(tp_five_SrRNA_fold_mfe_sum, fn_five_SrRNA_fold_mfe_sum)
# print(tpr_five_SrRNA_fold_mfe)

tpr_five_SrRNA_fold_bpp = tpr(tp_five_SrRNA_fold_bpp_sum, fn_five_SrRNA_fold_bpp_sum)
# print(tpr_five_SrRNA_fold_bpp)

tpr_five_SrRNA_drtr_mfe = tpr(tp_five_SrRNA_drtr_mfe_sum, fn_five_SrRNA_drtr_mfe_sum)
# print(tpr_five_SrRNA_drtr_mfe)

tpr_five_SrRNA_drtr_prob = tpr(tp_five_SrRNA_drtr_prob_sum, fn_five_SrRNA_drtr_prob_sum)
# print(tpr_five_SrRNA_drtr_prob)

tpr_five_SrRNA_drtr_bpp = tpr(tp_five_SrRNA_drtr_bpp_sum, fn_five_SrRNA_drtr_bpp_sum)
# print(tpr_five_SrRNA_drtr_bpp)


#PPV
ppv_five_SrRNA_fold_mfe = ppv(tp_five_SrRNA_fold_mfe_sum, fn_five_SrRNA_fold_mfe_sum)
# print(ppv_five_SrRNA_fold_mfe)

ppv_five_SrRNA_fold_bpp = ppv(tp_five_SrRNA_fold_bpp_sum, fn_five_SrRNA_fold_bpp_sum)
# print(ppv_five_SrRNA_fold_bpp)

ppv_five_SrRNA_drtr_mfe = ppv(tp_five_SrRNA_drtr_mfe_sum, fn_five_SrRNA_drtr_mfe_sum)
# print(ppv_five_SrRNA_drtr_mfe)

ppv_five_SrRNA_drtr_prob = ppv(tp_five_SrRNA_drtr_prob_sum, fn_five_SrRNA_drtr_prob_sum)
# print(ppv_five_SrRNA_drtr_prob)

ppv_five_SrRNA_drtr_bpp = ppv(tp_five_SrRNA_drtr_bpp_sum, fn_five_SrRNA_drtr_bpp_sum)
# print(ppv_five_SrRNA_drtr_bpp)

#F1
f_1_five_SrRNA_fold_mfe = f_1(tp_five_SrRNA_fold_mfe_sum, fp_five_SrRNA_fold_mfe_sum, fn_five_SrRNA_fold_mfe_sum)
# print(f_1_five_SrRNA_fold_mfe)

f_1_five_SrRNA_fold_bpp = f_1(tp_five_SrRNA_fold_bpp_sum, fp_five_SrRNA_fold_bpp_sum, fn_five_SrRNA_fold_bpp_sum)
# print(f_1_five_SrRNA_fold_bpp)

f_1_five_SrRNA_drtr_mfe = f_1(tp_five_SrRNA_drtr_mfe_sum, fp_five_SrRNA_drtr_mfe_sum, fn_five_SrRNA_drtr_mfe_sum)
# print(f_1_five_SrRNA_drtr_mfe)

f_1_five_SrRNA_drtr_prob = f_1(tp_five_SrRNA_drtr_prob_sum, fp_five_SrRNA_drtr_prob_sum, fn_five_SrRNA_drtr_prob_sum)
# print(f_1_five_SrRNA_drtr_prob)

f_1_five_SrRNA_drtr_bpp = f_1(tp_five_SrRNA_drtr_bpp_sum, fp_five_SrRNA_drtr_bpp_sum, fn_five_SrRNA_drtr_bpp_sum)
# print(f_1_five_SrRNA_drtr_bpp)

# #MCC
mcc_five_SrRNA_fold_mfe = mcc(tp_five_SrRNA_fold_mfe_sum, fp_five_SrRNA_fold_mfe_sum, tn_five_SrRNA_fold_mfe_sum, fn_five_SrRNA_fold_mfe_sum)
# print(mcc_five_SrRNA_fold_mfe)

mcc_five_SrRNA_fold_bpp = mcc(tp_five_SrRNA_fold_bpp_sum, fp_five_SrRNA_fold_bpp_sum, tn_five_SrRNA_fold_bpp_sum, fn_five_SrRNA_fold_bpp_sum)
# print(mcc_five_SrRNA_fold_bpp)

mcc_five_SrRNA_drtr_mfe = mcc(tp_five_SrRNA_drtr_mfe_sum, fp_five_SrRNA_drtr_mfe_sum, tn_five_SrRNA_drtr_mfe_sum, fn_five_SrRNA_drtr_mfe_sum)
# print(mcc_five_SrRNA_drtr_mfe)

mcc_five_SrRNA_drtr_prob = mcc(tp_five_SrRNA_drtr_prob_sum, fp_five_SrRNA_drtr_prob_sum, tn_five_SrRNA_drtr_prob_sum, fn_five_SrRNA_drtr_prob_sum)
# print(mcc_five_SrRNA_drtr_prob)

mcc_five_SrRNA_drtr_bpp = mcc(tp_five_SrRNA_drtr_bpp_sum, fp_five_SrRNA_drtr_bpp_sum, tn_five_SrRNA_drtr_bpp_sum, fn_five_SrRNA_drtr_bpp_sum)
# print(mcc_five_SrRNA_drtr_bpp)



#Cili.Telo. RNA
cili_telo_RNA=family_data.loc[family_data.family=='Cili.Telo. RNA']

tp_cili_telo_RNA_fold_mfe = (cili_telo_RNA[['RNAfold_mfe_tp']])
tp_cili_telo_RNA_fold_mfe_sum = boot(tp_cili_telo_RNA_fold_mfe)

fp_cili_telo_RNA_fold_mfe = (cili_telo_RNA[['RNAfold_mfe_fp']])
fp_cili_telo_RNA_fold_mfe_sum = boot(fp_cili_telo_RNA_fold_mfe)

tn_cili_telo_RNA_fold_mfe = (cili_telo_RNA[['RNAfold_mfe_tn']])
tn_cili_telo_RNA_fold_mfe_sum = boot(tn_cili_telo_RNA_fold_mfe)

fn_cili_telo_RNA_fold_mfe = (cili_telo_RNA[['RNAfold_mfe_fn']])
fn_cili_telo_RNA_fold_mfe_sum = boot(fn_cili_telo_RNA_fold_mfe)


tp_cili_telo_RNA_fold_bpp = (cili_telo_RNA[['RNAfold_bpp_tp']])
tp_cili_telo_RNA_fold_bpp_sum = boot(tp_cili_telo_RNA_fold_bpp)

fp_cili_telo_RNA_fold_bpp = (cili_telo_RNA[['RNAfold_bpp_fp']])
fp_cili_telo_RNA_fold_bpp_sum = boot(fp_cili_telo_RNA_fold_bpp)

tn_cili_telo_RNA_fold_bpp = (cili_telo_RNA[['RNAfold_bpp_tn']])
tn_cili_telo_RNA_fold_bpp_sum = boot(tn_cili_telo_RNA_fold_bpp)

fn_cili_telo_RNA_fold_bpp = (cili_telo_RNA[['RNAfold_bpp_fn']])
fn_cili_telo_RNA_fold_bpp_sum = boot(fn_cili_telo_RNA_fold_bpp)


tp_cili_telo_RNA_drtr_mfe = (cili_telo_RNA[['drtr_tp_mfe']])
tp_cili_telo_RNA_drtr_mfe_sum = boot(tp_cili_telo_RNA_drtr_mfe)

fp_cili_telo_RNA_drtr_mfe = (cili_telo_RNA[['drtr_fp_mfe']])
fp_cili_telo_RNA_drtr_mfe_sum = boot(fp_cili_telo_RNA_drtr_mfe)

tn_cili_telo_RNA_drtr_mfe = (cili_telo_RNA[['drtr_tn_mfe']])
tn_cili_telo_RNA_drtr_mfe_sum = boot(tn_cili_telo_RNA_drtr_mfe)

fn_cili_telo_RNA_drtr_mfe = (cili_telo_RNA[['drtr_fn_mfe']])
fn_cili_telo_RNA_drtr_mfe_sum = boot(fn_cili_telo_RNA_drtr_mfe)


tp_cili_telo_RNA_drtr_prob = (cili_telo_RNA[['drtr_tp_prob']])
tp_cili_telo_RNA_drtr_prob_sum = boot(tp_cili_telo_RNA_drtr_prob)

fp_cili_telo_RNA_drtr_prob = (cili_telo_RNA[['drtr_fp_prob']])
fp_cili_telo_RNA_drtr_prob_sum = boot(fp_cili_telo_RNA_drtr_prob)

tn_cili_telo_RNA_drtr_prob = (cili_telo_RNA[['drtr_tn_prob']])
tn_cili_telo_RNA_drtr_prob_sum = boot(tn_cili_telo_RNA_drtr_prob)

fn_cili_telo_RNA_drtr_prob = (cili_telo_RNA[['drtr_fn_prob']])
fn_cili_telo_RNA_drtr_prob_sum = boot(fn_cili_telo_RNA_drtr_prob)


tp_cili_telo_RNA_drtr_bpp = (cili_telo_RNA[['drtr_tp']])
tp_cili_telo_RNA_drtr_bpp_sum = boot(tp_cili_telo_RNA_drtr_bpp)

fp_cili_telo_RNA_drtr_bpp = (cili_telo_RNA[['drtr_fp']])
fp_cili_telo_RNA_drtr_bpp_sum = boot(fp_cili_telo_RNA_drtr_bpp)

tn_cili_telo_RNA_drtr_bpp = (cili_telo_RNA[['drtr_tn']])
tn_cili_telo_RNA_drtr_bpp_sum = boot(tn_cili_telo_RNA_drtr_bpp)

fn_cili_telo_RNA_drtr_bpp = (cili_telo_RNA[['drtr_fn']])
fn_cili_telo_RNA_drtr_bpp_sum = boot(fn_cili_telo_RNA_drtr_bpp)

#TPR
tpr_cili_telo_RNA_fold_mfe = tpr(tp_cili_telo_RNA_fold_mfe_sum, fn_cili_telo_RNA_fold_mfe_sum)
# print(tpr_cili_telo_RNA_fold_mfe)

tpr_cili_telo_RNA_fold_bpp = tpr(tp_cili_telo_RNA_fold_bpp_sum, fn_cili_telo_RNA_fold_bpp_sum)
# print(tpr_cili_telo_RNA_fold_bpp)

tpr_cili_telo_RNA_drtr_mfe = tpr(tp_cili_telo_RNA_drtr_mfe_sum, fn_cili_telo_RNA_drtr_mfe_sum)
# print(tpr_cili_telo_RNA_drtr_mfe)

tpr_cili_telo_RNA_drtr_prob = tpr(tp_cili_telo_RNA_drtr_prob_sum, fn_cili_telo_RNA_drtr_prob_sum)
# print(tpr_cili_telo_RNA_drtr_prob)

tpr_cili_telo_RNA_drtr_bpp = tpr(tp_cili_telo_RNA_drtr_bpp_sum, fn_cili_telo_RNA_drtr_bpp_sum)
# print(tpr_cili_telo_RNA_drtr_bpp)


#PPV
ppv_cili_telo_RNA_fold_mfe = ppv(tp_cili_telo_RNA_fold_mfe_sum, fn_cili_telo_RNA_fold_mfe_sum)
# print(ppv_cili_telo_RNA_fold_mfe)

ppv_cili_telo_RNA_fold_bpp = ppv(tp_cili_telo_RNA_fold_bpp_sum, fn_cili_telo_RNA_fold_bpp_sum)
# print(ppv_cili_telo_RNA_fold_bpp)

ppv_cili_telo_RNA_drtr_mfe = ppv(tp_cili_telo_RNA_drtr_mfe_sum, fn_cili_telo_RNA_drtr_mfe_sum)
# print(ppv_cili_telo_RNA_drtr_mfe)

ppv_cili_telo_RNA_drtr_prob = ppv(tp_cili_telo_RNA_drtr_prob_sum, fn_cili_telo_RNA_drtr_prob_sum)
# print(ppv_cili_telo_RNA_drtr_prob)

ppv_cili_telo_RNA_drtr_bpp = ppv(tp_cili_telo_RNA_drtr_bpp_sum, fn_cili_telo_RNA_drtr_bpp_sum)
# print(ppv_cili_telo_RNA_drtr_bpp)

#F1
f_1_cili_telo_RNA_fold_mfe = f_1(tp_cili_telo_RNA_fold_mfe_sum, fp_cili_telo_RNA_fold_mfe_sum, fn_cili_telo_RNA_fold_mfe_sum)
# print(f_1_cili_telo_RNA_fold_mfe)

f_1_cili_telo_RNA_fold_bpp = f_1(tp_cili_telo_RNA_fold_bpp_sum, fp_cili_telo_RNA_fold_bpp_sum, fn_cili_telo_RNA_fold_bpp_sum)
# print(f_1_cili_telo_RNA_fold_bpp)

f_1_cili_telo_RNA_drtr_mfe = f_1(tp_cili_telo_RNA_drtr_mfe_sum, fp_cili_telo_RNA_drtr_mfe_sum, fn_cili_telo_RNA_drtr_mfe_sum)
# print(f_1_cili_telo_RNA_drtr_mfe)

f_1_cili_telo_RNA_drtr_prob = f_1(tp_cili_telo_RNA_drtr_prob_sum, fp_cili_telo_RNA_drtr_prob_sum, fn_cili_telo_RNA_drtr_prob_sum)
# print(f_1_cili_telo_RNA_drtr_prob)

f_1_cili_telo_RNA_drtr_bpp = f_1(tp_cili_telo_RNA_drtr_bpp_sum, fp_cili_telo_RNA_drtr_bpp_sum, fn_cili_telo_RNA_drtr_bpp_sum)
# print(f_1_cili_telo_RNA_drtr_bpp)

# #MCC
mcc_cili_telo_RNA_fold_mfe = mcc(tp_cili_telo_RNA_fold_mfe_sum, fp_cili_telo_RNA_fold_mfe_sum, tn_cili_telo_RNA_fold_mfe_sum, fn_cili_telo_RNA_fold_mfe_sum)
# print(mcc_cili_telo_RNA_fold_mfe)

mcc_cili_telo_RNA_fold_bpp = mcc(tp_cili_telo_RNA_fold_bpp_sum, fp_cili_telo_RNA_fold_bpp_sum, tn_cili_telo_RNA_fold_bpp_sum, fn_cili_telo_RNA_fold_bpp_sum)
# print(mcc_cili_telo_RNA_fold_bpp)

mcc_cili_telo_RNA_drtr_mfe = mcc(tp_cili_telo_RNA_drtr_mfe_sum, fp_cili_telo_RNA_drtr_mfe_sum, tn_cili_telo_RNA_drtr_mfe_sum, fn_cili_telo_RNA_drtr_mfe_sum)
# print(mcc_cili_telo_RNA_drtr_mfe)

mcc_cili_telo_RNA_drtr_prob = mcc(tp_cili_telo_RNA_drtr_prob_sum, fp_cili_telo_RNA_drtr_prob_sum, tn_cili_telo_RNA_drtr_prob_sum, fn_cili_telo_RNA_drtr_prob_sum)
# print(mcc_cili_telo_RNA_drtr_prob)

mcc_cili_telo_RNA_drtr_bpp = mcc(tp_cili_telo_RNA_drtr_bpp_sum, fp_cili_telo_RNA_drtr_bpp_sum, tn_cili_telo_RNA_drtr_bpp_sum, fn_cili_telo_RNA_drtr_bpp_sum)
# print(mcc_cili_telo_RNA_drtr_bpp)

#Cis-reg.element
cis_regulatory_element=family_data.loc[family_data.family=='Cis-reg.element']

tp_cis_regulatory_element_fold_mfe = (cis_regulatory_element[['RNAfold_mfe_tp']])
tp_cis_regulatory_element_fold_mfe_sum = boot(tp_cis_regulatory_element_fold_mfe)

fp_cis_regulatory_element_fold_mfe = (cis_regulatory_element[['RNAfold_mfe_fp']])
fp_cis_regulatory_element_fold_mfe_sum = boot(fp_cis_regulatory_element_fold_mfe)

tn_cis_regulatory_element_fold_mfe = (cis_regulatory_element[['RNAfold_mfe_tn']])
tn_cis_regulatory_element_fold_mfe_sum = boot(tn_cis_regulatory_element_fold_mfe)

fn_cis_regulatory_element_fold_mfe = (cis_regulatory_element[['RNAfold_mfe_fn']])
fn_cis_regulatory_element_fold_mfe_sum = boot(fn_cis_regulatory_element_fold_mfe)


tp_cis_regulatory_element_fold_bpp = (cis_regulatory_element[['RNAfold_bpp_tp']])
tp_cis_regulatory_element_fold_bpp_sum = boot(tp_cis_regulatory_element_fold_bpp)

fp_cis_regulatory_element_fold_bpp = (cis_regulatory_element[['RNAfold_bpp_fp']])
fp_cis_regulatory_element_fold_bpp_sum = boot(fp_cis_regulatory_element_fold_bpp)

tn_cis_regulatory_element_fold_bpp = (cis_regulatory_element[['RNAfold_bpp_tn']])
tn_cis_regulatory_element_fold_bpp_sum = boot(tn_cis_regulatory_element_fold_bpp)

fn_cis_regulatory_element_fold_bpp = (cis_regulatory_element[['RNAfold_bpp_fn']])
fn_cis_regulatory_element_fold_bpp_sum = boot(fn_cis_regulatory_element_fold_bpp)


tp_cis_regulatory_element_drtr_mfe = (cis_regulatory_element[['drtr_tp_mfe']])
tp_cis_regulatory_element_drtr_mfe_sum = boot(tp_cis_regulatory_element_drtr_mfe)

fp_cis_regulatory_element_drtr_mfe = (cis_regulatory_element[['drtr_fp_mfe']])
fp_cis_regulatory_element_drtr_mfe_sum = boot(fp_cis_regulatory_element_drtr_mfe)

tn_cis_regulatory_element_drtr_mfe = (cis_regulatory_element[['drtr_tn_mfe']])
tn_cis_regulatory_element_drtr_mfe_sum = boot(tn_cis_regulatory_element_drtr_mfe)

fn_cis_regulatory_element_drtr_mfe = (cis_regulatory_element[['drtr_fn_mfe']])
fn_cis_regulatory_element_drtr_mfe_sum = boot(fn_cis_regulatory_element_drtr_mfe)


tp_cis_regulatory_element_drtr_prob = (cis_regulatory_element[['drtr_tp_prob']])
tp_cis_regulatory_element_drtr_prob_sum = boot(tp_cis_regulatory_element_drtr_prob)

fp_cis_regulatory_element_drtr_prob = (cis_regulatory_element[['drtr_fp_prob']])
fp_cis_regulatory_element_drtr_prob_sum = boot(fp_cis_regulatory_element_drtr_prob)

tn_cis_regulatory_element_drtr_prob = (cis_regulatory_element[['drtr_tn_prob']])
tn_cis_regulatory_element_drtr_prob_sum = boot(tn_cis_regulatory_element_drtr_prob)

fn_cis_regulatory_element_drtr_prob = (cis_regulatory_element[['drtr_fn_prob']])
fn_cis_regulatory_element_drtr_prob_sum = boot(fn_cis_regulatory_element_drtr_prob)


tp_cis_regulatory_element_drtr_bpp = (cis_regulatory_element[['drtr_tp']])
tp_cis_regulatory_element_drtr_bpp_sum = boot(tp_cis_regulatory_element_drtr_bpp)

fp_cis_regulatory_element_drtr_bpp = (cis_regulatory_element[['drtr_fp']])
fp_cis_regulatory_element_drtr_bpp_sum = boot(fp_cis_regulatory_element_drtr_bpp)

tn_cis_regulatory_element_drtr_bpp = (cis_regulatory_element[['drtr_tn']])
tn_cis_regulatory_element_drtr_bpp_sum = boot(tn_cis_regulatory_element_drtr_bpp)

fn_cis_regulatory_element_drtr_bpp = (cis_regulatory_element[['drtr_fn']])
fn_cis_regulatory_element_drtr_bpp_sum = boot(fn_cis_regulatory_element_drtr_bpp)

#TPR
tpr_cis_regulatory_element_fold_mfe = tpr(tp_cis_regulatory_element_fold_mfe_sum, fn_cis_regulatory_element_fold_mfe_sum)
# print(tpr_cis_regulatory_element_fold_mfe)

tpr_cis_regulatory_element_fold_bpp = tpr(tp_cis_regulatory_element_fold_bpp_sum, fn_cis_regulatory_element_fold_bpp_sum)
# print(tpr_cis_regulatory_element_fold_bpp)

tpr_cis_regulatory_element_drtr_mfe = tpr(tp_cis_regulatory_element_drtr_mfe_sum, fn_cis_regulatory_element_drtr_mfe_sum)
# print(tpr_cis_regulatory_element_drtr_mfe)

tpr_cis_regulatory_element_drtr_prob = tpr(tp_cis_regulatory_element_drtr_prob_sum, fn_cis_regulatory_element_drtr_prob_sum)
# print(tpr_cis_regulatory_element_drtr_prob)

tpr_cis_regulatory_element_drtr_bpp = tpr(tp_cis_regulatory_element_drtr_bpp_sum, fn_cis_regulatory_element_drtr_bpp_sum)
# print(tpr_cis_regulatory_element_drtr_bpp)


#PPV
ppv_cis_regulatory_element_fold_mfe = ppv(tp_cis_regulatory_element_fold_mfe_sum, fn_cis_regulatory_element_fold_mfe_sum)
# print(ppv_cis_regulatory_element_fold_mfe)

ppv_cis_regulatory_element_fold_bpp = ppv(tp_cis_regulatory_element_fold_bpp_sum, fn_cis_regulatory_element_fold_bpp_sum)
# print(ppv_cis_regulatory_element_fold_bpp)

ppv_cis_regulatory_element_drtr_mfe = ppv(tp_cis_regulatory_element_drtr_mfe_sum, fn_cis_regulatory_element_drtr_mfe_sum)
# print(ppv_cis_regulatory_element_drtr_mfe)

ppv_cis_regulatory_element_drtr_prob = ppv(tp_cis_regulatory_element_drtr_prob_sum, fn_cis_regulatory_element_drtr_prob_sum)
# print(ppv_cis_regulatory_element_drtr_prob)

ppv_cis_regulatory_element_drtr_bpp = ppv(tp_cis_regulatory_element_drtr_bpp_sum, fn_cis_regulatory_element_drtr_bpp_sum)
# print(ppv_cis_regulatory_element_drtr_bpp)

#F1
f_1_cis_regulatory_element_fold_mfe = f_1(tp_cis_regulatory_element_fold_mfe_sum, fp_cis_regulatory_element_fold_mfe_sum, fn_cis_regulatory_element_fold_mfe_sum)
# print(f_1_cis_regulatory_element_fold_mfe)

f_1_cis_regulatory_element_fold_bpp = f_1(tp_cis_regulatory_element_fold_bpp_sum, fp_cis_regulatory_element_fold_bpp_sum, fn_cis_regulatory_element_fold_bpp_sum)
# print(f_1_cis_regulatory_element_fold_bpp)

f_1_cis_regulatory_element_drtr_mfe = f_1(tp_cis_regulatory_element_drtr_mfe_sum, fp_cis_regulatory_element_drtr_mfe_sum, fn_cis_regulatory_element_drtr_mfe_sum)
# print(f_1_cis_regulatory_element_drtr_mfe)

f_1_cis_regulatory_element_drtr_prob = f_1(tp_cis_regulatory_element_drtr_prob_sum, fp_cis_regulatory_element_drtr_prob_sum, fn_cis_regulatory_element_drtr_prob_sum)
# print(f_1_cis_regulatory_element_drtr_prob)

f_1_cis_regulatory_element_drtr_bpp = f_1(tp_cis_regulatory_element_drtr_bpp_sum, fp_cis_regulatory_element_drtr_bpp_sum, fn_cis_regulatory_element_drtr_bpp_sum)
# print(f_1_cis_regulatory_element_drtr_bpp)

# #MCC
mcc_cis_regulatory_element_fold_mfe = mcc(tp_cis_regulatory_element_fold_mfe_sum, fp_cis_regulatory_element_fold_mfe_sum, tn_cis_regulatory_element_fold_mfe_sum, fn_cis_regulatory_element_fold_mfe_sum)
# print(mcc_cis_regulatory_element_fold_mfe)

mcc_cis_regulatory_element_fold_bpp = mcc(tp_cis_regulatory_element_fold_bpp_sum, fp_cis_regulatory_element_fold_bpp_sum, tn_cis_regulatory_element_fold_bpp_sum, fn_cis_regulatory_element_fold_bpp_sum)
# print(mcc_cis_regulatory_element_fold_bpp)

mcc_cis_regulatory_element_drtr_mfe = mcc(tp_cis_regulatory_element_drtr_mfe_sum, fp_cis_regulatory_element_drtr_mfe_sum, tn_cis_regulatory_element_drtr_mfe_sum, fn_cis_regulatory_element_drtr_mfe_sum)
# print(mcc_cis_regulatory_element_drtr_mfe)

mcc_cis_regulatory_element_drtr_prob = mcc(tp_cis_regulatory_element_drtr_prob_sum, fp_cis_regulatory_element_drtr_prob_sum, tn_cis_regulatory_element_drtr_prob_sum, fn_cis_regulatory_element_drtr_prob_sum)
# print(mcc_cis_regulatory_element_drtr_prob)

mcc_cis_regulatory_element_drtr_bpp = mcc(tp_cis_regulatory_element_drtr_bpp_sum, fp_cis_regulatory_element_drtr_bpp_sum, tn_cis_regulatory_element_drtr_bpp_sum, fn_cis_regulatory_element_drtr_bpp_sum)
# print(mcc_cis_regulatory_element_drtr_bpp)

#GIIIntron
GIIIntron=family_data.loc[family_data.family=='GIIIntron']
tp_GIIIntron_fold_mfe = (GIIIntron[['RNAfold_mfe_tp']])
tp_GIIIntron_fold_mfe_sum = boot(tp_GIIIntron_fold_mfe)

fp_GIIIntron_fold_mfe = (GIIIntron[['RNAfold_mfe_fp']])
fp_GIIIntron_fold_mfe_sum = boot(fp_GIIIntron_fold_mfe)

tn_GIIIntron_fold_mfe = (GIIIntron[['RNAfold_mfe_tn']])
tn_GIIIntron_fold_mfe_sum = boot(tn_GIIIntron_fold_mfe)

fn_GIIIntron_fold_mfe = (GIIIntron[['RNAfold_mfe_fn']])
fn_GIIIntron_fold_mfe_sum = boot(fn_GIIIntron_fold_mfe)


tp_GIIIntron_fold_bpp = (GIIIntron[['RNAfold_bpp_tp']])
tp_GIIIntron_fold_bpp_sum = boot(tp_GIIIntron_fold_bpp)

fp_GIIIntron_fold_bpp = (GIIIntron[['RNAfold_bpp_fp']])
fp_GIIIntron_fold_bpp_sum = boot(fp_GIIIntron_fold_bpp)

tn_GIIIntron_fold_bpp = (GIIIntron[['RNAfold_bpp_tn']])
tn_GIIIntron_fold_bpp_sum = boot(tn_GIIIntron_fold_bpp)

fn_GIIIntron_fold_bpp = (GIIIntron[['RNAfold_bpp_fn']])
fn_GIIIntron_fold_bpp_sum = boot(fn_GIIIntron_fold_bpp)


tp_GIIIntron_drtr_mfe = (GIIIntron[['drtr_tp_mfe']])
tp_GIIIntron_drtr_mfe_sum = boot(tp_GIIIntron_drtr_mfe)

fp_GIIIntron_drtr_mfe = (GIIIntron[['drtr_fp_mfe']])
fp_GIIIntron_drtr_mfe_sum = boot(fp_GIIIntron_drtr_mfe)

tn_GIIIntron_drtr_mfe = (GIIIntron[['drtr_tn_mfe']])
tn_GIIIntron_drtr_mfe_sum = boot(tn_GIIIntron_drtr_mfe)

fn_GIIIntron_drtr_mfe = (GIIIntron[['drtr_fn_mfe']])
fn_GIIIntron_drtr_mfe_sum = boot(fn_GIIIntron_drtr_mfe)


tp_GIIIntron_drtr_prob = (GIIIntron[['drtr_tp_prob']])
tp_GIIIntron_drtr_prob_sum = boot(tp_GIIIntron_drtr_prob)

fp_GIIIntron_drtr_prob = (GIIIntron[['drtr_fp_prob']])
fp_GIIIntron_drtr_prob_sum = boot(fp_GIIIntron_drtr_prob)

tn_GIIIntron_drtr_prob = (GIIIntron[['drtr_tn_prob']])
tn_GIIIntron_drtr_prob_sum = boot(tn_GIIIntron_drtr_prob)

fn_GIIIntron_drtr_prob = (GIIIntron[['drtr_fn_prob']])
fn_GIIIntron_drtr_prob_sum = boot(fn_GIIIntron_drtr_prob)


tp_GIIIntron_drtr_bpp = (GIIIntron[['drtr_tp']])
tp_GIIIntron_drtr_bpp_sum = boot(tp_GIIIntron_drtr_bpp)

fp_GIIIntron_drtr_bpp = (GIIIntron[['drtr_fp']])
fp_GIIIntron_drtr_bpp_sum = boot(fp_GIIIntron_drtr_bpp)

tn_GIIIntron_drtr_bpp = (GIIIntron[['drtr_tn']])
tn_GIIIntron_drtr_bpp_sum = boot(tn_GIIIntron_drtr_bpp)

fn_GIIIntron_drtr_bpp = (GIIIntron[['drtr_fn']])
fn_GIIIntron_drtr_bpp_sum = boot(fn_GIIIntron_drtr_bpp)

#TPR
tpr_GIIIntron_fold_mfe = tpr(tp_GIIIntron_fold_mfe_sum, fn_GIIIntron_fold_mfe_sum)
# print(tpr_GIIIntron_fold_mfe)

tpr_GIIIntron_fold_bpp = tpr(tp_GIIIntron_fold_bpp_sum, fn_GIIIntron_fold_bpp_sum)
# print(tpr_GIIIntron_fold_bpp)

tpr_GIIIntron_drtr_mfe = tpr(tp_GIIIntron_drtr_mfe_sum, fn_GIIIntron_drtr_mfe_sum)
# print(tpr_GIIIntron_drtr_mfe)

tpr_GIIIntron_drtr_prob = tpr(tp_GIIIntron_drtr_prob_sum, fn_GIIIntron_drtr_prob_sum)
# print(tpr_GIIIntron_drtr_prob)

tpr_GIIIntron_drtr_bpp = tpr(tp_GIIIntron_drtr_bpp_sum, fn_GIIIntron_drtr_bpp_sum)
# print(tpr_GIIIntron_drtr_bpp)


#PPV
ppv_GIIIntron_fold_mfe = ppv(tp_GIIIntron_fold_mfe_sum, fn_GIIIntron_fold_mfe_sum)
# print(ppv_GIIIntron_fold_mfe)

ppv_GIIIntron_fold_bpp = ppv(tp_GIIIntron_fold_bpp_sum, fn_GIIIntron_fold_bpp_sum)
# print(ppv_GIIIntron_fold_bpp)

ppv_GIIIntron_drtr_mfe = ppv(tp_GIIIntron_drtr_mfe_sum, fn_GIIIntron_drtr_mfe_sum)
# print(ppv_GIIIntron_drtr_mfe)

ppv_GIIIntron_drtr_prob = ppv(tp_GIIIntron_drtr_prob_sum, fn_GIIIntron_drtr_prob_sum)
# print(ppv_GIIIntron_drtr_prob)

ppv_GIIIntron_drtr_bpp = ppv(tp_GIIIntron_drtr_bpp_sum, fn_GIIIntron_drtr_bpp_sum)
# print(ppv_GIIIntron_drtr_bpp)

#F1
f_1_GIIIntron_fold_mfe = f_1(tp_GIIIntron_fold_mfe_sum, fp_GIIIntron_fold_mfe_sum, fn_GIIIntron_fold_mfe_sum)
# print(f_1_GIIIntron_fold_mfe)

f_1_GIIIntron_fold_bpp = f_1(tp_GIIIntron_fold_bpp_sum, fp_GIIIntron_fold_bpp_sum, fn_GIIIntron_fold_bpp_sum)
# print(f_1_GIIIntron_fold_bpp)

f_1_GIIIntron_drtr_mfe = f_1(tp_GIIIntron_drtr_mfe_sum, fp_GIIIntron_drtr_mfe_sum, fn_GIIIntron_drtr_mfe_sum)
# print(f_1_GIIIntron_drtr_mfe)

f_1_GIIIntron_drtr_prob = f_1(tp_GIIIntron_drtr_prob_sum, fp_GIIIntron_drtr_prob_sum, fn_GIIIntron_drtr_prob_sum)
# print(f_1_GIIIntron_drtr_prob)

f_1_GIIIntron_drtr_bpp = f_1(tp_GIIIntron_drtr_bpp_sum, fp_GIIIntron_drtr_bpp_sum, fn_GIIIntron_drtr_bpp_sum)
# print(f_1_GIIIntron_drtr_bpp)

# #MCC
mcc_GIIIntron_fold_mfe = mcc(tp_GIIIntron_fold_mfe_sum, fp_GIIIntron_fold_mfe_sum, tn_GIIIntron_fold_mfe_sum, fn_GIIIntron_fold_mfe_sum)
# print(mcc_GIIIntron_fold_mfe)

mcc_GIIIntron_fold_bpp = mcc(tp_GIIIntron_fold_bpp_sum, fp_GIIIntron_fold_bpp_sum, tn_GIIIntron_fold_bpp_sum, fn_GIIIntron_fold_bpp_sum)
# print(mcc_GIIIntron_fold_bpp)

mcc_GIIIntron_drtr_mfe = mcc(tp_GIIIntron_drtr_mfe_sum, fp_GIIIntron_drtr_mfe_sum, tn_GIIIntron_drtr_mfe_sum, fn_GIIIntron_drtr_mfe_sum)
# print(mcc_GIIIntron_drtr_mfe)

mcc_GIIIntron_drtr_prob = mcc(tp_GIIIntron_drtr_prob_sum, fp_GIIIntron_drtr_prob_sum, tn_GIIIntron_drtr_prob_sum, fn_GIIIntron_drtr_prob_sum)
# print(mcc_GIIIntron_drtr_prob)

mcc_GIIIntron_drtr_bpp = mcc(tp_GIIIntron_drtr_bpp_sum, fp_GIIIntron_drtr_bpp_sum, tn_GIIIntron_drtr_bpp_sum, fn_GIIIntron_drtr_bpp_sum)
# print(mcc_GIIIntron_drtr_bpp)



#GIIntron
GIIntron=family_data.loc[family_data.family=='GIIntron']

tp_GIIntron_fold_mfe = (GIIntron[['RNAfold_mfe_tp']])
tp_GIIntron_fold_mfe_sum = boot(tp_GIIntron_fold_mfe)

fp_GIIntron_fold_mfe = (GIIntron[['RNAfold_mfe_fp']])
fp_GIIntron_fold_mfe_sum = boot(fp_GIIntron_fold_mfe)

tn_GIIntron_fold_mfe = (GIIntron[['RNAfold_mfe_tn']])
tn_GIIntron_fold_mfe_sum = boot(tn_GIIntron_fold_mfe)

fn_GIIntron_fold_mfe = (GIIntron[['RNAfold_mfe_fn']])
fn_GIIntron_fold_mfe_sum = boot(fn_GIIntron_fold_mfe)


tp_GIIntron_fold_bpp = (GIIntron[['RNAfold_bpp_tp']])
tp_GIIntron_fold_bpp_sum = boot(tp_GIIntron_fold_bpp)

fp_GIIntron_fold_bpp = (GIIntron[['RNAfold_bpp_fp']])
fp_GIIntron_fold_bpp_sum = boot(fp_GIIntron_fold_bpp)

tn_GIIntron_fold_bpp = (GIIntron[['RNAfold_bpp_tn']])
tn_GIIntron_fold_bpp_sum = boot(tn_GIIntron_fold_bpp)

fn_GIIntron_fold_bpp = (GIIntron[['RNAfold_bpp_fn']])
fn_GIIntron_fold_bpp_sum = boot(fn_GIIntron_fold_bpp)


tp_GIIntron_drtr_mfe = (GIIntron[['drtr_tp_mfe']])
tp_GIIntron_drtr_mfe_sum = boot(tp_GIIntron_drtr_mfe)

fp_GIIntron_drtr_mfe = (GIIntron[['drtr_fp_mfe']])
fp_GIIntron_drtr_mfe_sum = boot(fp_GIIntron_drtr_mfe)

tn_GIIntron_drtr_mfe = (GIIntron[['drtr_tn_mfe']])
tn_GIIntron_drtr_mfe_sum = boot(tn_GIIntron_drtr_mfe)

fn_GIIntron_drtr_mfe = (GIIntron[['drtr_fn_mfe']])
fn_GIIntron_drtr_mfe_sum = boot(fn_GIIntron_drtr_mfe)


tp_GIIntron_drtr_prob = (GIIntron[['drtr_tp_prob']])
tp_GIIntron_drtr_prob_sum = boot(tp_GIIntron_drtr_prob)

fp_GIIntron_drtr_prob = (GIIntron[['drtr_fp_prob']])
fp_GIIntron_drtr_prob_sum = boot(fp_GIIntron_drtr_prob)

tn_GIIntron_drtr_prob = (GIIntron[['drtr_tn_prob']])
tn_GIIntron_drtr_prob_sum = boot(tn_GIIntron_drtr_prob)

fn_GIIntron_drtr_prob = (GIIntron[['drtr_fn_prob']])
fn_GIIntron_drtr_prob_sum = boot(fn_GIIntron_drtr_prob)


tp_GIIntron_drtr_bpp = (GIIntron[['drtr_tp']])
tp_GIIntron_drtr_bpp_sum = boot(tp_GIIntron_drtr_bpp)

fp_GIIntron_drtr_bpp = (GIIntron[['drtr_fp']])
fp_GIIntron_drtr_bpp_sum = boot(fp_GIIntron_drtr_bpp)

tn_GIIntron_drtr_bpp = (GIIntron[['drtr_tn']])
tn_GIIntron_drtr_bpp_sum = boot(tn_GIIntron_drtr_bpp)

fn_GIIntron_drtr_bpp = (GIIntron[['drtr_fn']])
fn_GIIntron_drtr_bpp_sum = boot(fn_GIIntron_drtr_bpp)

#TPR
tpr_GIIntron_fold_mfe = tpr(tp_GIIntron_fold_mfe_sum, fn_GIIntron_fold_mfe_sum)
# print(tpr_GIIntron_fold_mfe)

tpr_GIIntron_fold_bpp = tpr(tp_GIIntron_fold_bpp_sum, fn_GIIntron_fold_bpp_sum)
# print(tpr_GIIntron_fold_bpp)

tpr_GIIntron_drtr_mfe = tpr(tp_GIIntron_drtr_mfe_sum, fn_GIIntron_drtr_mfe_sum)
# print(tpr_GIIntron_drtr_mfe)

tpr_GIIntron_drtr_prob = tpr(tp_GIIntron_drtr_prob_sum, fn_GIIntron_drtr_prob_sum)
# print(tpr_GIIntron_drtr_prob)

tpr_GIIntron_drtr_bpp = tpr(tp_GIIntron_drtr_bpp_sum, fn_GIIntron_drtr_bpp_sum)
# print(tpr_GIIntron_drtr_bpp)


#PPV
ppv_GIIntron_fold_mfe = ppv(tp_GIIntron_fold_mfe_sum, fn_GIIntron_fold_mfe_sum)
# print(ppv_GIIntron_fold_mfe)

ppv_GIIntron_fold_bpp = ppv(tp_GIIntron_fold_bpp_sum, fn_GIIntron_fold_bpp_sum)
# print(ppv_GIIntron_fold_bpp)

ppv_GIIntron_drtr_mfe = ppv(tp_GIIntron_drtr_mfe_sum, fn_GIIntron_drtr_mfe_sum)
# print(ppv_GIIntron_drtr_mfe)

ppv_GIIntron_drtr_prob = ppv(tp_GIIntron_drtr_prob_sum, fn_GIIntron_drtr_prob_sum)
# print(ppv_GIIntron_drtr_prob)

ppv_GIIntron_drtr_bpp = ppv(tp_GIIntron_drtr_bpp_sum, fn_GIIntron_drtr_bpp_sum)
# print(ppv_GIIntron_drtr_bpp)

#F1
f_1_GIIntron_fold_mfe = f_1(tp_GIIntron_fold_mfe_sum, fp_GIIntron_fold_mfe_sum, fn_GIIntron_fold_mfe_sum)
# print(f_1_GIIntron_fold_mfe)

f_1_GIIntron_fold_bpp = f_1(tp_GIIntron_fold_bpp_sum, fp_GIIntron_fold_bpp_sum, fn_GIIntron_fold_bpp_sum)
# print(f_1_GIIntron_fold_bpp)

f_1_GIIntron_drtr_mfe = f_1(tp_GIIntron_drtr_mfe_sum, fp_GIIntron_drtr_mfe_sum, fn_GIIntron_drtr_mfe_sum)
# print(f_1_GIIntron_drtr_mfe)

f_1_GIIntron_drtr_prob = f_1(tp_GIIntron_drtr_prob_sum, fp_GIIntron_drtr_prob_sum, fn_GIIntron_drtr_prob_sum)
# print(f_1_GIIntron_drtr_prob)

f_1_GIIntron_drtr_bpp = f_1(tp_GIIntron_drtr_bpp_sum, fp_GIIntron_drtr_bpp_sum, fn_GIIntron_drtr_bpp_sum)
# print(f_1_GIIntron_drtr_bpp)

# #MCC
mcc_GIIntron_fold_mfe = mcc(tp_GIIntron_fold_mfe_sum, fp_GIIntron_fold_mfe_sum, tn_GIIntron_fold_mfe_sum, fn_GIIntron_fold_mfe_sum)
# print(mcc_GIIntron_fold_mfe)

mcc_GIIntron_fold_bpp = mcc(tp_GIIntron_fold_bpp_sum, fp_GIIntron_fold_bpp_sum, tn_GIIntron_fold_bpp_sum, fn_GIIntron_fold_bpp_sum)
# print(mcc_GIIntron_fold_bpp)

mcc_GIIntron_drtr_mfe = mcc(tp_GIIntron_drtr_mfe_sum, fp_GIIntron_drtr_mfe_sum, tn_GIIntron_drtr_mfe_sum, fn_GIIntron_drtr_mfe_sum)
# print(mcc_GIIntron_drtr_mfe)

mcc_GIIntron_drtr_prob = mcc(tp_GIIntron_drtr_prob_sum, fp_GIIntron_drtr_prob_sum, tn_GIIntron_drtr_prob_sum, fn_GIIntron_drtr_prob_sum)
# print(mcc_GIIntron_drtr_prob)

mcc_GIIntron_drtr_bpp = mcc(tp_GIIntron_drtr_bpp_sum, fp_GIIntron_drtr_bpp_sum, tn_GIIntron_drtr_bpp_sum, fn_GIIntron_drtr_bpp_sum)
# print(mcc_GIIntron_drtr_bpp)


#Ham.Ribozyme
ham_ribozyme = family_data.loc[family_data.family=='Ham.Ribozyme']

tp_ham_ribozyme_fold_mfe = (ham_ribozyme[['RNAfold_mfe_tp']])
tp_ham_ribozyme_fold_mfe_sum = boot(tp_ham_ribozyme_fold_mfe)

fp_ham_ribozyme_fold_mfe = (ham_ribozyme[['RNAfold_mfe_fp']])
fp_ham_ribozyme_fold_mfe_sum = boot(fp_ham_ribozyme_fold_mfe)

tn_ham_ribozyme_fold_mfe = (ham_ribozyme[['RNAfold_mfe_tn']])
tn_ham_ribozyme_fold_mfe_sum = boot(tn_ham_ribozyme_fold_mfe)

fn_ham_ribozyme_fold_mfe = (ham_ribozyme[['RNAfold_mfe_fn']])
fn_ham_ribozyme_fold_mfe_sum = boot(fn_ham_ribozyme_fold_mfe)


tp_ham_ribozyme_fold_bpp = (ham_ribozyme[['RNAfold_bpp_tp']])
tp_ham_ribozyme_fold_bpp_sum = boot(tp_ham_ribozyme_fold_bpp)

fp_ham_ribozyme_fold_bpp = (ham_ribozyme[['RNAfold_bpp_fp']])
fp_ham_ribozyme_fold_bpp_sum = boot(fp_ham_ribozyme_fold_bpp)

tn_ham_ribozyme_fold_bpp = (ham_ribozyme[['RNAfold_bpp_tn']])
tn_ham_ribozyme_fold_bpp_sum = boot(tn_ham_ribozyme_fold_bpp)

fn_ham_ribozyme_fold_bpp = (ham_ribozyme[['RNAfold_bpp_fn']])
fn_ham_ribozyme_fold_bpp_sum = boot(fn_ham_ribozyme_fold_bpp)


tp_ham_ribozyme_drtr_mfe = (ham_ribozyme[['drtr_tp_mfe']])
tp_ham_ribozyme_drtr_mfe_sum = boot(tp_ham_ribozyme_drtr_mfe)

fp_ham_ribozyme_drtr_mfe = (ham_ribozyme[['drtr_fp_mfe']])
fp_ham_ribozyme_drtr_mfe_sum = boot(fp_ham_ribozyme_drtr_mfe)

tn_ham_ribozyme_drtr_mfe = (ham_ribozyme[['drtr_tn_mfe']])
tn_ham_ribozyme_drtr_mfe_sum = boot(tn_ham_ribozyme_drtr_mfe)

fn_ham_ribozyme_drtr_mfe = (ham_ribozyme[['drtr_fn_mfe']])
fn_ham_ribozyme_drtr_mfe_sum = boot(fn_ham_ribozyme_drtr_mfe)


tp_ham_ribozyme_drtr_prob = (ham_ribozyme[['drtr_tp_prob']])
tp_ham_ribozyme_drtr_prob_sum = boot(tp_ham_ribozyme_drtr_prob)

fp_ham_ribozyme_drtr_prob = (ham_ribozyme[['drtr_fp_prob']])
fp_ham_ribozyme_drtr_prob_sum = boot(fp_ham_ribozyme_drtr_prob)

tn_ham_ribozyme_drtr_prob = (ham_ribozyme[['drtr_tn_prob']])
tn_ham_ribozyme_drtr_prob_sum = boot(tn_ham_ribozyme_drtr_prob)

fn_ham_ribozyme_drtr_prob = (ham_ribozyme[['drtr_fn_prob']])
fn_ham_ribozyme_drtr_prob_sum = boot(fn_ham_ribozyme_drtr_prob)


tp_ham_ribozyme_drtr_bpp = (ham_ribozyme[['drtr_tp']])
tp_ham_ribozyme_drtr_bpp_sum = boot(tp_ham_ribozyme_drtr_bpp)

fp_ham_ribozyme_drtr_bpp = (ham_ribozyme[['drtr_fp']])
fp_ham_ribozyme_drtr_bpp_sum = boot(fp_ham_ribozyme_drtr_bpp)

tn_ham_ribozyme_drtr_bpp = (ham_ribozyme[['drtr_tn']])
tn_ham_ribozyme_drtr_bpp_sum = boot(tn_ham_ribozyme_drtr_bpp)

fn_ham_ribozyme_drtr_bpp = (ham_ribozyme[['drtr_fn']])
fn_ham_ribozyme_drtr_bpp_sum = boot(fn_ham_ribozyme_drtr_bpp)

#TPR
tpr_ham_ribozyme_fold_mfe = tpr(tp_ham_ribozyme_fold_mfe_sum, fn_ham_ribozyme_fold_mfe_sum)
# print(tpr_ham_ribozyme_fold_mfe)

tpr_ham_ribozyme_fold_bpp = tpr(tp_ham_ribozyme_fold_bpp_sum, fn_ham_ribozyme_fold_bpp_sum)
# print(tpr_ham_ribozyme_fold_bpp)

tpr_ham_ribozyme_drtr_mfe = tpr(tp_ham_ribozyme_drtr_mfe_sum, fn_ham_ribozyme_drtr_mfe_sum)
# print(tpr_ham_ribozyme_drtr_mfe)

tpr_ham_ribozyme_drtr_prob = tpr(tp_ham_ribozyme_drtr_prob_sum, fn_ham_ribozyme_drtr_prob_sum)
# print(tpr_ham_ribozyme_drtr_prob)

tpr_ham_ribozyme_drtr_bpp = tpr(tp_ham_ribozyme_drtr_bpp_sum, fn_ham_ribozyme_drtr_bpp_sum)
# print(tpr_ham_ribozyme_drtr_bpp)


#PPV
ppv_ham_ribozyme_fold_mfe = ppv(tp_ham_ribozyme_fold_mfe_sum, fn_ham_ribozyme_fold_mfe_sum)
# print(ppv_ham_ribozyme_fold_mfe)

ppv_ham_ribozyme_fold_bpp = ppv(tp_ham_ribozyme_fold_bpp_sum, fn_ham_ribozyme_fold_bpp_sum)
# print(ppv_ham_ribozyme_fold_bpp)

ppv_ham_ribozyme_drtr_mfe = ppv(tp_ham_ribozyme_drtr_mfe_sum, fn_ham_ribozyme_drtr_mfe_sum)
# print(ppv_ham_ribozyme_drtr_mfe)

ppv_ham_ribozyme_drtr_prob = ppv(tp_ham_ribozyme_drtr_prob_sum, fn_ham_ribozyme_drtr_prob_sum)
# print(ppv_ham_ribozyme_drtr_prob)

ppv_ham_ribozyme_drtr_bpp = ppv(tp_ham_ribozyme_drtr_bpp_sum, fn_ham_ribozyme_drtr_bpp_sum)
# print(ppv_ham_ribozyme_drtr_bpp)

#F1
f_1_ham_ribozyme_fold_mfe = f_1(tp_ham_ribozyme_fold_mfe_sum, fp_ham_ribozyme_fold_mfe_sum, fn_ham_ribozyme_fold_mfe_sum)
# print(f_1_ham_ribozyme_fold_mfe)

f_1_ham_ribozyme_fold_bpp = f_1(tp_ham_ribozyme_fold_bpp_sum, fp_ham_ribozyme_fold_bpp_sum, fn_ham_ribozyme_fold_bpp_sum)
# print(f_1_ham_ribozyme_fold_bpp)

f_1_ham_ribozyme_drtr_mfe = f_1(tp_ham_ribozyme_drtr_mfe_sum, fp_ham_ribozyme_drtr_mfe_sum, fn_ham_ribozyme_drtr_mfe_sum)
# print(f_1_ham_ribozyme_drtr_mfe)

f_1_ham_ribozyme_drtr_prob = f_1(tp_ham_ribozyme_drtr_prob_sum, fp_ham_ribozyme_drtr_prob_sum, fn_ham_ribozyme_drtr_prob_sum)
# print(f_1_ham_ribozyme_drtr_prob)

f_1_ham_ribozyme_drtr_bpp = f_1(tp_ham_ribozyme_drtr_bpp_sum, fp_ham_ribozyme_drtr_bpp_sum, fn_ham_ribozyme_drtr_bpp_sum)
# print(f_1_ham_ribozyme_drtr_bpp)

# #MCC
mcc_ham_ribozyme_fold_mfe = mcc(tp_ham_ribozyme_fold_mfe_sum, fp_ham_ribozyme_fold_mfe_sum, tn_ham_ribozyme_fold_mfe_sum, fn_ham_ribozyme_fold_mfe_sum)
# print(mcc_ham_ribozyme_fold_mfe)

mcc_ham_ribozyme_fold_bpp = mcc(tp_ham_ribozyme_fold_bpp_sum, fp_ham_ribozyme_fold_bpp_sum, tn_ham_ribozyme_fold_bpp_sum, fn_ham_ribozyme_fold_bpp_sum)
# print(mcc_ham_ribozyme_fold_bpp)

mcc_ham_ribozyme_drtr_mfe = mcc(tp_ham_ribozyme_drtr_mfe_sum, fp_ham_ribozyme_drtr_mfe_sum, tn_ham_ribozyme_drtr_mfe_sum, fn_ham_ribozyme_drtr_mfe_sum)
# print(mcc_ham_ribozyme_drtr_mfe)

mcc_ham_ribozyme_drtr_prob = mcc(tp_ham_ribozyme_drtr_prob_sum, fp_ham_ribozyme_drtr_prob_sum, tn_ham_ribozyme_drtr_prob_sum, fn_ham_ribozyme_drtr_prob_sum)
# print(mcc_ham_ribozyme_drtr_prob)

mcc_ham_ribozyme_drtr_bpp = mcc(tp_ham_ribozyme_drtr_bpp_sum, fp_ham_ribozyme_drtr_bpp_sum, tn_ham_ribozyme_drtr_bpp_sum, fn_ham_ribozyme_drtr_bpp_sum)
# print(mcc_ham_ribozyme_drtr_bpp)


#HDVRibozyme
hdvr_ribozyme=family_data.loc[family_data.family=='HDVRibozyme']

tp_hdvr_ribozyme_fold_mfe = (hdvr_ribozyme[['RNAfold_mfe_tp']])
tp_hdvr_ribozyme_fold_mfe_sum = boot(tp_hdvr_ribozyme_fold_mfe)

fp_hdvr_ribozyme_fold_mfe = (hdvr_ribozyme[['RNAfold_mfe_fp']])
fp_hdvr_ribozyme_fold_mfe_sum = boot(fp_hdvr_ribozyme_fold_mfe)

tn_hdvr_ribozyme_fold_mfe = (hdvr_ribozyme[['RNAfold_mfe_tn']])
tn_hdvr_ribozyme_fold_mfe_sum = boot(tn_hdvr_ribozyme_fold_mfe)

fn_hdvr_ribozyme_fold_mfe = (hdvr_ribozyme[['RNAfold_mfe_fn']])
fn_hdvr_ribozyme_fold_mfe_sum = boot(fn_hdvr_ribozyme_fold_mfe)


tp_hdvr_ribozyme_fold_bpp = (hdvr_ribozyme[['RNAfold_bpp_tp']])
tp_hdvr_ribozyme_fold_bpp_sum = boot(tp_hdvr_ribozyme_fold_bpp)

fp_hdvr_ribozyme_fold_bpp = (hdvr_ribozyme[['RNAfold_bpp_fp']])
fp_hdvr_ribozyme_fold_bpp_sum = boot(fp_hdvr_ribozyme_fold_bpp)

tn_hdvr_ribozyme_fold_bpp = (hdvr_ribozyme[['RNAfold_bpp_tn']])
tn_hdvr_ribozyme_fold_bpp_sum = boot(tn_hdvr_ribozyme_fold_bpp)

fn_hdvr_ribozyme_fold_bpp = (hdvr_ribozyme[['RNAfold_bpp_fn']])
fn_hdvr_ribozyme_fold_bpp_sum = boot(fn_hdvr_ribozyme_fold_bpp)


tp_hdvr_ribozyme_drtr_mfe = (hdvr_ribozyme[['drtr_tp_mfe']])
tp_hdvr_ribozyme_drtr_mfe_sum = boot(tp_hdvr_ribozyme_drtr_mfe)

fp_hdvr_ribozyme_drtr_mfe = (hdvr_ribozyme[['drtr_fp_mfe']])
fp_hdvr_ribozyme_drtr_mfe_sum = boot(fp_hdvr_ribozyme_drtr_mfe)

tn_hdvr_ribozyme_drtr_mfe = (hdvr_ribozyme[['drtr_tn_mfe']])
tn_hdvr_ribozyme_drtr_mfe_sum = boot(tn_hdvr_ribozyme_drtr_mfe)

fn_hdvr_ribozyme_drtr_mfe = (hdvr_ribozyme[['drtr_fn_mfe']])
fn_hdvr_ribozyme_drtr_mfe_sum = boot(fn_hdvr_ribozyme_drtr_mfe)


tp_hdvr_ribozyme_drtr_prob = (hdvr_ribozyme[['drtr_tp_prob']])
tp_hdvr_ribozyme_drtr_prob_sum = boot(tp_hdvr_ribozyme_drtr_prob)

fp_hdvr_ribozyme_drtr_prob = (hdvr_ribozyme[['drtr_fp_prob']])
fp_hdvr_ribozyme_drtr_prob_sum = boot(fp_hdvr_ribozyme_drtr_prob)

tn_hdvr_ribozyme_drtr_prob = (hdvr_ribozyme[['drtr_tn_prob']])
tn_hdvr_ribozyme_drtr_prob_sum = boot(tn_hdvr_ribozyme_drtr_prob)

fn_hdvr_ribozyme_drtr_prob = (hdvr_ribozyme[['drtr_fn_prob']])
fn_hdvr_ribozyme_drtr_prob_sum = boot(fn_hdvr_ribozyme_drtr_prob)


tp_hdvr_ribozyme_drtr_bpp = (hdvr_ribozyme[['drtr_tp']])
tp_hdvr_ribozyme_drtr_bpp_sum = boot(tp_hdvr_ribozyme_drtr_bpp)

fp_hdvr_ribozyme_drtr_bpp = (hdvr_ribozyme[['drtr_fp']])
fp_hdvr_ribozyme_drtr_bpp_sum = boot(fp_hdvr_ribozyme_drtr_bpp)

tn_hdvr_ribozyme_drtr_bpp = (hdvr_ribozyme[['drtr_tn']])
tn_hdvr_ribozyme_drtr_bpp_sum = boot(tn_hdvr_ribozyme_drtr_bpp)

fn_hdvr_ribozyme_drtr_bpp = (hdvr_ribozyme[['drtr_fn']])
fn_hdvr_ribozyme_drtr_bpp_sum = boot(fn_hdvr_ribozyme_drtr_bpp)

#TPR
tpr_hdvr_ribozyme_fold_mfe = tpr(tp_hdvr_ribozyme_fold_mfe_sum, fn_hdvr_ribozyme_fold_mfe_sum)
# print(tpr_hdvr_ribozyme_fold_mfe)

tpr_hdvr_ribozyme_fold_bpp = tpr(tp_hdvr_ribozyme_fold_bpp_sum, fn_hdvr_ribozyme_fold_bpp_sum)
# print(tpr_hdvr_ribozyme_fold_bpp)

tpr_hdvr_ribozyme_drtr_mfe = tpr(tp_hdvr_ribozyme_drtr_mfe_sum, fn_hdvr_ribozyme_drtr_mfe_sum)
# print(tpr_hdvr_ribozyme_drtr_mfe)

tpr_hdvr_ribozyme_drtr_prob = tpr(tp_hdvr_ribozyme_drtr_prob_sum, fn_hdvr_ribozyme_drtr_prob_sum)
# print(tpr_hdvr_ribozyme_drtr_prob)

tpr_hdvr_ribozyme_drtr_bpp = tpr(tp_hdvr_ribozyme_drtr_bpp_sum, fn_hdvr_ribozyme_drtr_bpp_sum)
# print(tpr_hdvr_ribozyme_drtr_bpp)


#PPV
ppv_hdvr_ribozyme_fold_mfe = ppv(tp_hdvr_ribozyme_fold_mfe_sum, fn_hdvr_ribozyme_fold_mfe_sum)
# print(ppv_hdvr_ribozyme_fold_mfe)

ppv_hdvr_ribozyme_fold_bpp = ppv(tp_hdvr_ribozyme_fold_bpp_sum, fn_hdvr_ribozyme_fold_bpp_sum)
# print(ppv_hdvr_ribozyme_fold_bpp)

ppv_hdvr_ribozyme_drtr_mfe = ppv(tp_hdvr_ribozyme_drtr_mfe_sum, fn_hdvr_ribozyme_drtr_mfe_sum)
# print(ppv_hdvr_ribozyme_drtr_mfe)

ppv_hdvr_ribozyme_drtr_prob = ppv(tp_hdvr_ribozyme_drtr_prob_sum, fn_hdvr_ribozyme_drtr_prob_sum)
# print(ppv_hdvr_ribozyme_drtr_prob)

ppv_hdvr_ribozyme_drtr_bpp = ppv(tp_hdvr_ribozyme_drtr_bpp_sum, fn_hdvr_ribozyme_drtr_bpp_sum)
# print(ppv_hdvr_ribozyme_drtr_bpp)

#F1
f_1_hdvr_ribozyme_fold_mfe = f_1(tp_hdvr_ribozyme_fold_mfe_sum, fp_hdvr_ribozyme_fold_mfe_sum, fn_hdvr_ribozyme_fold_mfe_sum)
# print(f_1_hdvr_ribozyme_fold_mfe)

f_1_hdvr_ribozyme_fold_bpp = f_1(tp_hdvr_ribozyme_fold_bpp_sum, fp_hdvr_ribozyme_fold_bpp_sum, fn_hdvr_ribozyme_fold_bpp_sum)
# print(f_1_hdvr_ribozyme_fold_bpp)

f_1_hdvr_ribozyme_drtr_mfe = f_1(tp_hdvr_ribozyme_drtr_mfe_sum, fp_hdvr_ribozyme_drtr_mfe_sum, fn_hdvr_ribozyme_drtr_mfe_sum)
# print(f_1_hdvr_ribozyme_drtr_mfe)

f_1_hdvr_ribozyme_drtr_prob = f_1(tp_hdvr_ribozyme_drtr_prob_sum, fp_hdvr_ribozyme_drtr_prob_sum, fn_hdvr_ribozyme_drtr_prob_sum)
# print(f_1_hdvr_ribozyme_drtr_prob)

f_1_hdvr_ribozyme_drtr_bpp = f_1(tp_hdvr_ribozyme_drtr_bpp_sum, fp_hdvr_ribozyme_drtr_bpp_sum, fn_hdvr_ribozyme_drtr_bpp_sum)
# print(f_1_hdvr_ribozyme_drtr_bpp)

# #MCC
mcc_hdvr_ribozyme_fold_mfe = mcc(tp_hdvr_ribozyme_fold_mfe_sum, fp_hdvr_ribozyme_fold_mfe_sum, tn_hdvr_ribozyme_fold_mfe_sum, fn_hdvr_ribozyme_fold_mfe_sum)
# print(mcc_hdvr_ribozyme_fold_mfe)

mcc_hdvr_ribozyme_fold_bpp = mcc(tp_hdvr_ribozyme_fold_bpp_sum, fp_hdvr_ribozyme_fold_bpp_sum, tn_hdvr_ribozyme_fold_bpp_sum, fn_hdvr_ribozyme_fold_bpp_sum)
# print(mcc_hdvr_ribozyme_fold_bpp)

mcc_hdvr_ribozyme_drtr_mfe = mcc(tp_hdvr_ribozyme_drtr_mfe_sum, fp_hdvr_ribozyme_drtr_mfe_sum, tn_hdvr_ribozyme_drtr_mfe_sum, fn_hdvr_ribozyme_drtr_mfe_sum)
# print(mcc_hdvr_ribozyme_drtr_mfe)

mcc_hdvr_ribozyme_drtr_prob = mcc(tp_hdvr_ribozyme_drtr_prob_sum, fp_hdvr_ribozyme_drtr_prob_sum, tn_hdvr_ribozyme_drtr_prob_sum, fn_hdvr_ribozyme_drtr_prob_sum)
# print(mcc_hdvr_ribozyme_drtr_prob)

mcc_hdvr_ribozyme_drtr_bpp = mcc(tp_hdvr_ribozyme_drtr_bpp_sum, fp_hdvr_ribozyme_drtr_bpp_sum, tn_hdvr_ribozyme_drtr_bpp_sum, fn_hdvr_ribozyme_drtr_bpp_sum)
# print(mcc_hdvr_ribozyme_drtr_bpp)

#IRES
ires=family_data.loc[family_data.family=='IRES']

tp_ires_fold_mfe = (ires[['RNAfold_mfe_tp']])
tp_ires_fold_mfe_sum = boot(tp_ires_fold_mfe)

fp_ires_fold_mfe = (ires[['RNAfold_mfe_fp']])
fp_ires_fold_mfe_sum = boot(fp_ires_fold_mfe)

tn_ires_fold_mfe = (ires[['RNAfold_mfe_tn']])
tn_ires_fold_mfe_sum = boot(tn_ires_fold_mfe)

fn_ires_fold_mfe = (ires[['RNAfold_mfe_fn']])
fn_ires_fold_mfe_sum = boot(fn_ires_fold_mfe)


tp_ires_fold_bpp = (ires[['RNAfold_bpp_tp']])
tp_ires_fold_bpp_sum = boot(tp_ires_fold_bpp)

fp_ires_fold_bpp = (ires[['RNAfold_bpp_fp']])
fp_ires_fold_bpp_sum = boot(fp_ires_fold_bpp)

tn_ires_fold_bpp = (ires[['RNAfold_bpp_tn']])
tn_ires_fold_bpp_sum = boot(tn_ires_fold_bpp)

fn_ires_fold_bpp = (ires[['RNAfold_bpp_fn']])
fn_ires_fold_bpp_sum = boot(fn_ires_fold_bpp)


tp_ires_drtr_mfe = (ires[['drtr_tp_mfe']])
tp_ires_drtr_mfe_sum = boot(tp_ires_drtr_mfe)

fp_ires_drtr_mfe = (ires[['drtr_fp_mfe']])
fp_ires_drtr_mfe_sum = boot(fp_ires_drtr_mfe)

tn_ires_drtr_mfe = (ires[['drtr_tn_mfe']])
tn_ires_drtr_mfe_sum = boot(tn_ires_drtr_mfe)

fn_ires_drtr_mfe = (ires[['drtr_fn_mfe']])
fn_ires_drtr_mfe_sum = boot(fn_ires_drtr_mfe)


tp_ires_drtr_prob = (ires[['drtr_tp_prob']])
tp_ires_drtr_prob_sum = boot(tp_ires_drtr_prob)

fp_ires_drtr_prob = (ires[['drtr_fp_prob']])
fp_ires_drtr_prob_sum = boot(fp_ires_drtr_prob)

tn_ires_drtr_prob = (ires[['drtr_tn_prob']])
tn_ires_drtr_prob_sum = boot(tn_ires_drtr_prob)

fn_ires_drtr_prob = (ires[['drtr_fn_prob']])
fn_ires_drtr_prob_sum = boot(fn_ires_drtr_prob)


tp_ires_drtr_bpp = (ires[['drtr_tp']])
tp_ires_drtr_bpp_sum = boot(tp_ires_drtr_bpp)

fp_ires_drtr_bpp = (ires[['drtr_fp']])
fp_ires_drtr_bpp_sum = boot(fp_ires_drtr_bpp)

tn_ires_drtr_bpp = (ires[['drtr_tn']])
tn_ires_drtr_bpp_sum = boot(tn_ires_drtr_bpp)

fn_ires_drtr_bpp = (ires[['drtr_fn']])
fn_ires_drtr_bpp_sum = boot(fn_ires_drtr_bpp)

#TPR
tpr_ires_fold_mfe = tpr(tp_ires_fold_mfe_sum, fn_ires_fold_mfe_sum)
# print(tpr_ires_fold_mfe)

tpr_ires_fold_bpp = tpr(tp_ires_fold_bpp_sum, fn_ires_fold_bpp_sum)
# print(tpr_ires_fold_bpp)

tpr_ires_drtr_mfe = tpr(tp_ires_drtr_mfe_sum, fn_ires_drtr_mfe_sum)
# print(tpr_ires_drtr_mfe)

tpr_ires_drtr_prob = tpr(tp_ires_drtr_prob_sum, fn_ires_drtr_prob_sum)
# print(tpr_ires_drtr_prob)

tpr_ires_drtr_bpp = tpr(tp_ires_drtr_bpp_sum, fn_ires_drtr_bpp_sum)
# print(tpr_ires_drtr_bpp)


#PPV
ppv_ires_fold_mfe = ppv(tp_ires_fold_mfe_sum, fn_ires_fold_mfe_sum)
# print(ppv_ires_fold_mfe)

ppv_ires_fold_bpp = ppv(tp_ires_fold_bpp_sum, fn_ires_fold_bpp_sum)
# print(ppv_ires_fold_bpp)

ppv_ires_drtr_mfe = ppv(tp_ires_drtr_mfe_sum, fn_ires_drtr_mfe_sum)
# print(ppv_ires_drtr_mfe)

ppv_ires_drtr_prob = ppv(tp_ires_drtr_prob_sum, fn_ires_drtr_prob_sum)
# print(ppv_ires_drtr_prob)

ppv_ires_drtr_bpp = ppv(tp_ires_drtr_bpp_sum, fn_ires_drtr_bpp_sum)
# print(ppv_ires_drtr_bpp)

#F1
f_1_ires_fold_mfe = f_1(tp_ires_fold_mfe_sum, fp_ires_fold_mfe_sum, fn_ires_fold_mfe_sum)
# print(f_1_ires_fold_mfe)

f_1_ires_fold_bpp = f_1(tp_ires_fold_bpp_sum, fp_ires_fold_bpp_sum, fn_ires_fold_bpp_sum)
# print(f_1_ires_fold_bpp)

f_1_ires_drtr_mfe = f_1(tp_ires_drtr_mfe_sum, fp_ires_drtr_mfe_sum, fn_ires_drtr_mfe_sum)
# print(f_1_ires_drtr_mfe)

f_1_ires_drtr_prob = f_1(tp_ires_drtr_prob_sum, fp_ires_drtr_prob_sum, fn_ires_drtr_prob_sum)
# print(f_1_ires_drtr_prob)

f_1_ires_drtr_bpp = f_1(tp_ires_drtr_bpp_sum, fp_ires_drtr_bpp_sum, fn_ires_drtr_bpp_sum)
# print(f_1_ires_drtr_bpp)

# #MCC
mcc_ires_fold_mfe = mcc(tp_ires_fold_mfe_sum, fp_ires_fold_mfe_sum, tn_ires_fold_mfe_sum, fn_ires_fold_mfe_sum)
# print(mcc_ires_fold_mfe)

mcc_ires_fold_bpp = mcc(tp_ires_fold_bpp_sum, fp_ires_fold_bpp_sum, tn_ires_fold_bpp_sum, fn_ires_fold_bpp_sum)
# print(mcc_ires_fold_bpp)

mcc_ires_drtr_mfe = mcc(tp_ires_drtr_mfe_sum, fp_ires_drtr_mfe_sum, tn_ires_drtr_mfe_sum, fn_ires_drtr_mfe_sum)
# print(mcc_ires_drtr_mfe)

mcc_ires_drtr_prob = mcc(tp_ires_drtr_prob_sum, fp_ires_drtr_prob_sum, tn_ires_drtr_prob_sum, fn_ires_drtr_prob_sum)
# print(mcc_ires_drtr_prob)

mcc_ires_drtr_bpp = mcc(tp_ires_drtr_bpp_sum, fp_ires_drtr_bpp_sum, tn_ires_drtr_bpp_sum, fn_ires_drtr_bpp_sum)
# print(mcc_ires_drtr_bpp)


#OtherRibozyme
other_ribozyme=family_data.loc[family_data.family=='OtherRibozyme']

tp_other_ribozyme_fold_mfe = (other_ribozyme[['RNAfold_mfe_tp']])
tp_other_ribozyme_fold_mfe_sum = boot(tp_other_ribozyme_fold_mfe)

fp_other_ribozyme_fold_mfe = (other_ribozyme[['RNAfold_mfe_fp']])
fp_other_ribozyme_fold_mfe_sum = boot(fp_other_ribozyme_fold_mfe)

tn_other_ribozyme_fold_mfe = (other_ribozyme[['RNAfold_mfe_tn']])
tn_other_ribozyme_fold_mfe_sum = boot(tn_other_ribozyme_fold_mfe)

fn_other_ribozyme_fold_mfe = (other_ribozyme[['RNAfold_mfe_fn']])
fn_other_ribozyme_fold_mfe_sum = boot(fn_other_ribozyme_fold_mfe)


tp_other_ribozyme_fold_bpp = (other_ribozyme[['RNAfold_bpp_tp']])
tp_other_ribozyme_fold_bpp_sum = boot(tp_other_ribozyme_fold_bpp)

fp_other_ribozyme_fold_bpp = (other_ribozyme[['RNAfold_bpp_fp']])
fp_other_ribozyme_fold_bpp_sum = boot(fp_other_ribozyme_fold_bpp)

tn_other_ribozyme_fold_bpp = (other_ribozyme[['RNAfold_bpp_tn']])
tn_other_ribozyme_fold_bpp_sum = boot(tn_other_ribozyme_fold_bpp)

fn_other_ribozyme_fold_bpp = (other_ribozyme[['RNAfold_bpp_fn']])
fn_other_ribozyme_fold_bpp_sum = boot(fn_other_ribozyme_fold_bpp)


tp_other_ribozyme_drtr_mfe = (other_ribozyme[['drtr_tp_mfe']])
tp_other_ribozyme_drtr_mfe_sum = boot(tp_other_ribozyme_drtr_mfe)

fp_other_ribozyme_drtr_mfe = (other_ribozyme[['drtr_fp_mfe']])
fp_other_ribozyme_drtr_mfe_sum = boot(fp_other_ribozyme_drtr_mfe)

tn_other_ribozyme_drtr_mfe = (other_ribozyme[['drtr_tn_mfe']])
tn_other_ribozyme_drtr_mfe_sum = boot(tn_other_ribozyme_drtr_mfe)

fn_other_ribozyme_drtr_mfe = (other_ribozyme[['drtr_fn_mfe']])
fn_other_ribozyme_drtr_mfe_sum = boot(fn_other_ribozyme_drtr_mfe)


tp_other_ribozyme_drtr_prob = (other_ribozyme[['drtr_tp_prob']])
tp_other_ribozyme_drtr_prob_sum = boot(tp_other_ribozyme_drtr_prob)

fp_other_ribozyme_drtr_prob = (other_ribozyme[['drtr_fp_prob']])
fp_other_ribozyme_drtr_prob_sum = boot(fp_other_ribozyme_drtr_prob)

tn_other_ribozyme_drtr_prob = (other_ribozyme[['drtr_tn_prob']])
tn_other_ribozyme_drtr_prob_sum = boot(tn_other_ribozyme_drtr_prob)

fn_other_ribozyme_drtr_prob = (other_ribozyme[['drtr_fn_prob']])
fn_other_ribozyme_drtr_prob_sum = boot(fn_other_ribozyme_drtr_prob)


tp_other_ribozyme_drtr_bpp = (other_ribozyme[['drtr_tp']])
tp_other_ribozyme_drtr_bpp_sum = boot(tp_other_ribozyme_drtr_bpp)

fp_other_ribozyme_drtr_bpp = (other_ribozyme[['drtr_fp']])
fp_other_ribozyme_drtr_bpp_sum = boot(fp_other_ribozyme_drtr_bpp)

tn_other_ribozyme_drtr_bpp = (other_ribozyme[['drtr_tn']])
tn_other_ribozyme_drtr_bpp_sum = boot(tn_other_ribozyme_drtr_bpp)

fn_other_ribozyme_drtr_bpp = (other_ribozyme[['drtr_fn']])
fn_other_ribozyme_drtr_bpp_sum = boot(fn_other_ribozyme_drtr_bpp)

#TPR
tpr_other_ribozyme_fold_mfe = tpr(tp_other_ribozyme_fold_mfe_sum, fn_other_ribozyme_fold_mfe_sum)
# print(tpr_other_ribozyme_fold_mfe)

tpr_other_ribozyme_fold_bpp = tpr(tp_other_ribozyme_fold_bpp_sum, fn_other_ribozyme_fold_bpp_sum)
# print(tpr_other_ribozyme_fold_bpp)

tpr_other_ribozyme_drtr_mfe = tpr(tp_other_ribozyme_drtr_mfe_sum, fn_other_ribozyme_drtr_mfe_sum)
# print(tpr_other_ribozyme_drtr_mfe)

tpr_other_ribozyme_drtr_prob = tpr(tp_other_ribozyme_drtr_prob_sum, fn_other_ribozyme_drtr_prob_sum)
# print(tpr_other_ribozyme_drtr_prob)

tpr_other_ribozyme_drtr_bpp = tpr(tp_other_ribozyme_drtr_bpp_sum, fn_other_ribozyme_drtr_bpp_sum)
# print(tpr_other_ribozyme_drtr_bpp)


#PPV
ppv_other_ribozyme_fold_mfe = ppv(tp_other_ribozyme_fold_mfe_sum, fn_other_ribozyme_fold_mfe_sum)
# print(ppv_other_ribozyme_fold_mfe)

ppv_other_ribozyme_fold_bpp = ppv(tp_other_ribozyme_fold_bpp_sum, fn_other_ribozyme_fold_bpp_sum)
# print(ppv_other_ribozyme_fold_bpp)

ppv_other_ribozyme_drtr_mfe = ppv(tp_other_ribozyme_drtr_mfe_sum, fn_other_ribozyme_drtr_mfe_sum)
# print(ppv_other_ribozyme_drtr_mfe)

ppv_other_ribozyme_drtr_prob = ppv(tp_other_ribozyme_drtr_prob_sum, fn_other_ribozyme_drtr_prob_sum)
# print(ppv_other_ribozyme_drtr_prob)

ppv_other_ribozyme_drtr_bpp = ppv(tp_other_ribozyme_drtr_bpp_sum, fn_other_ribozyme_drtr_bpp_sum)
# print(ppv_other_ribozyme_drtr_bpp)

#F1
f_1_other_ribozyme_fold_mfe = f_1(tp_other_ribozyme_fold_mfe_sum, fp_other_ribozyme_fold_mfe_sum, fn_other_ribozyme_fold_mfe_sum)
# print(f_1_other_ribozyme_fold_mfe)

f_1_other_ribozyme_fold_bpp = f_1(tp_other_ribozyme_fold_bpp_sum, fp_other_ribozyme_fold_bpp_sum, fn_other_ribozyme_fold_bpp_sum)
# print(f_1_other_ribozyme_fold_bpp)

f_1_other_ribozyme_drtr_mfe = f_1(tp_other_ribozyme_drtr_mfe_sum, fp_other_ribozyme_drtr_mfe_sum, fn_other_ribozyme_drtr_mfe_sum)
# print(f_1_other_ribozyme_drtr_mfe)

f_1_other_ribozyme_drtr_prob = f_1(tp_other_ribozyme_drtr_prob_sum, fp_other_ribozyme_drtr_prob_sum, fn_other_ribozyme_drtr_prob_sum)
# print(f_1_other_ribozyme_drtr_prob)

f_1_other_ribozyme_drtr_bpp = f_1(tp_other_ribozyme_drtr_bpp_sum, fp_other_ribozyme_drtr_bpp_sum, fn_other_ribozyme_drtr_bpp_sum)
# print(f_1_other_ribozyme_drtr_bpp)

# #MCC
mcc_other_ribozyme_fold_mfe = mcc(tp_other_ribozyme_fold_mfe_sum, fp_other_ribozyme_fold_mfe_sum, tn_other_ribozyme_fold_mfe_sum, fn_other_ribozyme_fold_mfe_sum)
# print(mcc_other_ribozyme_fold_mfe)

mcc_other_ribozyme_fold_bpp = mcc(tp_other_ribozyme_fold_bpp_sum, fp_other_ribozyme_fold_bpp_sum, tn_other_ribozyme_fold_bpp_sum, fn_other_ribozyme_fold_bpp_sum)
# print(mcc_other_ribozyme_fold_bpp)

mcc_other_ribozyme_drtr_mfe = mcc(tp_other_ribozyme_drtr_mfe_sum, fp_other_ribozyme_drtr_mfe_sum, tn_other_ribozyme_drtr_mfe_sum, fn_other_ribozyme_drtr_mfe_sum)
# print(mcc_other_ribozyme_drtr_mfe)

mcc_other_ribozyme_drtr_prob = mcc(tp_other_ribozyme_drtr_prob_sum, fp_other_ribozyme_drtr_prob_sum, tn_other_ribozyme_drtr_prob_sum, fn_other_ribozyme_drtr_prob_sum)
# print(mcc_other_ribozyme_drtr_prob)

mcc_other_ribozyme_drtr_bpp = mcc(tp_other_ribozyme_drtr_bpp_sum, fp_other_ribozyme_drtr_bpp_sum, tn_other_ribozyme_drtr_bpp_sum, fn_other_ribozyme_drtr_bpp_sum)
# print(mcc_other_ribozyme_drtr_bpp)



#OtherRNA
other_RNA=family_data.loc[family_data.family=='OtherRNA']

tp_other_RNA_fold_mfe = (other_RNA[['RNAfold_mfe_tp']])
tp_other_RNA_fold_mfe_sum = boot(tp_other_RNA_fold_mfe)

fp_other_RNA_fold_mfe = (other_RNA[['RNAfold_mfe_fp']])
fp_other_RNA_fold_mfe_sum = boot(fp_other_RNA_fold_mfe)

tn_other_RNA_fold_mfe = (other_RNA[['RNAfold_mfe_tn']])
tn_other_RNA_fold_mfe_sum = boot(tn_other_RNA_fold_mfe)

fn_other_RNA_fold_mfe = (other_RNA[['RNAfold_mfe_fn']])
fn_other_RNA_fold_mfe_sum = boot(fn_other_RNA_fold_mfe)


tp_other_RNA_fold_bpp = (other_RNA[['RNAfold_bpp_tp']])
tp_other_RNA_fold_bpp_sum = boot(tp_other_RNA_fold_bpp)

fp_other_RNA_fold_bpp = (other_RNA[['RNAfold_bpp_fp']])
fp_other_RNA_fold_bpp_sum = boot(fp_other_RNA_fold_bpp)

tn_other_RNA_fold_bpp = (other_RNA[['RNAfold_bpp_tn']])
tn_other_RNA_fold_bpp_sum = boot(tn_other_RNA_fold_bpp)

fn_other_RNA_fold_bpp = (other_RNA[['RNAfold_bpp_fn']])
fn_other_RNA_fold_bpp_sum = boot(fn_other_RNA_fold_bpp)


tp_other_RNA_drtr_mfe = (other_RNA[['drtr_tp_mfe']])
tp_other_RNA_drtr_mfe_sum = boot(tp_other_RNA_drtr_mfe)

fp_other_RNA_drtr_mfe = (other_RNA[['drtr_fp_mfe']])
fp_other_RNA_drtr_mfe_sum = boot(fp_other_RNA_drtr_mfe)

tn_other_RNA_drtr_mfe = (other_RNA[['drtr_tn_mfe']])
tn_other_RNA_drtr_mfe_sum = boot(tn_other_RNA_drtr_mfe)

fn_other_RNA_drtr_mfe = (other_RNA[['drtr_fn_mfe']])
fn_other_RNA_drtr_mfe_sum = boot(fn_other_RNA_drtr_mfe)


tp_other_RNA_drtr_prob = (other_RNA[['drtr_tp_prob']])
tp_other_RNA_drtr_prob_sum = boot(tp_other_RNA_drtr_prob)

fp_other_RNA_drtr_prob = (other_RNA[['drtr_fp_prob']])
fp_other_RNA_drtr_prob_sum = boot(fp_other_RNA_drtr_prob)

tn_other_RNA_drtr_prob = (other_RNA[['drtr_tn_prob']])
tn_other_RNA_drtr_prob_sum = boot(tn_other_RNA_drtr_prob)

fn_other_RNA_drtr_prob = (other_RNA[['drtr_fn_prob']])
fn_other_RNA_drtr_prob_sum = boot(fn_other_RNA_drtr_prob)


tp_other_RNA_drtr_bpp = (other_RNA[['drtr_tp']])
tp_other_RNA_drtr_bpp_sum = boot(tp_other_RNA_drtr_bpp)

fp_other_RNA_drtr_bpp = (other_RNA[['drtr_fp']])
fp_other_RNA_drtr_bpp_sum = boot(fp_other_RNA_drtr_bpp)

tn_other_RNA_drtr_bpp = (other_RNA[['drtr_tn']])
tn_other_RNA_drtr_bpp_sum = boot(tn_other_RNA_drtr_bpp)

fn_other_RNA_drtr_bpp = (other_RNA[['drtr_fn']])
fn_other_RNA_drtr_bpp_sum = boot(fn_other_RNA_drtr_bpp)

#TPR
tpr_other_RNA_fold_mfe = tpr(tp_other_RNA_fold_mfe_sum, fn_other_RNA_fold_mfe_sum)
# print(tpr_other_RNA_fold_mfe)

tpr_other_RNA_fold_bpp = tpr(tp_other_RNA_fold_bpp_sum, fn_other_RNA_fold_bpp_sum)
# print(tpr_other_RNA_fold_bpp)

tpr_other_RNA_drtr_mfe = tpr(tp_other_RNA_drtr_mfe_sum, fn_other_RNA_drtr_mfe_sum)
# print(tpr_other_RNA_drtr_mfe)

tpr_other_RNA_drtr_prob = tpr(tp_other_RNA_drtr_prob_sum, fn_other_RNA_drtr_prob_sum)
# print(tpr_other_RNA_drtr_prob)

tpr_other_RNA_drtr_bpp = tpr(tp_other_RNA_drtr_bpp_sum, fn_other_RNA_drtr_bpp_sum)
# print(tpr_other_RNA_drtr_bpp)


#PPV
ppv_other_RNA_fold_mfe = ppv(tp_other_RNA_fold_mfe_sum, fn_other_RNA_fold_mfe_sum)
# print(ppv_other_RNA_fold_mfe)

ppv_other_RNA_fold_bpp = ppv(tp_other_RNA_fold_bpp_sum, fn_other_RNA_fold_bpp_sum)
# print(ppv_other_RNA_fold_bpp)

ppv_other_RNA_drtr_mfe = ppv(tp_other_RNA_drtr_mfe_sum, fn_other_RNA_drtr_mfe_sum)
# print(ppv_other_RNA_drtr_mfe)

ppv_other_RNA_drtr_prob = ppv(tp_other_RNA_drtr_prob_sum, fn_other_RNA_drtr_prob_sum)
# print(ppv_other_RNA_drtr_prob)

ppv_other_RNA_drtr_bpp = ppv(tp_other_RNA_drtr_bpp_sum, fn_other_RNA_drtr_bpp_sum)
# print(ppv_other_RNA_drtr_bpp)

#F1
f_1_other_RNA_fold_mfe = f_1(tp_other_RNA_fold_mfe_sum, fp_other_RNA_fold_mfe_sum, fn_other_RNA_fold_mfe_sum)
# print(f_1_other_RNA_fold_mfe)

f_1_other_RNA_fold_bpp = f_1(tp_other_RNA_fold_bpp_sum, fp_other_RNA_fold_bpp_sum, fn_other_RNA_fold_bpp_sum)
# print(f_1_other_RNA_fold_bpp)

f_1_other_RNA_drtr_mfe = f_1(tp_other_RNA_drtr_mfe_sum, fp_other_RNA_drtr_mfe_sum, fn_other_RNA_drtr_mfe_sum)
# print(f_1_other_RNA_drtr_mfe)

f_1_other_RNA_drtr_prob = f_1(tp_other_RNA_drtr_prob_sum, fp_other_RNA_drtr_prob_sum, fn_other_RNA_drtr_prob_sum)
# print(f_1_other_RNA_drtr_prob)

f_1_other_RNA_drtr_bpp = f_1(tp_other_RNA_drtr_bpp_sum, fp_other_RNA_drtr_bpp_sum, fn_other_RNA_drtr_bpp_sum)
# print(f_1_other_RNA_drtr_bpp)

# #MCC
mcc_other_RNA_fold_mfe = mcc(tp_other_RNA_fold_mfe_sum, fp_other_RNA_fold_mfe_sum, tn_other_RNA_fold_mfe_sum, fn_other_RNA_fold_mfe_sum)
# print(mcc_other_RNA_fold_mfe)

mcc_other_RNA_fold_bpp = mcc(tp_other_RNA_fold_bpp_sum, fp_other_RNA_fold_bpp_sum, tn_other_RNA_fold_bpp_sum, fn_other_RNA_fold_bpp_sum)
# print(mcc_other_RNA_fold_bpp)

mcc_other_RNA_drtr_mfe = mcc(tp_other_RNA_drtr_mfe_sum, fp_other_RNA_drtr_mfe_sum, tn_other_RNA_drtr_mfe_sum, fn_other_RNA_drtr_mfe_sum)
# print(mcc_other_RNA_drtr_mfe)

mcc_other_RNA_drtr_prob = mcc(tp_other_RNA_drtr_prob_sum, fp_other_RNA_drtr_prob_sum, tn_other_RNA_drtr_prob_sum, fn_other_RNA_drtr_prob_sum)
# print(mcc_other_RNA_drtr_prob)

mcc_other_RNA_drtr_bpp = mcc(tp_other_RNA_drtr_bpp_sum, fp_other_RNA_drtr_bpp_sum, tn_other_RNA_drtr_bpp_sum, fn_other_RNA_drtr_bpp_sum)
# print(mcc_other_RNA_drtr_bpp)


#OtherrRNA
other_rRNA=family_data.loc[family_data.family=='OtherrRNA']

tp_other_rRNA_fold_mfe = (other_rRNA[['RNAfold_mfe_tp']])
tp_other_rRNA_fold_mfe_sum = boot(tp_other_rRNA_fold_mfe)

fp_other_rRNA_fold_mfe = (other_rRNA[['RNAfold_mfe_fp']])
fp_other_rRNA_fold_mfe_sum = boot(fp_other_rRNA_fold_mfe)

tn_other_rRNA_fold_mfe = (other_rRNA[['RNAfold_mfe_tn']])
tn_other_rRNA_fold_mfe_sum = boot(tn_other_rRNA_fold_mfe)

fn_other_rRNA_fold_mfe = (other_rRNA[['RNAfold_mfe_fn']])
fn_other_rRNA_fold_mfe_sum = boot(fn_other_rRNA_fold_mfe)


tp_other_rRNA_fold_bpp = (other_rRNA[['RNAfold_bpp_tp']])
tp_other_rRNA_fold_bpp_sum = boot(tp_other_rRNA_fold_bpp)

fp_other_rRNA_fold_bpp = (other_rRNA[['RNAfold_bpp_fp']])
fp_other_rRNA_fold_bpp_sum = boot(fp_other_rRNA_fold_bpp)

tn_other_rRNA_fold_bpp = (other_rRNA[['RNAfold_bpp_tn']])
tn_other_rRNA_fold_bpp_sum = boot(tn_other_rRNA_fold_bpp)

fn_other_rRNA_fold_bpp = (other_rRNA[['RNAfold_bpp_fn']])
fn_other_rRNA_fold_bpp_sum = boot(fn_other_rRNA_fold_bpp)


tp_other_rRNA_drtr_mfe = (other_rRNA[['drtr_tp_mfe']])
tp_other_rRNA_drtr_mfe_sum = boot(tp_other_rRNA_drtr_mfe)

fp_other_rRNA_drtr_mfe = (other_rRNA[['drtr_fp_mfe']])
fp_other_rRNA_drtr_mfe_sum = boot(fp_other_rRNA_drtr_mfe)

tn_other_rRNA_drtr_mfe = (other_rRNA[['drtr_tn_mfe']])
tn_other_rRNA_drtr_mfe_sum = boot(tn_other_rRNA_drtr_mfe)

fn_other_rRNA_drtr_mfe = (other_rRNA[['drtr_fn_mfe']])
fn_other_rRNA_drtr_mfe_sum = boot(fn_other_rRNA_drtr_mfe)


tp_other_rRNA_drtr_prob = (other_rRNA[['drtr_tp_prob']])
tp_other_rRNA_drtr_prob_sum = boot(tp_other_rRNA_drtr_prob)

fp_other_rRNA_drtr_prob = (other_rRNA[['drtr_fp_prob']])
fp_other_rRNA_drtr_prob_sum = boot(fp_other_rRNA_drtr_prob)

tn_other_rRNA_drtr_prob = (other_rRNA[['drtr_tn_prob']])
tn_other_rRNA_drtr_prob_sum = boot(tn_other_rRNA_drtr_prob)

fn_other_rRNA_drtr_prob = (other_rRNA[['drtr_fn_prob']])
fn_other_rRNA_drtr_prob_sum = boot(fn_other_rRNA_drtr_prob)


tp_other_rRNA_drtr_bpp = (other_rRNA[['drtr_tp']])
tp_other_rRNA_drtr_bpp_sum = boot(tp_other_rRNA_drtr_bpp)

fp_other_rRNA_drtr_bpp = (other_rRNA[['drtr_fp']])
fp_other_rRNA_drtr_bpp_sum = boot(fp_other_rRNA_drtr_bpp)

tn_other_rRNA_drtr_bpp = (other_rRNA[['drtr_tn']])
tn_other_rRNA_drtr_bpp_sum = boot(tn_other_rRNA_drtr_bpp)

fn_other_rRNA_drtr_bpp = (other_rRNA[['drtr_fn']])
fn_other_rRNA_drtr_bpp_sum = boot(fn_other_rRNA_drtr_bpp)

#TPR
tpr_other_rRNA_fold_mfe = tpr(tp_other_rRNA_fold_mfe_sum, fn_other_rRNA_fold_mfe_sum)
# print(tpr_other_rRNA_fold_mfe)

tpr_other_rRNA_fold_bpp = tpr(tp_other_rRNA_fold_bpp_sum, fn_other_rRNA_fold_bpp_sum)
# print(tpr_other_rRNA_fold_bpp)

tpr_other_rRNA_drtr_mfe = tpr(tp_other_rRNA_drtr_mfe_sum, fn_other_rRNA_drtr_mfe_sum)
# print(tpr_other_rRNA_drtr_mfe)

tpr_other_rRNA_drtr_prob = tpr(tp_other_rRNA_drtr_prob_sum, fn_other_rRNA_drtr_prob_sum)
# print(tpr_other_rRNA_drtr_prob)

tpr_other_rRNA_drtr_bpp = tpr(tp_other_rRNA_drtr_bpp_sum, fn_other_rRNA_drtr_bpp_sum)
# print(tpr_other_rRNA_drtr_bpp)


#PPV
ppv_other_rRNA_fold_mfe = ppv(tp_other_rRNA_fold_mfe_sum, fn_other_rRNA_fold_mfe_sum)
# print(ppv_other_rRNA_fold_mfe)

ppv_other_rRNA_fold_bpp = ppv(tp_other_rRNA_fold_bpp_sum, fn_other_rRNA_fold_bpp_sum)
# print(ppv_other_rRNA_fold_bpp)

ppv_other_rRNA_drtr_mfe = ppv(tp_other_rRNA_drtr_mfe_sum, fn_other_rRNA_drtr_mfe_sum)
# print(ppv_other_rRNA_drtr_mfe)

ppv_other_rRNA_drtr_prob = ppv(tp_other_rRNA_drtr_prob_sum, fn_other_rRNA_drtr_prob_sum)
# print(ppv_other_rRNA_drtr_prob)

ppv_other_rRNA_drtr_bpp = ppv(tp_other_rRNA_drtr_bpp_sum, fn_other_rRNA_drtr_bpp_sum)
# print(ppv_other_rRNA_drtr_bpp)

#F1
f_1_other_rRNA_fold_mfe = f_1(tp_other_rRNA_fold_mfe_sum, fp_other_rRNA_fold_mfe_sum, fn_other_rRNA_fold_mfe_sum)
# print(f_1_other_rRNA_fold_mfe)

f_1_other_rRNA_fold_bpp = f_1(tp_other_rRNA_fold_bpp_sum, fp_other_rRNA_fold_bpp_sum, fn_other_rRNA_fold_bpp_sum)
# print(f_1_other_rRNA_fold_bpp)

f_1_other_rRNA_drtr_mfe = f_1(tp_other_rRNA_drtr_mfe_sum, fp_other_rRNA_drtr_mfe_sum, fn_other_rRNA_drtr_mfe_sum)
# print(f_1_other_rRNA_drtr_mfe)

f_1_other_rRNA_drtr_prob = f_1(tp_other_rRNA_drtr_prob_sum, fp_other_rRNA_drtr_prob_sum, fn_other_rRNA_drtr_prob_sum)
# print(f_1_other_rRNA_drtr_prob)

f_1_other_rRNA_drtr_bpp = f_1(tp_other_rRNA_drtr_bpp_sum, fp_other_rRNA_drtr_bpp_sum, fn_other_rRNA_drtr_bpp_sum)
# print(f_1_other_rRNA_drtr_bpp)

# #MCC
mcc_other_rRNA_fold_mfe = mcc(tp_other_rRNA_fold_mfe_sum, fp_other_rRNA_fold_mfe_sum, tn_other_rRNA_fold_mfe_sum, fn_other_rRNA_fold_mfe_sum)
# print(mcc_other_rRNA_fold_mfe)

mcc_other_rRNA_fold_bpp = mcc(tp_other_rRNA_fold_bpp_sum, fp_other_rRNA_fold_bpp_sum, tn_other_rRNA_fold_bpp_sum, fn_other_rRNA_fold_bpp_sum)
# print(mcc_other_rRNA_fold_bpp)

mcc_other_rRNA_drtr_mfe = mcc(tp_other_rRNA_drtr_mfe_sum, fp_other_rRNA_drtr_mfe_sum, tn_other_rRNA_drtr_mfe_sum, fn_other_rRNA_drtr_mfe_sum)
# print(mcc_other_rRNA_drtr_mfe)

mcc_other_rRNA_drtr_prob = mcc(tp_other_rRNA_drtr_prob_sum, fp_other_rRNA_drtr_prob_sum, tn_other_rRNA_drtr_prob_sum, fn_other_rRNA_drtr_prob_sum)
# print(mcc_other_rRNA_drtr_prob)

mcc_other_rRNA_drtr_bpp = mcc(tp_other_rRNA_drtr_bpp_sum, fp_other_rRNA_drtr_bpp_sum, tn_other_rRNA_drtr_bpp_sum, fn_other_rRNA_drtr_bpp_sum)
# print(mcc_other_rRNA_drtr_bpp)



#RNAIII
RNAIII=family_data.loc[family_data.family=='RNAIII']

tp_RNAIII_fold_mfe = (RNAIII[['RNAfold_mfe_tp']])
tp_RNAIII_fold_mfe_sum = boot(tp_RNAIII_fold_mfe)

fp_RNAIII_fold_mfe = (RNAIII[['RNAfold_mfe_fp']])
fp_RNAIII_fold_mfe_sum = boot(fp_RNAIII_fold_mfe)

tn_RNAIII_fold_mfe = (RNAIII[['RNAfold_mfe_tn']])
tn_RNAIII_fold_mfe_sum = boot(tn_RNAIII_fold_mfe)

fn_RNAIII_fold_mfe = (RNAIII[['RNAfold_mfe_fn']])
fn_RNAIII_fold_mfe_sum = boot(fn_RNAIII_fold_mfe)


tp_RNAIII_fold_bpp = (RNAIII[['RNAfold_bpp_tp']])
tp_RNAIII_fold_bpp_sum = boot(tp_RNAIII_fold_bpp)

fp_RNAIII_fold_bpp = (RNAIII[['RNAfold_bpp_fp']])
fp_RNAIII_fold_bpp_sum = boot(fp_RNAIII_fold_bpp)

tn_RNAIII_fold_bpp = (RNAIII[['RNAfold_bpp_tn']])
tn_RNAIII_fold_bpp_sum = boot(tn_RNAIII_fold_bpp)

fn_RNAIII_fold_bpp = (RNAIII[['RNAfold_bpp_fn']])
fn_RNAIII_fold_bpp_sum = boot(fn_RNAIII_fold_bpp)


tp_RNAIII_drtr_mfe = (RNAIII[['drtr_tp_mfe']])
tp_RNAIII_drtr_mfe_sum = boot(tp_RNAIII_drtr_mfe)

fp_RNAIII_drtr_mfe = (RNAIII[['drtr_fp_mfe']])
fp_RNAIII_drtr_mfe_sum = boot(fp_RNAIII_drtr_mfe)

tn_RNAIII_drtr_mfe = (RNAIII[['drtr_tn_mfe']])
tn_RNAIII_drtr_mfe_sum = boot(tn_RNAIII_drtr_mfe)

fn_RNAIII_drtr_mfe = (RNAIII[['drtr_fn_mfe']])
fn_RNAIII_drtr_mfe_sum = boot(fn_RNAIII_drtr_mfe)


tp_RNAIII_drtr_prob = (RNAIII[['drtr_tp_prob']])
tp_RNAIII_drtr_prob_sum = boot(tp_RNAIII_drtr_prob)

fp_RNAIII_drtr_prob = (RNAIII[['drtr_fp_prob']])
fp_RNAIII_drtr_prob_sum = boot(fp_RNAIII_drtr_prob)

tn_RNAIII_drtr_prob = (RNAIII[['drtr_tn_prob']])
tn_RNAIII_drtr_prob_sum = boot(tn_RNAIII_drtr_prob)

fn_RNAIII_drtr_prob = (RNAIII[['drtr_fn_prob']])
fn_RNAIII_drtr_prob_sum = boot(fn_RNAIII_drtr_prob)


tp_RNAIII_drtr_bpp = (RNAIII[['drtr_tp']])
tp_RNAIII_drtr_bpp_sum = boot(tp_RNAIII_drtr_bpp)

fp_RNAIII_drtr_bpp = (RNAIII[['drtr_fp']])
fp_RNAIII_drtr_bpp_sum = boot(fp_RNAIII_drtr_bpp)

tn_RNAIII_drtr_bpp = (RNAIII[['drtr_tn']])
tn_RNAIII_drtr_bpp_sum = boot(tn_RNAIII_drtr_bpp)

fn_RNAIII_drtr_bpp = (RNAIII[['drtr_fn']])
fn_RNAIII_drtr_bpp_sum = boot(fn_RNAIII_drtr_bpp)

#TPR
tpr_RNAIII_fold_mfe = tpr(tp_RNAIII_fold_mfe_sum, fn_RNAIII_fold_mfe_sum)
# print(tpr_RNAIII_fold_mfe)

tpr_RNAIII_fold_bpp = tpr(tp_RNAIII_fold_bpp_sum, fn_RNAIII_fold_bpp_sum)
# print(tpr_RNAIII_fold_bpp)

tpr_RNAIII_drtr_mfe = tpr(tp_RNAIII_drtr_mfe_sum, fn_RNAIII_drtr_mfe_sum)
# print(tpr_RNAIII_drtr_mfe)

tpr_RNAIII_drtr_prob = tpr(tp_RNAIII_drtr_prob_sum, fn_RNAIII_drtr_prob_sum)
# print(tpr_RNAIII_drtr_prob)

tpr_RNAIII_drtr_bpp = tpr(tp_RNAIII_drtr_bpp_sum, fn_RNAIII_drtr_bpp_sum)
# print(tpr_RNAIII_drtr_bpp)


#PPV
ppv_RNAIII_fold_mfe = ppv(tp_RNAIII_fold_mfe_sum, fn_RNAIII_fold_mfe_sum)
# print(ppv_RNAIII_fold_mfe)

ppv_RNAIII_fold_bpp = ppv(tp_RNAIII_fold_bpp_sum, fn_RNAIII_fold_bpp_sum)
# print(ppv_RNAIII_fold_bpp)

ppv_RNAIII_drtr_mfe = ppv(tp_RNAIII_drtr_mfe_sum, fn_RNAIII_drtr_mfe_sum)
# print(ppv_RNAIII_drtr_mfe)

ppv_RNAIII_drtr_prob = ppv(tp_RNAIII_drtr_prob_sum, fn_RNAIII_drtr_prob_sum)
# print(ppv_RNAIII_drtr_prob)

ppv_RNAIII_drtr_bpp = ppv(tp_RNAIII_drtr_bpp_sum, fn_RNAIII_drtr_bpp_sum)
# print(ppv_RNAIII_drtr_bpp)

#F1
f_1_RNAIII_fold_mfe = f_1(tp_RNAIII_fold_mfe_sum, fp_RNAIII_fold_mfe_sum, fn_RNAIII_fold_mfe_sum)
# print(f_1_RNAIII_fold_mfe)

f_1_RNAIII_fold_bpp = f_1(tp_RNAIII_fold_bpp_sum, fp_RNAIII_fold_bpp_sum, fn_RNAIII_fold_bpp_sum)
# print(f_1_RNAIII_fold_bpp)

f_1_RNAIII_drtr_mfe = f_1(tp_RNAIII_drtr_mfe_sum, fp_RNAIII_drtr_mfe_sum, fn_RNAIII_drtr_mfe_sum)
# print(f_1_RNAIII_drtr_mfe)

f_1_RNAIII_drtr_prob = f_1(tp_RNAIII_drtr_prob_sum, fp_RNAIII_drtr_prob_sum, fn_RNAIII_drtr_prob_sum)
# print(f_1_RNAIII_drtr_prob)

f_1_RNAIII_drtr_bpp = f_1(tp_RNAIII_drtr_bpp_sum, fp_RNAIII_drtr_bpp_sum, fn_RNAIII_drtr_bpp_sum)
# print(f_1_RNAIII_drtr_bpp)

# #MCC
mcc_RNAIII_fold_mfe = mcc(tp_RNAIII_fold_mfe_sum, fp_RNAIII_fold_mfe_sum, tn_RNAIII_fold_mfe_sum, fn_RNAIII_fold_mfe_sum)
# print(mcc_RNAIII_fold_mfe)

mcc_RNAIII_fold_bpp = mcc(tp_RNAIII_fold_bpp_sum, fp_RNAIII_fold_bpp_sum, tn_RNAIII_fold_bpp_sum, fn_RNAIII_fold_bpp_sum)
# print(mcc_RNAIII_fold_bpp)

mcc_RNAIII_drtr_mfe = mcc(tp_RNAIII_drtr_mfe_sum, fp_RNAIII_drtr_mfe_sum, tn_RNAIII_drtr_mfe_sum, fn_RNAIII_drtr_mfe_sum)
# print(mcc_RNAIII_drtr_mfe)

mcc_RNAIII_drtr_prob = mcc(tp_RNAIII_drtr_prob_sum, fp_RNAIII_drtr_prob_sum, tn_RNAIII_drtr_prob_sum, fn_RNAIII_drtr_prob_sum)
# print(mcc_RNAIII_drtr_prob)

mcc_RNAIII_drtr_bpp = mcc(tp_RNAIII_drtr_bpp_sum, fp_RNAIII_drtr_bpp_sum, tn_RNAIII_drtr_bpp_sum, fn_RNAIII_drtr_bpp_sum)
# print(mcc_RNAIII_drtr_bpp)


#RNaseE5UTR
RNaseE5UTR=family_data.loc[family_data.family=='RNaseE5UTR']

tp_RNaseE5UTR_fold_mfe = (RNaseE5UTR[['RNAfold_mfe_tp']])
tp_RNaseE5UTR_fold_mfe_sum = boot(tp_RNaseE5UTR_fold_mfe)

fp_RNaseE5UTR_fold_mfe = (RNaseE5UTR[['RNAfold_mfe_fp']])
fp_RNaseE5UTR_fold_mfe_sum = boot(fp_RNaseE5UTR_fold_mfe)

tn_RNaseE5UTR_fold_mfe = (RNaseE5UTR[['RNAfold_mfe_tn']])
tn_RNaseE5UTR_fold_mfe_sum = boot(tn_RNaseE5UTR_fold_mfe)

fn_RNaseE5UTR_fold_mfe = (RNaseE5UTR[['RNAfold_mfe_fn']])
fn_RNaseE5UTR_fold_mfe_sum = boot(fn_RNaseE5UTR_fold_mfe)


tp_RNaseE5UTR_fold_bpp = (RNaseE5UTR[['RNAfold_bpp_tp']])
tp_RNaseE5UTR_fold_bpp_sum = boot(tp_RNaseE5UTR_fold_bpp)

fp_RNaseE5UTR_fold_bpp = (RNaseE5UTR[['RNAfold_bpp_fp']])
fp_RNaseE5UTR_fold_bpp_sum = boot(fp_RNaseE5UTR_fold_bpp)

tn_RNaseE5UTR_fold_bpp = (RNaseE5UTR[['RNAfold_bpp_tn']])
tn_RNaseE5UTR_fold_bpp_sum = boot(tn_RNaseE5UTR_fold_bpp)

fn_RNaseE5UTR_fold_bpp = (RNaseE5UTR[['RNAfold_bpp_fn']])
fn_RNaseE5UTR_fold_bpp_sum = boot(fn_RNaseE5UTR_fold_bpp)


tp_RNaseE5UTR_drtr_mfe = (RNaseE5UTR[['drtr_tp_mfe']])
tp_RNaseE5UTR_drtr_mfe_sum = boot(tp_RNaseE5UTR_drtr_mfe)

fp_RNaseE5UTR_drtr_mfe = (RNaseE5UTR[['drtr_fp_mfe']])
fp_RNaseE5UTR_drtr_mfe_sum = boot(fp_RNaseE5UTR_drtr_mfe)

tn_RNaseE5UTR_drtr_mfe = (RNaseE5UTR[['drtr_tn_mfe']])
tn_RNaseE5UTR_drtr_mfe_sum = boot(tn_RNaseE5UTR_drtr_mfe)

fn_RNaseE5UTR_drtr_mfe = (RNaseE5UTR[['drtr_fn_mfe']])
fn_RNaseE5UTR_drtr_mfe_sum = boot(fn_RNaseE5UTR_drtr_mfe)


tp_RNaseE5UTR_drtr_prob = (RNaseE5UTR[['drtr_tp_prob']])
tp_RNaseE5UTR_drtr_prob_sum = boot(tp_RNaseE5UTR_drtr_prob)

fp_RNaseE5UTR_drtr_prob = (RNaseE5UTR[['drtr_fp_prob']])
fp_RNaseE5UTR_drtr_prob_sum = boot(fp_RNaseE5UTR_drtr_prob)

tn_RNaseE5UTR_drtr_prob = (RNaseE5UTR[['drtr_tn_prob']])
tn_RNaseE5UTR_drtr_prob_sum = boot(tn_RNaseE5UTR_drtr_prob)

fn_RNaseE5UTR_drtr_prob = (RNaseE5UTR[['drtr_fn_prob']])
fn_RNaseE5UTR_drtr_prob_sum = boot(fn_RNaseE5UTR_drtr_prob)


tp_RNaseE5UTR_drtr_bpp = (RNaseE5UTR[['drtr_tp']])
tp_RNaseE5UTR_drtr_bpp_sum = boot(tp_RNaseE5UTR_drtr_bpp)

fp_RNaseE5UTR_drtr_bpp = (RNaseE5UTR[['drtr_fp']])
fp_RNaseE5UTR_drtr_bpp_sum = boot(fp_RNaseE5UTR_drtr_bpp)

tn_RNaseE5UTR_drtr_bpp = (RNaseE5UTR[['drtr_tn']])
tn_RNaseE5UTR_drtr_bpp_sum = boot(tn_RNaseE5UTR_drtr_bpp)

fn_RNaseE5UTR_drtr_bpp = (RNaseE5UTR[['drtr_fn']])
fn_RNaseE5UTR_drtr_bpp_sum = boot(fn_RNaseE5UTR_drtr_bpp)

#TPR
tpr_RNaseE5UTR_fold_mfe = tpr(tp_RNaseE5UTR_fold_mfe_sum, fn_RNaseE5UTR_fold_mfe_sum)
# print(tpr_RNaseE5UTR_fold_mfe)

tpr_RNaseE5UTR_fold_bpp = tpr(tp_RNaseE5UTR_fold_bpp_sum, fn_RNaseE5UTR_fold_bpp_sum)
# print(tpr_RNaseE5UTR_fold_bpp)

tpr_RNaseE5UTR_drtr_mfe = tpr(tp_RNaseE5UTR_drtr_mfe_sum, fn_RNaseE5UTR_drtr_mfe_sum)
# print(tpr_RNaseE5UTR_drtr_mfe)

tpr_RNaseE5UTR_drtr_prob = tpr(tp_RNaseE5UTR_drtr_prob_sum, fn_RNaseE5UTR_drtr_prob_sum)
# print(tpr_RNaseE5UTR_drtr_prob)

tpr_RNaseE5UTR_drtr_bpp = tpr(tp_RNaseE5UTR_drtr_bpp_sum, fn_RNaseE5UTR_drtr_bpp_sum)
# print(tpr_RNaseE5UTR_drtr_bpp)


#PPV
ppv_RNaseE5UTR_fold_mfe = ppv(tp_RNaseE5UTR_fold_mfe_sum, fn_RNaseE5UTR_fold_mfe_sum)
# print(ppv_RNaseE5UTR_fold_mfe)

ppv_RNaseE5UTR_fold_bpp = ppv(tp_RNaseE5UTR_fold_bpp_sum, fn_RNaseE5UTR_fold_bpp_sum)
# print(ppv_RNaseE5UTR_fold_bpp)

ppv_RNaseE5UTR_drtr_mfe = ppv(tp_RNaseE5UTR_drtr_mfe_sum, fn_RNaseE5UTR_drtr_mfe_sum)
# print(ppv_RNaseE5UTR_drtr_mfe)

ppv_RNaseE5UTR_drtr_prob = ppv(tp_RNaseE5UTR_drtr_prob_sum, fn_RNaseE5UTR_drtr_prob_sum)
# print(ppv_RNaseE5UTR_drtr_prob)

ppv_RNaseE5UTR_drtr_bpp = ppv(tp_RNaseE5UTR_drtr_bpp_sum, fn_RNaseE5UTR_drtr_bpp_sum)
# print(ppv_RNaseE5UTR_drtr_bpp)

#F1
f_1_RNaseE5UTR_fold_mfe = f_1(tp_RNaseE5UTR_fold_mfe_sum, fp_RNaseE5UTR_fold_mfe_sum, fn_RNaseE5UTR_fold_mfe_sum)
# print(f_1_RNaseE5UTR_fold_mfe)

f_1_RNaseE5UTR_fold_bpp = f_1(tp_RNaseE5UTR_fold_bpp_sum, fp_RNaseE5UTR_fold_bpp_sum, fn_RNaseE5UTR_fold_bpp_sum)
# print(f_1_RNaseE5UTR_fold_bpp)

f_1_RNaseE5UTR_drtr_mfe = f_1(tp_RNaseE5UTR_drtr_mfe_sum, fp_RNaseE5UTR_drtr_mfe_sum, fn_RNaseE5UTR_drtr_mfe_sum)
# print(f_1_RNaseE5UTR_drtr_mfe)

f_1_RNaseE5UTR_drtr_prob = f_1(tp_RNaseE5UTR_drtr_prob_sum, fp_RNaseE5UTR_drtr_prob_sum, fn_RNaseE5UTR_drtr_prob_sum)
# print(f_1_RNaseE5UTR_drtr_prob)

f_1_RNaseE5UTR_drtr_bpp = f_1(tp_RNaseE5UTR_drtr_bpp_sum, fp_RNaseE5UTR_drtr_bpp_sum, fn_RNaseE5UTR_drtr_bpp_sum)
# print(f_1_RNaseE5UTR_drtr_bpp)

# #MCC
mcc_RNaseE5UTR_fold_mfe = mcc(tp_RNaseE5UTR_fold_mfe_sum, fp_RNaseE5UTR_fold_mfe_sum, tn_RNaseE5UTR_fold_mfe_sum, fn_RNaseE5UTR_fold_mfe_sum)
# print(mcc_RNaseE5UTR_fold_mfe)

mcc_RNaseE5UTR_fold_bpp = mcc(tp_RNaseE5UTR_fold_bpp_sum, fp_RNaseE5UTR_fold_bpp_sum, tn_RNaseE5UTR_fold_bpp_sum, fn_RNaseE5UTR_fold_bpp_sum)
# print(mcc_RNaseE5UTR_fold_bpp)

mcc_RNaseE5UTR_drtr_mfe = mcc(tp_RNaseE5UTR_drtr_mfe_sum, fp_RNaseE5UTR_drtr_mfe_sum, tn_RNaseE5UTR_drtr_mfe_sum, fn_RNaseE5UTR_drtr_mfe_sum)
# print(mcc_RNaseE5UTR_drtr_mfe)

mcc_RNaseE5UTR_drtr_prob = mcc(tp_RNaseE5UTR_drtr_prob_sum, fp_RNaseE5UTR_drtr_prob_sum, tn_RNaseE5UTR_drtr_prob_sum, fn_RNaseE5UTR_drtr_prob_sum)
# print(mcc_RNaseE5UTR_drtr_prob)

mcc_RNaseE5UTR_drtr_bpp = mcc(tp_RNaseE5UTR_drtr_bpp_sum, fp_RNaseE5UTR_drtr_bpp_sum, tn_RNaseE5UTR_drtr_bpp_sum, fn_RNaseE5UTR_drtr_bpp_sum)
# print(mcc_RNaseE5UTR_drtr_bpp)




#RNaseMRPRNA
RNaseMRPRNA=family_data.loc[family_data.family=='RNaseMRPRNA']

tp_RNaseMRPRNA_fold_mfe = (RNaseMRPRNA[['RNAfold_mfe_tp']])
tp_RNaseMRPRNA_fold_mfe_sum = boot(tp_RNaseMRPRNA_fold_mfe)

fp_RNaseMRPRNA_fold_mfe = (RNaseMRPRNA[['RNAfold_mfe_fp']])
fp_RNaseMRPRNA_fold_mfe_sum = boot(fp_RNaseMRPRNA_fold_mfe)

tn_RNaseMRPRNA_fold_mfe = (RNaseMRPRNA[['RNAfold_mfe_tn']])
tn_RNaseMRPRNA_fold_mfe_sum = boot(tn_RNaseMRPRNA_fold_mfe)

fn_RNaseMRPRNA_fold_mfe = (RNaseMRPRNA[['RNAfold_mfe_fn']])
fn_RNaseMRPRNA_fold_mfe_sum = boot(fn_RNaseMRPRNA_fold_mfe)


tp_RNaseMRPRNA_fold_bpp = (RNaseMRPRNA[['RNAfold_bpp_tp']])
tp_RNaseMRPRNA_fold_bpp_sum = boot(tp_RNaseMRPRNA_fold_bpp)

fp_RNaseMRPRNA_fold_bpp = (RNaseMRPRNA[['RNAfold_bpp_fp']])
fp_RNaseMRPRNA_fold_bpp_sum = boot(fp_RNaseMRPRNA_fold_bpp)

tn_RNaseMRPRNA_fold_bpp = (RNaseMRPRNA[['RNAfold_bpp_tn']])
tn_RNaseMRPRNA_fold_bpp_sum = boot(tn_RNaseMRPRNA_fold_bpp)

fn_RNaseMRPRNA_fold_bpp = (RNaseMRPRNA[['RNAfold_bpp_fn']])
fn_RNaseMRPRNA_fold_bpp_sum = boot(fn_RNaseMRPRNA_fold_bpp)


tp_RNaseMRPRNA_drtr_mfe = (RNaseMRPRNA[['drtr_tp_mfe']])
tp_RNaseMRPRNA_drtr_mfe_sum = boot(tp_RNaseMRPRNA_drtr_mfe)

fp_RNaseMRPRNA_drtr_mfe = (RNaseMRPRNA[['drtr_fp_mfe']])
fp_RNaseMRPRNA_drtr_mfe_sum = boot(fp_RNaseMRPRNA_drtr_mfe)

tn_RNaseMRPRNA_drtr_mfe = (RNaseMRPRNA[['drtr_tn_mfe']])
tn_RNaseMRPRNA_drtr_mfe_sum = boot(tn_RNaseMRPRNA_drtr_mfe)

fn_RNaseMRPRNA_drtr_mfe = (RNaseMRPRNA[['drtr_fn_mfe']])
fn_RNaseMRPRNA_drtr_mfe_sum = boot(fn_RNaseMRPRNA_drtr_mfe)


tp_RNaseMRPRNA_drtr_prob = (RNaseMRPRNA[['drtr_tp_prob']])
tp_RNaseMRPRNA_drtr_prob_sum = boot(tp_RNaseMRPRNA_drtr_prob)

fp_RNaseMRPRNA_drtr_prob = (RNaseMRPRNA[['drtr_fp_prob']])
fp_RNaseMRPRNA_drtr_prob_sum = boot(fp_RNaseMRPRNA_drtr_prob)

tn_RNaseMRPRNA_drtr_prob = (RNaseMRPRNA[['drtr_tn_prob']])
tn_RNaseMRPRNA_drtr_prob_sum = boot(tn_RNaseMRPRNA_drtr_prob)

fn_RNaseMRPRNA_drtr_prob = (RNaseMRPRNA[['drtr_fn_prob']])
fn_RNaseMRPRNA_drtr_prob_sum = boot(fn_RNaseMRPRNA_drtr_prob)


tp_RNaseMRPRNA_drtr_bpp = (RNaseMRPRNA[['drtr_tp']])
tp_RNaseMRPRNA_drtr_bpp_sum = boot(tp_RNaseMRPRNA_drtr_bpp)

fp_RNaseMRPRNA_drtr_bpp = (RNaseMRPRNA[['drtr_fp']])
fp_RNaseMRPRNA_drtr_bpp_sum = boot(fp_RNaseMRPRNA_drtr_bpp)

tn_RNaseMRPRNA_drtr_bpp = (RNaseMRPRNA[['drtr_tn']])
tn_RNaseMRPRNA_drtr_bpp_sum = boot(tn_RNaseMRPRNA_drtr_bpp)

fn_RNaseMRPRNA_drtr_bpp = (RNaseMRPRNA[['drtr_fn']])
fn_RNaseMRPRNA_drtr_bpp_sum = boot(fn_RNaseMRPRNA_drtr_bpp)

#TPR
tpr_RNaseMRPRNA_fold_mfe = tpr(tp_RNaseMRPRNA_fold_mfe_sum, fn_RNaseMRPRNA_fold_mfe_sum)
# print(tpr_RNaseMRPRNA_fold_mfe)

tpr_RNaseMRPRNA_fold_bpp = tpr(tp_RNaseMRPRNA_fold_bpp_sum, fn_RNaseMRPRNA_fold_bpp_sum)
# print(tpr_RNaseMRPRNA_fold_bpp)

tpr_RNaseMRPRNA_drtr_mfe = tpr(tp_RNaseMRPRNA_drtr_mfe_sum, fn_RNaseMRPRNA_drtr_mfe_sum)
# print(tpr_RNaseMRPRNA_drtr_mfe)

tpr_RNaseMRPRNA_drtr_prob = tpr(tp_RNaseMRPRNA_drtr_prob_sum, fn_RNaseMRPRNA_drtr_prob_sum)
# print(tpr_RNaseMRPRNA_drtr_prob)

tpr_RNaseMRPRNA_drtr_bpp = tpr(tp_RNaseMRPRNA_drtr_bpp_sum, fn_RNaseMRPRNA_drtr_bpp_sum)
# print(tpr_RNaseMRPRNA_drtr_bpp)


#PPV
ppv_RNaseMRPRNA_fold_mfe = ppv(tp_RNaseMRPRNA_fold_mfe_sum, fn_RNaseMRPRNA_fold_mfe_sum)
# print(ppv_RNaseMRPRNA_fold_mfe)

ppv_RNaseMRPRNA_fold_bpp = ppv(tp_RNaseMRPRNA_fold_bpp_sum, fn_RNaseMRPRNA_fold_bpp_sum)
# print(ppv_RNaseMRPRNA_fold_bpp)

ppv_RNaseMRPRNA_drtr_mfe = ppv(tp_RNaseMRPRNA_drtr_mfe_sum, fn_RNaseMRPRNA_drtr_mfe_sum)
# print(ppv_RNaseMRPRNA_drtr_mfe)

ppv_RNaseMRPRNA_drtr_prob = ppv(tp_RNaseMRPRNA_drtr_prob_sum, fn_RNaseMRPRNA_drtr_prob_sum)
# print(ppv_RNaseMRPRNA_drtr_prob)

ppv_RNaseMRPRNA_drtr_bpp = ppv(tp_RNaseMRPRNA_drtr_bpp_sum, fn_RNaseMRPRNA_drtr_bpp_sum)
# print(ppv_RNaseMRPRNA_drtr_bpp)

#F1
f_1_RNaseMRPRNA_fold_mfe = f_1(tp_RNaseMRPRNA_fold_mfe_sum, fp_RNaseMRPRNA_fold_mfe_sum, fn_RNaseMRPRNA_fold_mfe_sum)
# print(f_1_RNaseMRPRNA_fold_mfe)

f_1_RNaseMRPRNA_fold_bpp = f_1(tp_RNaseMRPRNA_fold_bpp_sum, fp_RNaseMRPRNA_fold_bpp_sum, fn_RNaseMRPRNA_fold_bpp_sum)
# print(f_1_RNaseMRPRNA_fold_bpp)

f_1_RNaseMRPRNA_drtr_mfe = f_1(tp_RNaseMRPRNA_drtr_mfe_sum, fp_RNaseMRPRNA_drtr_mfe_sum, fn_RNaseMRPRNA_drtr_mfe_sum)
# print(f_1_RNaseMRPRNA_drtr_mfe)

f_1_RNaseMRPRNA_drtr_prob = f_1(tp_RNaseMRPRNA_drtr_prob_sum, fp_RNaseMRPRNA_drtr_prob_sum, fn_RNaseMRPRNA_drtr_prob_sum)
# print(f_1_RNaseMRPRNA_drtr_prob)

f_1_RNaseMRPRNA_drtr_bpp = f_1(tp_RNaseMRPRNA_drtr_bpp_sum, fp_RNaseMRPRNA_drtr_bpp_sum, fn_RNaseMRPRNA_drtr_bpp_sum)
# print(f_1_RNaseMRPRNA_drtr_bpp)

# #MCC
mcc_RNaseMRPRNA_fold_mfe = mcc(tp_RNaseMRPRNA_fold_mfe_sum, fp_RNaseMRPRNA_fold_mfe_sum, tn_RNaseMRPRNA_fold_mfe_sum, fn_RNaseMRPRNA_fold_mfe_sum)
# print(mcc_RNaseMRPRNA_fold_mfe)

mcc_RNaseMRPRNA_fold_bpp = mcc(tp_RNaseMRPRNA_fold_bpp_sum, fp_RNaseMRPRNA_fold_bpp_sum, tn_RNaseMRPRNA_fold_bpp_sum, fn_RNaseMRPRNA_fold_bpp_sum)
# print(mcc_RNaseMRPRNA_fold_bpp)

mcc_RNaseMRPRNA_drtr_mfe = mcc(tp_RNaseMRPRNA_drtr_mfe_sum, fp_RNaseMRPRNA_drtr_mfe_sum, tn_RNaseMRPRNA_drtr_mfe_sum, fn_RNaseMRPRNA_drtr_mfe_sum)
# print(mcc_RNaseMRPRNA_drtr_mfe)

mcc_RNaseMRPRNA_drtr_prob = mcc(tp_RNaseMRPRNA_drtr_prob_sum, fp_RNaseMRPRNA_drtr_prob_sum, tn_RNaseMRPRNA_drtr_prob_sum, fn_RNaseMRPRNA_drtr_prob_sum)
# print(mcc_RNaseMRPRNA_drtr_prob)

mcc_RNaseMRPRNA_drtr_bpp = mcc(tp_RNaseMRPRNA_drtr_bpp_sum, fp_RNaseMRPRNA_drtr_bpp_sum, tn_RNaseMRPRNA_drtr_bpp_sum, fn_RNaseMRPRNA_drtr_bpp_sum)
# print(mcc_RNaseMRPRNA_drtr_bpp)


#RNasePRNA
RNasePRNA=family_data.loc[family_data.family=='RNasePRNA']

tp_RNasePRNA_fold_mfe = (RNasePRNA[['RNAfold_mfe_tp']])
tp_RNasePRNA_fold_mfe_sum = boot(tp_RNasePRNA_fold_mfe)

fp_RNasePRNA_fold_mfe = (RNasePRNA[['RNAfold_mfe_fp']])
fp_RNasePRNA_fold_mfe_sum = boot(fp_RNasePRNA_fold_mfe)

tn_RNasePRNA_fold_mfe = (RNasePRNA[['RNAfold_mfe_tn']])
tn_RNasePRNA_fold_mfe_sum = boot(tn_RNasePRNA_fold_mfe)

fn_RNasePRNA_fold_mfe = (RNasePRNA[['RNAfold_mfe_fn']])
fn_RNasePRNA_fold_mfe_sum = boot(fn_RNasePRNA_fold_mfe)


tp_RNasePRNA_fold_bpp = (RNasePRNA[['RNAfold_bpp_tp']])
tp_RNasePRNA_fold_bpp_sum = boot(tp_RNasePRNA_fold_bpp)

fp_RNasePRNA_fold_bpp = (RNasePRNA[['RNAfold_bpp_fp']])
fp_RNasePRNA_fold_bpp_sum = boot(fp_RNasePRNA_fold_bpp)

tn_RNasePRNA_fold_bpp = (RNasePRNA[['RNAfold_bpp_tn']])
tn_RNasePRNA_fold_bpp_sum = boot(tn_RNasePRNA_fold_bpp)

fn_RNasePRNA_fold_bpp = (RNasePRNA[['RNAfold_bpp_fn']])
fn_RNasePRNA_fold_bpp_sum = boot(fn_RNasePRNA_fold_bpp)


tp_RNasePRNA_drtr_mfe = (RNasePRNA[['drtr_tp_mfe']])
tp_RNasePRNA_drtr_mfe_sum = boot(tp_RNasePRNA_drtr_mfe)

fp_RNasePRNA_drtr_mfe = (RNasePRNA[['drtr_fp_mfe']])
fp_RNasePRNA_drtr_mfe_sum = boot(fp_RNasePRNA_drtr_mfe)

tn_RNasePRNA_drtr_mfe = (RNasePRNA[['drtr_tn_mfe']])
tn_RNasePRNA_drtr_mfe_sum = boot(tn_RNasePRNA_drtr_mfe)

fn_RNasePRNA_drtr_mfe = (RNasePRNA[['drtr_fn_mfe']])
fn_RNasePRNA_drtr_mfe_sum = boot(fn_RNasePRNA_drtr_mfe)


tp_RNasePRNA_drtr_prob = (RNasePRNA[['drtr_tp_prob']])
tp_RNasePRNA_drtr_prob_sum = boot(tp_RNasePRNA_drtr_prob)

fp_RNasePRNA_drtr_prob = (RNasePRNA[['drtr_fp_prob']])
fp_RNasePRNA_drtr_prob_sum = boot(fp_RNasePRNA_drtr_prob)

tn_RNasePRNA_drtr_prob = (RNasePRNA[['drtr_tn_prob']])
tn_RNasePRNA_drtr_prob_sum = boot(tn_RNasePRNA_drtr_prob)

fn_RNasePRNA_drtr_prob = (RNasePRNA[['drtr_fn_prob']])
fn_RNasePRNA_drtr_prob_sum = boot(fn_RNasePRNA_drtr_prob)


tp_RNasePRNA_drtr_bpp = (RNasePRNA[['drtr_tp']])
tp_RNasePRNA_drtr_bpp_sum = boot(tp_RNasePRNA_drtr_bpp)

fp_RNasePRNA_drtr_bpp = (RNasePRNA[['drtr_fp']])
fp_RNasePRNA_drtr_bpp_sum = boot(fp_RNasePRNA_drtr_bpp)

tn_RNasePRNA_drtr_bpp = (RNasePRNA[['drtr_tn']])
tn_RNasePRNA_drtr_bpp_sum = boot(tn_RNasePRNA_drtr_bpp)

fn_RNasePRNA_drtr_bpp = (RNasePRNA[['drtr_fn']])
fn_RNasePRNA_drtr_bpp_sum = boot(fn_RNasePRNA_drtr_bpp)

#TPR
tpr_RNasePRNA_fold_mfe = tpr(tp_RNasePRNA_fold_mfe_sum, fn_RNasePRNA_fold_mfe_sum)
# print(tpr_RNasePRNA_fold_mfe)

tpr_RNasePRNA_fold_bpp = tpr(tp_RNasePRNA_fold_bpp_sum, fn_RNasePRNA_fold_bpp_sum)
# print(tpr_RNasePRNA_fold_bpp)

tpr_RNasePRNA_drtr_mfe = tpr(tp_RNasePRNA_drtr_mfe_sum, fn_RNasePRNA_drtr_mfe_sum)
# print(tpr_RNasePRNA_drtr_mfe)

tpr_RNasePRNA_drtr_prob = tpr(tp_RNasePRNA_drtr_prob_sum, fn_RNasePRNA_drtr_prob_sum)
# print(tpr_RNasePRNA_drtr_prob)

tpr_RNasePRNA_drtr_bpp = tpr(tp_RNasePRNA_drtr_bpp_sum, fn_RNasePRNA_drtr_bpp_sum)
# print(tpr_RNasePRNA_drtr_bpp)


#PPV
ppv_RNasePRNA_fold_mfe = ppv(tp_RNasePRNA_fold_mfe_sum, fn_RNasePRNA_fold_mfe_sum)
# print(ppv_RNasePRNA_fold_mfe)

ppv_RNasePRNA_fold_bpp = ppv(tp_RNasePRNA_fold_bpp_sum, fn_RNasePRNA_fold_bpp_sum)
# print(ppv_RNasePRNA_fold_bpp)

ppv_RNasePRNA_drtr_mfe = ppv(tp_RNasePRNA_drtr_mfe_sum, fn_RNasePRNA_drtr_mfe_sum)
# print(ppv_RNasePRNA_drtr_mfe)

ppv_RNasePRNA_drtr_prob = ppv(tp_RNasePRNA_drtr_prob_sum, fn_RNasePRNA_drtr_prob_sum)
# print(ppv_RNasePRNA_drtr_prob)

ppv_RNasePRNA_drtr_bpp = ppv(tp_RNasePRNA_drtr_bpp_sum, fn_RNasePRNA_drtr_bpp_sum)
# print(ppv_RNasePRNA_drtr_bpp)

#F1
f_1_RNasePRNA_fold_mfe = f_1(tp_RNasePRNA_fold_mfe_sum, fp_RNasePRNA_fold_mfe_sum, fn_RNasePRNA_fold_mfe_sum)
# print(f_1_RNasePRNA_fold_mfe)

f_1_RNasePRNA_fold_bpp = f_1(tp_RNasePRNA_fold_bpp_sum, fp_RNasePRNA_fold_bpp_sum, fn_RNasePRNA_fold_bpp_sum)
# print(f_1_RNasePRNA_fold_bpp)

f_1_RNasePRNA_drtr_mfe = f_1(tp_RNasePRNA_drtr_mfe_sum, fp_RNasePRNA_drtr_mfe_sum, fn_RNasePRNA_drtr_mfe_sum)
# print(f_1_RNasePRNA_drtr_mfe)

f_1_RNasePRNA_drtr_prob = f_1(tp_RNasePRNA_drtr_prob_sum, fp_RNasePRNA_drtr_prob_sum, fn_RNasePRNA_drtr_prob_sum)
# print(f_1_RNasePRNA_drtr_prob)

f_1_RNasePRNA_drtr_bpp = f_1(tp_RNasePRNA_drtr_bpp_sum, fp_RNasePRNA_drtr_bpp_sum, fn_RNasePRNA_drtr_bpp_sum)
# print(f_1_RNasePRNA_drtr_bpp)

# #MCC
mcc_RNasePRNA_fold_mfe = mcc(tp_RNasePRNA_fold_mfe_sum, fp_RNasePRNA_fold_mfe_sum, tn_RNasePRNA_fold_mfe_sum, fn_RNasePRNA_fold_mfe_sum)
# print(mcc_RNasePRNA_fold_mfe)

mcc_RNasePRNA_fold_bpp = mcc(tp_RNasePRNA_fold_bpp_sum, fp_RNasePRNA_fold_bpp_sum, tn_RNasePRNA_fold_bpp_sum, fn_RNasePRNA_fold_bpp_sum)
# print(mcc_RNasePRNA_fold_bpp)

mcc_RNasePRNA_drtr_mfe = mcc(tp_RNasePRNA_drtr_mfe_sum, fp_RNasePRNA_drtr_mfe_sum, tn_RNasePRNA_drtr_mfe_sum, fn_RNasePRNA_drtr_mfe_sum)
# print(mcc_RNasePRNA_drtr_mfe)

mcc_RNasePRNA_drtr_prob = mcc(tp_RNasePRNA_drtr_prob_sum, fp_RNasePRNA_drtr_prob_sum, tn_RNasePRNA_drtr_prob_sum, fn_RNasePRNA_drtr_prob_sum)
# print(mcc_RNasePRNA_drtr_prob)

mcc_RNasePRNA_drtr_bpp = mcc(tp_RNasePRNA_drtr_bpp_sum, fp_RNasePRNA_drtr_bpp_sum, tn_RNasePRNA_drtr_bpp_sum, fn_RNasePRNA_drtr_bpp_sum)
# print(mcc_RNasePRNA_drtr_bpp)



# snRNA
snRNA=family_data.loc[family_data.family=='snRNA']

tp_snRNA_fold_mfe = (snRNA[['RNAfold_mfe_tp']])
tp_snRNA_fold_mfe_sum = boot(tp_snRNA_fold_mfe)

fp_snRNA_fold_mfe = (snRNA[['RNAfold_mfe_fp']])
fp_snRNA_fold_mfe_sum = boot(fp_snRNA_fold_mfe)

tn_snRNA_fold_mfe = (snRNA[['RNAfold_mfe_tn']])
tn_snRNA_fold_mfe_sum = boot(tn_snRNA_fold_mfe)

fn_snRNA_fold_mfe = (snRNA[['RNAfold_mfe_fn']])
fn_snRNA_fold_mfe_sum = boot(fn_snRNA_fold_mfe)


tp_snRNA_fold_bpp = (snRNA[['RNAfold_bpp_tp']])
tp_snRNA_fold_bpp_sum = boot(tp_snRNA_fold_bpp)

fp_snRNA_fold_bpp = (snRNA[['RNAfold_bpp_fp']])
fp_snRNA_fold_bpp_sum = boot(fp_snRNA_fold_bpp)

tn_snRNA_fold_bpp = (snRNA[['RNAfold_bpp_tn']])
tn_snRNA_fold_bpp_sum = boot(tn_snRNA_fold_bpp)

fn_snRNA_fold_bpp = (snRNA[['RNAfold_bpp_fn']])
fn_snRNA_fold_bpp_sum = boot(fn_snRNA_fold_bpp)


tp_snRNA_drtr_mfe = (snRNA[['drtr_tp_mfe']])
tp_snRNA_drtr_mfe_sum = boot(tp_snRNA_drtr_mfe)

fp_snRNA_drtr_mfe = (snRNA[['drtr_fp_mfe']])
fp_snRNA_drtr_mfe_sum = boot(fp_snRNA_drtr_mfe)

tn_snRNA_drtr_mfe = (snRNA[['drtr_tn_mfe']])
tn_snRNA_drtr_mfe_sum = boot(tn_snRNA_drtr_mfe)

fn_snRNA_drtr_mfe = (snRNA[['drtr_fn_mfe']])
fn_snRNA_drtr_mfe_sum = boot(fn_snRNA_drtr_mfe)


tp_snRNA_drtr_prob = (snRNA[['drtr_tp_prob']])
tp_snRNA_drtr_prob_sum = boot(tp_snRNA_drtr_prob)

fp_snRNA_drtr_prob = (snRNA[['drtr_fp_prob']])
fp_snRNA_drtr_prob_sum = boot(fp_snRNA_drtr_prob)

tn_snRNA_drtr_prob = (snRNA[['drtr_tn_prob']])
tn_snRNA_drtr_prob_sum = boot(tn_snRNA_drtr_prob)

fn_snRNA_drtr_prob = (snRNA[['drtr_fn_prob']])
fn_snRNA_drtr_prob_sum = boot(fn_snRNA_drtr_prob)


tp_snRNA_drtr_bpp = (snRNA[['drtr_tp']])
tp_snRNA_drtr_bpp_sum = boot(tp_snRNA_drtr_bpp)

fp_snRNA_drtr_bpp = (snRNA[['drtr_fp']])
fp_snRNA_drtr_bpp_sum = boot(fp_snRNA_drtr_bpp)

tn_snRNA_drtr_bpp = (snRNA[['drtr_tn']])
tn_snRNA_drtr_bpp_sum = boot(tn_snRNA_drtr_bpp)

fn_snRNA_drtr_bpp = (snRNA[['drtr_fn']])
fn_snRNA_drtr_bpp_sum = boot(fn_snRNA_drtr_bpp)

#TPR
tpr_snRNA_fold_mfe = tpr(tp_snRNA_fold_mfe_sum, fn_snRNA_fold_mfe_sum)
# print(tpr_snRNA_fold_mfe)

tpr_snRNA_fold_bpp = tpr(tp_snRNA_fold_bpp_sum, fn_snRNA_fold_bpp_sum)
# print(tpr_snRNA_fold_bpp)

tpr_snRNA_drtr_mfe = tpr(tp_snRNA_drtr_mfe_sum, fn_snRNA_drtr_mfe_sum)
# print(tpr_snRNA_drtr_mfe)

tpr_snRNA_drtr_prob = tpr(tp_snRNA_drtr_prob_sum, fn_snRNA_drtr_prob_sum)
# print(tpr_snRNA_drtr_prob)

tpr_snRNA_drtr_bpp = tpr(tp_snRNA_drtr_bpp_sum, fn_snRNA_drtr_bpp_sum)
# print(tpr_snRNA_drtr_bpp)


#PPV
ppv_snRNA_fold_mfe = ppv(tp_snRNA_fold_mfe_sum, fn_snRNA_fold_mfe_sum)
# print(ppv_snRNA_fold_mfe)

ppv_snRNA_fold_bpp = ppv(tp_snRNA_fold_bpp_sum, fn_snRNA_fold_bpp_sum)
# print(ppv_snRNA_fold_bpp)

ppv_snRNA_drtr_mfe = ppv(tp_snRNA_drtr_mfe_sum, fn_snRNA_drtr_mfe_sum)
# print(ppv_snRNA_drtr_mfe)

ppv_snRNA_drtr_prob = ppv(tp_snRNA_drtr_prob_sum, fn_snRNA_drtr_prob_sum)
# print(ppv_snRNA_drtr_prob)

ppv_snRNA_drtr_bpp = ppv(tp_snRNA_drtr_bpp_sum, fn_snRNA_drtr_bpp_sum)
# print(ppv_snRNA_drtr_bpp)

#F1
f_1_snRNA_fold_mfe = f_1(tp_snRNA_fold_mfe_sum, fp_snRNA_fold_mfe_sum, fn_snRNA_fold_mfe_sum)
# print(f_1_snRNA_fold_mfe)

f_1_snRNA_fold_bpp = f_1(tp_snRNA_fold_bpp_sum, fp_snRNA_fold_bpp_sum, fn_snRNA_fold_bpp_sum)
# print(f_1_snRNA_fold_bpp)

f_1_snRNA_drtr_mfe = f_1(tp_snRNA_drtr_mfe_sum, fp_snRNA_drtr_mfe_sum, fn_snRNA_drtr_mfe_sum)
# print(f_1_snRNA_drtr_mfe)

f_1_snRNA_drtr_prob = f_1(tp_snRNA_drtr_prob_sum, fp_snRNA_drtr_prob_sum, fn_snRNA_drtr_prob_sum)
# print(f_1_snRNA_drtr_prob)

f_1_snRNA_drtr_bpp = f_1(tp_snRNA_drtr_bpp_sum, fp_snRNA_drtr_bpp_sum, fn_snRNA_drtr_bpp_sum)
# print(f_1_snRNA_drtr_bpp)

# #MCC
mcc_snRNA_fold_mfe = mcc(tp_snRNA_fold_mfe_sum, fp_snRNA_fold_mfe_sum, tn_snRNA_fold_mfe_sum, fn_snRNA_fold_mfe_sum)
# print(mcc_snRNA_fold_mfe)

mcc_snRNA_fold_bpp = mcc(tp_snRNA_fold_bpp_sum, fp_snRNA_fold_bpp_sum, tn_snRNA_fold_bpp_sum, fn_snRNA_fold_bpp_sum)
# print(mcc_snRNA_fold_bpp)

mcc_snRNA_drtr_mfe = mcc(tp_snRNA_drtr_mfe_sum, fp_snRNA_drtr_mfe_sum, tn_snRNA_drtr_mfe_sum, fn_snRNA_drtr_mfe_sum)
# print(mcc_snRNA_drtr_mfe)

mcc_snRNA_drtr_prob = mcc(tp_snRNA_drtr_prob_sum, fp_snRNA_drtr_prob_sum, tn_snRNA_drtr_prob_sum, fn_snRNA_drtr_prob_sum)
# print(mcc_snRNA_drtr_prob)

mcc_snRNA_drtr_bpp = mcc(tp_snRNA_drtr_bpp_sum, fp_snRNA_drtr_bpp_sum, tn_snRNA_drtr_bpp_sum, fn_snRNA_drtr_bpp_sum)
# print(mcc_snRNA_drtr_bpp)


#SRPRNA
SRPRNA=family_data.loc[family_data.family=='SRPRNA']

tp_SRPRNA_fold_mfe = (SRPRNA[['RNAfold_mfe_tp']])
tp_SRPRNA_fold_mfe_sum = boot(tp_SRPRNA_fold_mfe)

fp_SRPRNA_fold_mfe = (SRPRNA[['RNAfold_mfe_fp']])
fp_SRPRNA_fold_mfe_sum = boot(fp_SRPRNA_fold_mfe)

tn_SRPRNA_fold_mfe = (SRPRNA[['RNAfold_mfe_tn']])
tn_SRPRNA_fold_mfe_sum = boot(tn_SRPRNA_fold_mfe)

fn_SRPRNA_fold_mfe = (SRPRNA[['RNAfold_mfe_fn']])
fn_SRPRNA_fold_mfe_sum = boot(fn_SRPRNA_fold_mfe)


tp_SRPRNA_fold_bpp = (SRPRNA[['RNAfold_bpp_tp']])
tp_SRPRNA_fold_bpp_sum = boot(tp_SRPRNA_fold_bpp)

fp_SRPRNA_fold_bpp = (SRPRNA[['RNAfold_bpp_fp']])
fp_SRPRNA_fold_bpp_sum = boot(fp_SRPRNA_fold_bpp)

tn_SRPRNA_fold_bpp = (SRPRNA[['RNAfold_bpp_tn']])
tn_SRPRNA_fold_bpp_sum = boot(tn_SRPRNA_fold_bpp)

fn_SRPRNA_fold_bpp = (SRPRNA[['RNAfold_bpp_fn']])
fn_SRPRNA_fold_bpp_sum = boot(fn_SRPRNA_fold_bpp)


tp_SRPRNA_drtr_mfe = (SRPRNA[['drtr_tp_mfe']])
tp_SRPRNA_drtr_mfe_sum = boot(tp_SRPRNA_drtr_mfe)

fp_SRPRNA_drtr_mfe = (SRPRNA[['drtr_fp_mfe']])
fp_SRPRNA_drtr_mfe_sum = boot(fp_SRPRNA_drtr_mfe)

tn_SRPRNA_drtr_mfe = (SRPRNA[['drtr_tn_mfe']])
tn_SRPRNA_drtr_mfe_sum = boot(tn_SRPRNA_drtr_mfe)

fn_SRPRNA_drtr_mfe = (SRPRNA[['drtr_fn_mfe']])
fn_SRPRNA_drtr_mfe_sum = boot(fn_SRPRNA_drtr_mfe)


tp_SRPRNA_drtr_prob = (SRPRNA[['drtr_tp_prob']])
tp_SRPRNA_drtr_prob_sum = boot(tp_SRPRNA_drtr_prob)

fp_SRPRNA_drtr_prob = (SRPRNA[['drtr_fp_prob']])
fp_SRPRNA_drtr_prob_sum = boot(fp_SRPRNA_drtr_prob)

tn_SRPRNA_drtr_prob = (SRPRNA[['drtr_tn_prob']])
tn_SRPRNA_drtr_prob_sum = boot(tn_SRPRNA_drtr_prob)

fn_SRPRNA_drtr_prob = (SRPRNA[['drtr_fn_prob']])
fn_SRPRNA_drtr_prob_sum = boot(fn_SRPRNA_drtr_prob)


tp_SRPRNA_drtr_bpp = (SRPRNA[['drtr_tp']])
tp_SRPRNA_drtr_bpp_sum = boot(tp_SRPRNA_drtr_bpp)

fp_SRPRNA_drtr_bpp = (SRPRNA[['drtr_fp']])
fp_SRPRNA_drtr_bpp_sum = boot(fp_SRPRNA_drtr_bpp)

tn_SRPRNA_drtr_bpp = (SRPRNA[['drtr_tn']])
tn_SRPRNA_drtr_bpp_sum = boot(tn_SRPRNA_drtr_bpp)

fn_SRPRNA_drtr_bpp = (SRPRNA[['drtr_fn']])
fn_SRPRNA_drtr_bpp_sum = boot(fn_SRPRNA_drtr_bpp)

#TPR
tpr_SRPRNA_fold_mfe = tpr(tp_SRPRNA_fold_mfe_sum, fn_SRPRNA_fold_mfe_sum)
# print(tpr_SRPRNA_fold_mfe)

tpr_SRPRNA_fold_bpp = tpr(tp_SRPRNA_fold_bpp_sum, fn_SRPRNA_fold_bpp_sum)
# print(tpr_SRPRNA_fold_bpp)

tpr_SRPRNA_drtr_mfe = tpr(tp_SRPRNA_drtr_mfe_sum, fn_SRPRNA_drtr_mfe_sum)
# print(tpr_SRPRNA_drtr_mfe)

tpr_SRPRNA_drtr_prob = tpr(tp_SRPRNA_drtr_prob_sum, fn_SRPRNA_drtr_prob_sum)
# print(tpr_SRPRNA_drtr_prob)

tpr_SRPRNA_drtr_bpp = tpr(tp_SRPRNA_drtr_bpp_sum, fn_SRPRNA_drtr_bpp_sum)
# print(tpr_SRPRNA_drtr_bpp)


#PPV
ppv_SRPRNA_fold_mfe = ppv(tp_SRPRNA_fold_mfe_sum, fn_SRPRNA_fold_mfe_sum)
# print(ppv_SRPRNA_fold_mfe)

ppv_SRPRNA_fold_bpp = ppv(tp_SRPRNA_fold_bpp_sum, fn_SRPRNA_fold_bpp_sum)
# print(ppv_SRPRNA_fold_bpp)

ppv_SRPRNA_drtr_mfe = ppv(tp_SRPRNA_drtr_mfe_sum, fn_SRPRNA_drtr_mfe_sum)
# print(ppv_SRPRNA_drtr_mfe)

ppv_SRPRNA_drtr_prob = ppv(tp_SRPRNA_drtr_prob_sum, fn_SRPRNA_drtr_prob_sum)
# print(ppv_SRPRNA_drtr_prob)

ppv_SRPRNA_drtr_bpp = ppv(tp_SRPRNA_drtr_bpp_sum, fn_SRPRNA_drtr_bpp_sum)
# print(ppv_SRPRNA_drtr_bpp)

#F1
f_1_SRPRNA_fold_mfe = f_1(tp_SRPRNA_fold_mfe_sum, fp_SRPRNA_fold_mfe_sum, fn_SRPRNA_fold_mfe_sum)
# print(f_1_SRPRNA_fold_mfe)

f_1_SRPRNA_fold_bpp = f_1(tp_SRPRNA_fold_bpp_sum, fp_SRPRNA_fold_bpp_sum, fn_SRPRNA_fold_bpp_sum)
# print(f_1_SRPRNA_fold_bpp)

f_1_SRPRNA_drtr_mfe = f_1(tp_SRPRNA_drtr_mfe_sum, fp_SRPRNA_drtr_mfe_sum, fn_SRPRNA_drtr_mfe_sum)
# print(f_1_SRPRNA_drtr_mfe)

f_1_SRPRNA_drtr_prob = f_1(tp_SRPRNA_drtr_prob_sum, fp_SRPRNA_drtr_prob_sum, fn_SRPRNA_drtr_prob_sum)
# print(f_1_SRPRNA_drtr_prob)

f_1_SRPRNA_drtr_bpp = f_1(tp_SRPRNA_drtr_bpp_sum, fp_SRPRNA_drtr_bpp_sum, fn_SRPRNA_drtr_bpp_sum)
# print(f_1_SRPRNA_drtr_bpp)

# #MCC
mcc_SRPRNA_fold_mfe = mcc(tp_SRPRNA_fold_mfe_sum, fp_SRPRNA_fold_mfe_sum, tn_SRPRNA_fold_mfe_sum, fn_SRPRNA_fold_mfe_sum)
# print(mcc_SRPRNA_fold_mfe)

mcc_SRPRNA_fold_bpp = mcc(tp_SRPRNA_fold_bpp_sum, fp_SRPRNA_fold_bpp_sum, tn_SRPRNA_fold_bpp_sum, fn_SRPRNA_fold_bpp_sum)
# print(mcc_SRPRNA_fold_bpp)

mcc_SRPRNA_drtr_mfe = mcc(tp_SRPRNA_drtr_mfe_sum, fp_SRPRNA_drtr_mfe_sum, tn_SRPRNA_drtr_mfe_sum, fn_SRPRNA_drtr_mfe_sum)
# print(mcc_SRPRNA_drtr_mfe)

mcc_SRPRNA_drtr_prob = mcc(tp_SRPRNA_drtr_prob_sum, fp_SRPRNA_drtr_prob_sum, tn_SRPRNA_drtr_prob_sum, fn_SRPRNA_drtr_prob_sum)
# print(mcc_SRPRNA_drtr_prob)

mcc_SRPRNA_drtr_bpp = mcc(tp_SRPRNA_drtr_bpp_sum, fp_SRPRNA_drtr_bpp_sum, tn_SRPRNA_drtr_bpp_sum, fn_SRPRNA_drtr_bpp_sum)
# print(mcc_SRPRNA_drtr_bpp)


#SyntheticRNA
SyntheticRNA=family_data.loc[family_data.family=='SyntheticRNA']

tp_SyntheticRNA_fold_mfe = (SyntheticRNA[['RNAfold_mfe_tp']])
tp_SyntheticRNA_fold_mfe_sum = boot(tp_SyntheticRNA_fold_mfe)

fp_SyntheticRNA_fold_mfe = (SyntheticRNA[['RNAfold_mfe_fp']])
fp_SyntheticRNA_fold_mfe_sum = boot(fp_SyntheticRNA_fold_mfe)

tn_SyntheticRNA_fold_mfe = (SyntheticRNA[['RNAfold_mfe_tn']])
tn_SyntheticRNA_fold_mfe_sum = boot(tn_SyntheticRNA_fold_mfe)

fn_SyntheticRNA_fold_mfe = (SyntheticRNA[['RNAfold_mfe_fn']])
fn_SyntheticRNA_fold_mfe_sum = boot(fn_SyntheticRNA_fold_mfe)


tp_SyntheticRNA_fold_bpp = (SyntheticRNA[['RNAfold_bpp_tp']])
tp_SyntheticRNA_fold_bpp_sum = boot(tp_SyntheticRNA_fold_bpp)

fp_SyntheticRNA_fold_bpp = (SyntheticRNA[['RNAfold_bpp_fp']])
fp_SyntheticRNA_fold_bpp_sum = boot(fp_SyntheticRNA_fold_bpp)

tn_SyntheticRNA_fold_bpp = (SyntheticRNA[['RNAfold_bpp_tn']])
tn_SyntheticRNA_fold_bpp_sum = boot(tn_SyntheticRNA_fold_bpp)

fn_SyntheticRNA_fold_bpp = (SyntheticRNA[['RNAfold_bpp_fn']])
fn_SyntheticRNA_fold_bpp_sum = boot(fn_SyntheticRNA_fold_bpp)


tp_SyntheticRNA_drtr_mfe = (SyntheticRNA[['drtr_tp_mfe']])
tp_SyntheticRNA_drtr_mfe_sum = boot(tp_SyntheticRNA_drtr_mfe)

fp_SyntheticRNA_drtr_mfe = (SyntheticRNA[['drtr_fp_mfe']])
fp_SyntheticRNA_drtr_mfe_sum = boot(fp_SyntheticRNA_drtr_mfe)

tn_SyntheticRNA_drtr_mfe = (SyntheticRNA[['drtr_tn_mfe']])
tn_SyntheticRNA_drtr_mfe_sum = boot(tn_SyntheticRNA_drtr_mfe)

fn_SyntheticRNA_drtr_mfe = (SyntheticRNA[['drtr_fn_mfe']])
fn_SyntheticRNA_drtr_mfe_sum = boot(fn_SyntheticRNA_drtr_mfe)


tp_SyntheticRNA_drtr_prob = (SyntheticRNA[['drtr_tp_prob']])
tp_SyntheticRNA_drtr_prob_sum = boot(tp_SyntheticRNA_drtr_prob)

fp_SyntheticRNA_drtr_prob = (SyntheticRNA[['drtr_fp_prob']])
fp_SyntheticRNA_drtr_prob_sum = boot(fp_SyntheticRNA_drtr_prob)

tn_SyntheticRNA_drtr_prob = (SyntheticRNA[['drtr_tn_prob']])
tn_SyntheticRNA_drtr_prob_sum = boot(tn_SyntheticRNA_drtr_prob)

fn_SyntheticRNA_drtr_prob = (SyntheticRNA[['drtr_fn_prob']])
fn_SyntheticRNA_drtr_prob_sum = boot(fn_SyntheticRNA_drtr_prob)


tp_SyntheticRNA_drtr_bpp = (SyntheticRNA[['drtr_tp']])
tp_SyntheticRNA_drtr_bpp_sum = boot(tp_SyntheticRNA_drtr_bpp)

fp_SyntheticRNA_drtr_bpp = (SyntheticRNA[['drtr_fp']])
fp_SyntheticRNA_drtr_bpp_sum = boot(fp_SyntheticRNA_drtr_bpp)

tn_SyntheticRNA_drtr_bpp = (SyntheticRNA[['drtr_tn']])
tn_SyntheticRNA_drtr_bpp_sum = boot(tn_SyntheticRNA_drtr_bpp)

fn_SyntheticRNA_drtr_bpp = (SyntheticRNA[['drtr_fn']])
fn_SyntheticRNA_drtr_bpp_sum = boot(fn_SyntheticRNA_drtr_bpp)

#TPR
tpr_SyntheticRNA_fold_mfe = tpr(tp_SyntheticRNA_fold_mfe_sum, fn_SyntheticRNA_fold_mfe_sum)
# print(tpr_SyntheticRNA_fold_mfe)

tpr_SyntheticRNA_fold_bpp = tpr(tp_SyntheticRNA_fold_bpp_sum, fn_SyntheticRNA_fold_bpp_sum)
# print(tpr_SyntheticRNA_fold_bpp)

tpr_SyntheticRNA_drtr_mfe = tpr(tp_SyntheticRNA_drtr_mfe_sum, fn_SyntheticRNA_drtr_mfe_sum)
# print(tpr_SyntheticRNA_drtr_mfe)

tpr_SyntheticRNA_drtr_prob = tpr(tp_SyntheticRNA_drtr_prob_sum, fn_SyntheticRNA_drtr_prob_sum)
# print(tpr_SyntheticRNA_drtr_prob)

tpr_SyntheticRNA_drtr_bpp = tpr(tp_SyntheticRNA_drtr_bpp_sum, fn_SyntheticRNA_drtr_bpp_sum)
# print(tpr_SyntheticRNA_drtr_bpp)


#PPV
ppv_SyntheticRNA_fold_mfe = ppv(tp_SyntheticRNA_fold_mfe_sum, fn_SyntheticRNA_fold_mfe_sum)
# print(ppv_SyntheticRNA_fold_mfe)

ppv_SyntheticRNA_fold_bpp = ppv(tp_SyntheticRNA_fold_bpp_sum, fn_SyntheticRNA_fold_bpp_sum)
# print(ppv_SyntheticRNA_fold_bpp)

ppv_SyntheticRNA_drtr_mfe = ppv(tp_SyntheticRNA_drtr_mfe_sum, fn_SyntheticRNA_drtr_mfe_sum)
# print(ppv_SyntheticRNA_drtr_mfe)

ppv_SyntheticRNA_drtr_prob = ppv(tp_SyntheticRNA_drtr_prob_sum, fn_SyntheticRNA_drtr_prob_sum)
# print(ppv_SyntheticRNA_drtr_prob)

ppv_SyntheticRNA_drtr_bpp = ppv(tp_SyntheticRNA_drtr_bpp_sum, fn_SyntheticRNA_drtr_bpp_sum)
# print(ppv_SyntheticRNA_drtr_bpp)

#F1
f_1_SyntheticRNA_fold_mfe = f_1(tp_SyntheticRNA_fold_mfe_sum, fp_SyntheticRNA_fold_mfe_sum, fn_SyntheticRNA_fold_mfe_sum)
# print(f_1_SyntheticRNA_fold_mfe)

f_1_SyntheticRNA_fold_bpp = f_1(tp_SyntheticRNA_fold_bpp_sum, fp_SyntheticRNA_fold_bpp_sum, fn_SyntheticRNA_fold_bpp_sum)
# print(f_1_SyntheticRNA_fold_bpp)

f_1_SyntheticRNA_drtr_mfe = f_1(tp_SyntheticRNA_drtr_mfe_sum, fp_SyntheticRNA_drtr_mfe_sum, fn_SyntheticRNA_drtr_mfe_sum)
# print(f_1_SyntheticRNA_drtr_mfe)

f_1_SyntheticRNA_drtr_prob = f_1(tp_SyntheticRNA_drtr_prob_sum, fp_SyntheticRNA_drtr_prob_sum, fn_SyntheticRNA_drtr_prob_sum)
# print(f_1_SyntheticRNA_drtr_prob)

f_1_SyntheticRNA_drtr_bpp = f_1(tp_SyntheticRNA_drtr_bpp_sum, fp_SyntheticRNA_drtr_bpp_sum, fn_SyntheticRNA_drtr_bpp_sum)
# print(f_1_SyntheticRNA_drtr_bpp)

# #MCC
mcc_SyntheticRNA_fold_mfe = mcc(tp_SyntheticRNA_fold_mfe_sum, fp_SyntheticRNA_fold_mfe_sum, tn_SyntheticRNA_fold_mfe_sum, fn_SyntheticRNA_fold_mfe_sum)
# print(mcc_SyntheticRNA_fold_mfe)

mcc_SyntheticRNA_fold_bpp = mcc(tp_SyntheticRNA_fold_bpp_sum, fp_SyntheticRNA_fold_bpp_sum, tn_SyntheticRNA_fold_bpp_sum, fn_SyntheticRNA_fold_bpp_sum)
# print(mcc_SyntheticRNA_fold_bpp)

mcc_SyntheticRNA_drtr_mfe = mcc(tp_SyntheticRNA_drtr_mfe_sum, fp_SyntheticRNA_drtr_mfe_sum, tn_SyntheticRNA_drtr_mfe_sum, fn_SyntheticRNA_drtr_mfe_sum)
# print(mcc_SyntheticRNA_drtr_mfe)

mcc_SyntheticRNA_drtr_prob = mcc(tp_SyntheticRNA_drtr_prob_sum, fp_SyntheticRNA_drtr_prob_sum, tn_SyntheticRNA_drtr_prob_sum, fn_SyntheticRNA_drtr_prob_sum)
# print(mcc_SyntheticRNA_drtr_prob)

mcc_SyntheticRNA_drtr_bpp = mcc(tp_SyntheticRNA_drtr_bpp_sum, fp_SyntheticRNA_drtr_bpp_sum, tn_SyntheticRNA_drtr_bpp_sum, fn_SyntheticRNA_drtr_bpp_sum)
# print(mcc_SyntheticRNA_drtr_bpp)


#tmRNA
tmRNA=family_data.loc[family_data.family=='tmRNA']

tp_tmRNA_fold_mfe = (tmRNA[['RNAfold_mfe_tp']])
tp_tmRNA_fold_mfe_sum = boot(tp_tmRNA_fold_mfe)

fp_tmRNA_fold_mfe = (tmRNA[['RNAfold_mfe_fp']])
fp_tmRNA_fold_mfe_sum = boot(fp_tmRNA_fold_mfe)

tn_tmRNA_fold_mfe = (tmRNA[['RNAfold_mfe_tn']])
tn_tmRNA_fold_mfe_sum = boot(tn_tmRNA_fold_mfe)

fn_tmRNA_fold_mfe = (tmRNA[['RNAfold_mfe_fn']])
fn_tmRNA_fold_mfe_sum = boot(fn_tmRNA_fold_mfe)


tp_tmRNA_fold_bpp = (tmRNA[['RNAfold_bpp_tp']])
tp_tmRNA_fold_bpp_sum = boot(tp_tmRNA_fold_bpp)

fp_tmRNA_fold_bpp = (tmRNA[['RNAfold_bpp_fp']])
fp_tmRNA_fold_bpp_sum = boot(fp_tmRNA_fold_bpp)

tn_tmRNA_fold_bpp = (tmRNA[['RNAfold_bpp_tn']])
tn_tmRNA_fold_bpp_sum = boot(tn_tmRNA_fold_bpp)

fn_tmRNA_fold_bpp = (tmRNA[['RNAfold_bpp_fn']])
fn_tmRNA_fold_bpp_sum = boot(fn_tmRNA_fold_bpp)


tp_tmRNA_drtr_mfe = (tmRNA[['drtr_tp_mfe']])
tp_tmRNA_drtr_mfe_sum = boot(tp_tmRNA_drtr_mfe)

fp_tmRNA_drtr_mfe = (tmRNA[['drtr_fp_mfe']])
fp_tmRNA_drtr_mfe_sum = boot(fp_tmRNA_drtr_mfe)

tn_tmRNA_drtr_mfe = (tmRNA[['drtr_tn_mfe']])
tn_tmRNA_drtr_mfe_sum = boot(tn_tmRNA_drtr_mfe)

fn_tmRNA_drtr_mfe = (tmRNA[['drtr_fn_mfe']])
fn_tmRNA_drtr_mfe_sum = boot(fn_tmRNA_drtr_mfe)


tp_tmRNA_drtr_prob = (tmRNA[['drtr_tp_prob']])
tp_tmRNA_drtr_prob_sum = boot(tp_tmRNA_drtr_prob)

fp_tmRNA_drtr_prob = (tmRNA[['drtr_fp_prob']])
fp_tmRNA_drtr_prob_sum = boot(fp_tmRNA_drtr_prob)

tn_tmRNA_drtr_prob = (tmRNA[['drtr_tn_prob']])
tn_tmRNA_drtr_prob_sum = boot(tn_tmRNA_drtr_prob)

fn_tmRNA_drtr_prob = (tmRNA[['drtr_fn_prob']])
fn_tmRNA_drtr_prob_sum = boot(fn_tmRNA_drtr_prob)


tp_tmRNA_drtr_bpp = (tmRNA[['drtr_tp']])
tp_tmRNA_drtr_bpp_sum = boot(tp_tmRNA_drtr_bpp)

fp_tmRNA_drtr_bpp = (tmRNA[['drtr_fp']])
fp_tmRNA_drtr_bpp_sum = boot(fp_tmRNA_drtr_bpp)

tn_tmRNA_drtr_bpp = (tmRNA[['drtr_tn']])
tn_tmRNA_drtr_bpp_sum = boot(tn_tmRNA_drtr_bpp)

fn_tmRNA_drtr_bpp = (tmRNA[['drtr_fn']])
fn_tmRNA_drtr_bpp_sum = boot(fn_tmRNA_drtr_bpp)

#TPR
tpr_tmRNA_fold_mfe = tpr(tp_tmRNA_fold_mfe_sum, fn_tmRNA_fold_mfe_sum)
# print(tpr_tmRNA_fold_mfe)

tpr_tmRNA_fold_bpp = tpr(tp_tmRNA_fold_bpp_sum, fn_tmRNA_fold_bpp_sum)
# print(tpr_tmRNA_fold_bpp)

tpr_tmRNA_drtr_mfe = tpr(tp_tmRNA_drtr_mfe_sum, fn_tmRNA_drtr_mfe_sum)
# print(tpr_tmRNA_drtr_mfe)

tpr_tmRNA_drtr_prob = tpr(tp_tmRNA_drtr_prob_sum, fn_tmRNA_drtr_prob_sum)
# print(tpr_tmRNA_drtr_prob)

tpr_tmRNA_drtr_bpp = tpr(tp_tmRNA_drtr_bpp_sum, fn_tmRNA_drtr_bpp_sum)
# print(tpr_tmRNA_drtr_bpp)


#PPV
ppv_tmRNA_fold_mfe = ppv(tp_tmRNA_fold_mfe_sum, fn_tmRNA_fold_mfe_sum)
# print(ppv_tmRNA_fold_mfe)

ppv_tmRNA_fold_bpp = ppv(tp_tmRNA_fold_bpp_sum, fn_tmRNA_fold_bpp_sum)
# print(ppv_tmRNA_fold_bpp)

ppv_tmRNA_drtr_mfe = ppv(tp_tmRNA_drtr_mfe_sum, fn_tmRNA_drtr_mfe_sum)
# print(ppv_tmRNA_drtr_mfe)

ppv_tmRNA_drtr_prob = ppv(tp_tmRNA_drtr_prob_sum, fn_tmRNA_drtr_prob_sum)
# print(ppv_tmRNA_drtr_prob)

ppv_tmRNA_drtr_bpp = ppv(tp_tmRNA_drtr_bpp_sum, fn_tmRNA_drtr_bpp_sum)
# print(ppv_tmRNA_drtr_bpp)

#F1
f_1_tmRNA_fold_mfe = f_1(tp_tmRNA_fold_mfe_sum, fp_tmRNA_fold_mfe_sum, fn_tmRNA_fold_mfe_sum)
# print(f_1_tmRNA_fold_mfe)

f_1_tmRNA_fold_bpp = f_1(tp_tmRNA_fold_bpp_sum, fp_tmRNA_fold_bpp_sum, fn_tmRNA_fold_bpp_sum)
# print(f_1_tmRNA_fold_bpp)

f_1_tmRNA_drtr_mfe = f_1(tp_tmRNA_drtr_mfe_sum, fp_tmRNA_drtr_mfe_sum, fn_tmRNA_drtr_mfe_sum)
# print(f_1_tmRNA_drtr_mfe)

f_1_tmRNA_drtr_prob = f_1(tp_tmRNA_drtr_prob_sum, fp_tmRNA_drtr_prob_sum, fn_tmRNA_drtr_prob_sum)
# print(f_1_tmRNA_drtr_prob)

f_1_tmRNA_drtr_bpp = f_1(tp_tmRNA_drtr_bpp_sum, fp_tmRNA_drtr_bpp_sum, fn_tmRNA_drtr_bpp_sum)
# print(f_1_tmRNA_drtr_bpp)

# #MCC
mcc_tmRNA_fold_mfe = mcc(tp_tmRNA_fold_mfe_sum, fp_tmRNA_fold_mfe_sum, tn_tmRNA_fold_mfe_sum, fn_tmRNA_fold_mfe_sum)
# print(mcc_tmRNA_fold_mfe)

mcc_tmRNA_fold_bpp = mcc(tp_tmRNA_fold_bpp_sum, fp_tmRNA_fold_bpp_sum, tn_tmRNA_fold_bpp_sum, fn_tmRNA_fold_bpp_sum)
# print(mcc_tmRNA_fold_bpp)

mcc_tmRNA_drtr_mfe = mcc(tp_tmRNA_drtr_mfe_sum, fp_tmRNA_drtr_mfe_sum, tn_tmRNA_drtr_mfe_sum, fn_tmRNA_drtr_mfe_sum)
# print(mcc_tmRNA_drtr_mfe)

mcc_tmRNA_drtr_prob = mcc(tp_tmRNA_drtr_prob_sum, fp_tmRNA_drtr_prob_sum, tn_tmRNA_drtr_prob_sum, fn_tmRNA_drtr_prob_sum)
# print(mcc_tmRNA_drtr_prob)

mcc_tmRNA_drtr_bpp = mcc(tp_tmRNA_drtr_bpp_sum, fp_tmRNA_drtr_bpp_sum, tn_tmRNA_drtr_bpp_sum, fn_tmRNA_drtr_bpp_sum)
# print(mcc_tmRNA_drtr_bpp)


#tRNA
tRNA=family_data.loc[family_data.family=='tRNA']

tp_tRNA_fold_mfe = (tRNA[['RNAfold_mfe_tp']])
tp_tRNA_fold_mfe_sum = boot(tp_tRNA_fold_mfe)

fp_tRNA_fold_mfe = (tRNA[['RNAfold_mfe_fp']])
fp_tRNA_fold_mfe_sum = boot(fp_tRNA_fold_mfe)

tn_tRNA_fold_mfe = (tRNA[['RNAfold_mfe_tn']])
tn_tRNA_fold_mfe_sum = boot(tn_tRNA_fold_mfe)

fn_tRNA_fold_mfe = (tRNA[['RNAfold_mfe_fn']])
fn_tRNA_fold_mfe_sum = boot(fn_tRNA_fold_mfe)


tp_tRNA_fold_bpp = (tRNA[['RNAfold_bpp_tp']])
tp_tRNA_fold_bpp_sum = boot(tp_tRNA_fold_bpp)

fp_tRNA_fold_bpp = (tRNA[['RNAfold_bpp_fp']])
fp_tRNA_fold_bpp_sum = boot(fp_tRNA_fold_bpp)

tn_tRNA_fold_bpp = (tRNA[['RNAfold_bpp_tn']])
tn_tRNA_fold_bpp_sum = boot(tn_tRNA_fold_bpp)

fn_tRNA_fold_bpp = (tRNA[['RNAfold_bpp_fn']])
fn_tRNA_fold_bpp_sum = boot(fn_tRNA_fold_bpp)


tp_tRNA_drtr_mfe = (tRNA[['drtr_tp_mfe']])
tp_tRNA_drtr_mfe_sum = boot(tp_tRNA_drtr_mfe)

fp_tRNA_drtr_mfe = (tRNA[['drtr_fp_mfe']])
fp_tRNA_drtr_mfe_sum = boot(fp_tRNA_drtr_mfe)

tn_tRNA_drtr_mfe = (tRNA[['drtr_tn_mfe']])
tn_tRNA_drtr_mfe_sum = boot(tn_tRNA_drtr_mfe)

fn_tRNA_drtr_mfe = (tRNA[['drtr_fn_mfe']])
fn_tRNA_drtr_mfe_sum = boot(fn_tRNA_drtr_mfe)


tp_tRNA_drtr_prob = (tRNA[['drtr_tp_prob']])
tp_tRNA_drtr_prob_sum = boot(tp_tRNA_drtr_prob)

fp_tRNA_drtr_prob = (tRNA[['drtr_fp_prob']])
fp_tRNA_drtr_prob_sum = boot(fp_tRNA_drtr_prob)

tn_tRNA_drtr_prob = (tRNA[['drtr_tn_prob']])
tn_tRNA_drtr_prob_sum = boot(tn_tRNA_drtr_prob)

fn_tRNA_drtr_prob = (tRNA[['drtr_fn_prob']])
fn_tRNA_drtr_prob_sum = boot(fn_tRNA_drtr_prob)


tp_tRNA_drtr_bpp = (tRNA[['drtr_tp']])
tp_tRNA_drtr_bpp_sum = boot(tp_tRNA_drtr_bpp)

fp_tRNA_drtr_bpp = (tRNA[['drtr_fp']])
fp_tRNA_drtr_bpp_sum = boot(fp_tRNA_drtr_bpp)

tn_tRNA_drtr_bpp = (tRNA[['drtr_tn']])
tn_tRNA_drtr_bpp_sum = boot(tn_tRNA_drtr_bpp)

fn_tRNA_drtr_bpp = (tRNA[['drtr_fn']])
fn_tRNA_drtr_bpp_sum = boot(fn_tRNA_drtr_bpp)

#TPR
tpr_tRNA_fold_mfe = tpr(tp_tRNA_fold_mfe_sum, fn_tRNA_fold_mfe_sum)
# print(tpr_tRNA_fold_mfe)

tpr_tRNA_fold_bpp = tpr(tp_tRNA_fold_bpp_sum, fn_tRNA_fold_bpp_sum)
# print(tpr_tRNA_fold_bpp)

tpr_tRNA_drtr_mfe = tpr(tp_tRNA_drtr_mfe_sum, fn_tRNA_drtr_mfe_sum)
# print(tpr_tRNA_drtr_mfe)

tpr_tRNA_drtr_prob = tpr(tp_tRNA_drtr_prob_sum, fn_tRNA_drtr_prob_sum)
# print(tpr_tRNA_drtr_prob)

tpr_tRNA_drtr_bpp = tpr(tp_tRNA_drtr_bpp_sum, fn_tRNA_drtr_bpp_sum)
# print(tpr_tRNA_drtr_bpp)


#PPV
ppv_tRNA_fold_mfe = ppv(tp_tRNA_fold_mfe_sum, fn_tRNA_fold_mfe_sum)
# print(ppv_tRNA_fold_mfe)

ppv_tRNA_fold_bpp = ppv(tp_tRNA_fold_bpp_sum, fn_tRNA_fold_bpp_sum)
# print(ppv_tRNA_fold_bpp)

ppv_tRNA_drtr_mfe = ppv(tp_tRNA_drtr_mfe_sum, fn_tRNA_drtr_mfe_sum)
# print(ppv_tRNA_drtr_mfe)

ppv_tRNA_drtr_prob = ppv(tp_tRNA_drtr_prob_sum, fn_tRNA_drtr_prob_sum)
# print(ppv_tRNA_drtr_prob)

ppv_tRNA_drtr_bpp = ppv(tp_tRNA_drtr_bpp_sum, fn_tRNA_drtr_bpp_sum)
# print(ppv_tRNA_drtr_bpp)

#F1
f_1_tRNA_fold_mfe = f_1(tp_tRNA_fold_mfe_sum, fp_tRNA_fold_mfe_sum, fn_tRNA_fold_mfe_sum)
# print(f_1_tRNA_fold_mfe)

f_1_tRNA_fold_bpp = f_1(tp_tRNA_fold_bpp_sum, fp_tRNA_fold_bpp_sum, fn_tRNA_fold_bpp_sum)
# print(f_1_tRNA_fold_bpp)

f_1_tRNA_drtr_mfe = f_1(tp_tRNA_drtr_mfe_sum, fp_tRNA_drtr_mfe_sum, fn_tRNA_drtr_mfe_sum)
# print(f_1_tRNA_drtr_mfe)

f_1_tRNA_drtr_prob = f_1(tp_tRNA_drtr_prob_sum, fp_tRNA_drtr_prob_sum, fn_tRNA_drtr_prob_sum)
# print(f_1_tRNA_drtr_prob)

f_1_tRNA_drtr_bpp = f_1(tp_tRNA_drtr_bpp_sum, fp_tRNA_drtr_bpp_sum, fn_tRNA_drtr_bpp_sum)
# print(f_1_tRNA_drtr_bpp)

# #MCC
mcc_tRNA_fold_mfe = mcc(tp_tRNA_fold_mfe_sum, fp_tRNA_fold_mfe_sum, tn_tRNA_fold_mfe_sum, fn_tRNA_fold_mfe_sum)
# print(mcc_tRNA_fold_mfe)

mcc_tRNA_fold_bpp = mcc(tp_tRNA_fold_bpp_sum, fp_tRNA_fold_bpp_sum, tn_tRNA_fold_bpp_sum, fn_tRNA_fold_bpp_sum)
# print(mcc_tRNA_fold_bpp)

mcc_tRNA_drtr_mfe = mcc(tp_tRNA_drtr_mfe_sum, fp_tRNA_drtr_mfe_sum, tn_tRNA_drtr_mfe_sum, fn_tRNA_drtr_mfe_sum)
# print(mcc_tRNA_drtr_mfe)

mcc_tRNA_drtr_prob = mcc(tp_tRNA_drtr_prob_sum, fp_tRNA_drtr_prob_sum, tn_tRNA_drtr_prob_sum, fn_tRNA_drtr_prob_sum)
# print(mcc_tRNA_drtr_prob)

mcc_tRNA_drtr_bpp = mcc(tp_tRNA_drtr_bpp_sum, fp_tRNA_drtr_bpp_sum, tn_tRNA_drtr_bpp_sum, fn_tRNA_drtr_bpp_sum)
# print(mcc_tRNA_drtr_bpp)


#Vert.Telo. RNA
Vert_Telo_RNA=family_data.loc[family_data.family=='Vert.Telo. RNA']

tp_Vert_Telo_RNA_fold_mfe = (Vert_Telo_RNA[['RNAfold_mfe_tp']])
tp_Vert_Telo_RNA_fold_mfe_sum = boot(tp_Vert_Telo_RNA_fold_mfe)

fp_Vert_Telo_RNA_fold_mfe = (Vert_Telo_RNA[['RNAfold_mfe_fp']])
fp_Vert_Telo_RNA_fold_mfe_sum = boot(fp_Vert_Telo_RNA_fold_mfe)

tn_Vert_Telo_RNA_fold_mfe = (Vert_Telo_RNA[['RNAfold_mfe_tn']])
tn_Vert_Telo_RNA_fold_mfe_sum = boot(tn_Vert_Telo_RNA_fold_mfe)

fn_Vert_Telo_RNA_fold_mfe = (Vert_Telo_RNA[['RNAfold_mfe_fn']])
fn_Vert_Telo_RNA_fold_mfe_sum = boot(fn_Vert_Telo_RNA_fold_mfe)


tp_Vert_Telo_RNA_fold_bpp = (Vert_Telo_RNA[['RNAfold_bpp_tp']])
tp_Vert_Telo_RNA_fold_bpp_sum = boot(tp_Vert_Telo_RNA_fold_bpp)

fp_Vert_Telo_RNA_fold_bpp = (Vert_Telo_RNA[['RNAfold_bpp_fp']])
fp_Vert_Telo_RNA_fold_bpp_sum = boot(fp_Vert_Telo_RNA_fold_bpp)

tn_Vert_Telo_RNA_fold_bpp = (Vert_Telo_RNA[['RNAfold_bpp_tn']])
tn_Vert_Telo_RNA_fold_bpp_sum = boot(tn_Vert_Telo_RNA_fold_bpp)

fn_Vert_Telo_RNA_fold_bpp = (Vert_Telo_RNA[['RNAfold_bpp_fn']])
fn_Vert_Telo_RNA_fold_bpp_sum = boot(fn_Vert_Telo_RNA_fold_bpp)


tp_Vert_Telo_RNA_drtr_mfe = (Vert_Telo_RNA[['drtr_tp_mfe']])
tp_Vert_Telo_RNA_drtr_mfe_sum = boot(tp_Vert_Telo_RNA_drtr_mfe)

fp_Vert_Telo_RNA_drtr_mfe = (Vert_Telo_RNA[['drtr_fp_mfe']])
fp_Vert_Telo_RNA_drtr_mfe_sum = boot(fp_Vert_Telo_RNA_drtr_mfe)

tn_Vert_Telo_RNA_drtr_mfe = (Vert_Telo_RNA[['drtr_tn_mfe']])
tn_Vert_Telo_RNA_drtr_mfe_sum = boot(tn_Vert_Telo_RNA_drtr_mfe)

fn_Vert_Telo_RNA_drtr_mfe = (Vert_Telo_RNA[['drtr_fn_mfe']])
fn_Vert_Telo_RNA_drtr_mfe_sum = boot(fn_Vert_Telo_RNA_drtr_mfe)


tp_Vert_Telo_RNA_drtr_prob = (Vert_Telo_RNA[['drtr_tp_prob']])
tp_Vert_Telo_RNA_drtr_prob_sum = boot(tp_Vert_Telo_RNA_drtr_prob)

fp_Vert_Telo_RNA_drtr_prob = (Vert_Telo_RNA[['drtr_fp_prob']])
fp_Vert_Telo_RNA_drtr_prob_sum = boot(fp_Vert_Telo_RNA_drtr_prob)

tn_Vert_Telo_RNA_drtr_prob = (Vert_Telo_RNA[['drtr_tn_prob']])
tn_Vert_Telo_RNA_drtr_prob_sum = boot(tn_Vert_Telo_RNA_drtr_prob)

fn_Vert_Telo_RNA_drtr_prob = (Vert_Telo_RNA[['drtr_fn_prob']])
fn_Vert_Telo_RNA_drtr_prob_sum = boot(fn_Vert_Telo_RNA_drtr_prob)


tp_Vert_Telo_RNA_drtr_bpp = (Vert_Telo_RNA[['drtr_tp']])
tp_Vert_Telo_RNA_drtr_bpp_sum = boot(tp_Vert_Telo_RNA_drtr_bpp)

fp_Vert_Telo_RNA_drtr_bpp = (Vert_Telo_RNA[['drtr_fp']])
fp_Vert_Telo_RNA_drtr_bpp_sum = boot(fp_Vert_Telo_RNA_drtr_bpp)

tn_Vert_Telo_RNA_drtr_bpp = (Vert_Telo_RNA[['drtr_tn']])
tn_Vert_Telo_RNA_drtr_bpp_sum = boot(tn_Vert_Telo_RNA_drtr_bpp)

fn_Vert_Telo_RNA_drtr_bpp = (Vert_Telo_RNA[['drtr_fn']])
fn_Vert_Telo_RNA_drtr_bpp_sum = boot(fn_Vert_Telo_RNA_drtr_bpp)

#TPR
tpr_Vert_Telo_RNA_fold_mfe = tpr(tp_Vert_Telo_RNA_fold_mfe_sum, fn_Vert_Telo_RNA_fold_mfe_sum)
# print(tpr_Vert_Telo_RNA_fold_mfe)

tpr_Vert_Telo_RNA_fold_bpp = tpr(tp_Vert_Telo_RNA_fold_bpp_sum, fn_Vert_Telo_RNA_fold_bpp_sum)
# print(tpr_Vert_Telo_RNA_fold_bpp)

tpr_Vert_Telo_RNA_drtr_mfe = tpr(tp_Vert_Telo_RNA_drtr_mfe_sum, fn_Vert_Telo_RNA_drtr_mfe_sum)
# print(tpr_Vert_Telo_RNA_drtr_mfe)

tpr_Vert_Telo_RNA_drtr_prob = tpr(tp_Vert_Telo_RNA_drtr_prob_sum, fn_Vert_Telo_RNA_drtr_prob_sum)
# print(tpr_Vert_Telo_RNA_drtr_prob)

tpr_Vert_Telo_RNA_drtr_bpp = tpr(tp_Vert_Telo_RNA_drtr_bpp_sum, fn_Vert_Telo_RNA_drtr_bpp_sum)
# print(tpr_Vert_Telo_RNA_drtr_bpp)


#PPV
ppv_Vert_Telo_RNA_fold_mfe = ppv(tp_Vert_Telo_RNA_fold_mfe_sum, fn_Vert_Telo_RNA_fold_mfe_sum)
# print(ppv_Vert_Telo_RNA_fold_mfe)

ppv_Vert_Telo_RNA_fold_bpp = ppv(tp_Vert_Telo_RNA_fold_bpp_sum, fn_Vert_Telo_RNA_fold_bpp_sum)
# print(ppv_Vert_Telo_RNA_fold_bpp)

ppv_Vert_Telo_RNA_drtr_mfe = ppv(tp_Vert_Telo_RNA_drtr_mfe_sum, fn_Vert_Telo_RNA_drtr_mfe_sum)
# print(ppv_Vert_Telo_RNA_drtr_mfe)

ppv_Vert_Telo_RNA_drtr_prob = ppv(tp_Vert_Telo_RNA_drtr_prob_sum, fn_Vert_Telo_RNA_drtr_prob_sum)
# print(ppv_Vert_Telo_RNA_drtr_prob)

ppv_Vert_Telo_RNA_drtr_bpp = ppv(tp_Vert_Telo_RNA_drtr_bpp_sum, fn_Vert_Telo_RNA_drtr_bpp_sum)
# print(ppv_Vert_Telo_RNA_drtr_bpp)

#F1
f_1_Vert_Telo_RNA_fold_mfe = f_1(tp_Vert_Telo_RNA_fold_mfe_sum, fp_Vert_Telo_RNA_fold_mfe_sum, fn_Vert_Telo_RNA_fold_mfe_sum)
# print(f_1_Vert_Telo_RNA_fold_mfe)

f_1_Vert_Telo_RNA_fold_bpp = f_1(tp_Vert_Telo_RNA_fold_bpp_sum, fp_Vert_Telo_RNA_fold_bpp_sum, fn_Vert_Telo_RNA_fold_bpp_sum)
# print(f_1_Vert_Telo_RNA_fold_bpp)

f_1_Vert_Telo_RNA_drtr_mfe = f_1(tp_Vert_Telo_RNA_drtr_mfe_sum, fp_Vert_Telo_RNA_drtr_mfe_sum, fn_Vert_Telo_RNA_drtr_mfe_sum)
# print(f_1_Vert_Telo_RNA_drtr_mfe)

f_1_Vert_Telo_RNA_drtr_prob = f_1(tp_Vert_Telo_RNA_drtr_prob_sum, fp_Vert_Telo_RNA_drtr_prob_sum, fn_Vert_Telo_RNA_drtr_prob_sum)
# print(f_1_Vert_Telo_RNA_drtr_prob)

f_1_Vert_Telo_RNA_drtr_bpp = f_1(tp_Vert_Telo_RNA_drtr_bpp_sum, fp_Vert_Telo_RNA_drtr_bpp_sum, fn_Vert_Telo_RNA_drtr_bpp_sum)
# print(f_1_Vert_Telo_RNA_drtr_bpp)

# #MCC
mcc_Vert_Telo_RNA_fold_mfe = mcc(tp_Vert_Telo_RNA_fold_mfe_sum, fp_Vert_Telo_RNA_fold_mfe_sum, tn_Vert_Telo_RNA_fold_mfe_sum, fn_Vert_Telo_RNA_fold_mfe_sum)
# print(mcc_Vert_Telo_RNA_fold_mfe)

mcc_Vert_Telo_RNA_fold_bpp = mcc(tp_Vert_Telo_RNA_fold_bpp_sum, fp_Vert_Telo_RNA_fold_bpp_sum, tn_Vert_Telo_RNA_fold_bpp_sum, fn_Vert_Telo_RNA_fold_bpp_sum)
# print(mcc_Vert_Telo_RNA_fold_bpp)

mcc_Vert_Telo_RNA_drtr_mfe = mcc(tp_Vert_Telo_RNA_drtr_mfe_sum, fp_Vert_Telo_RNA_drtr_mfe_sum, tn_Vert_Telo_RNA_drtr_mfe_sum, fn_Vert_Telo_RNA_drtr_mfe_sum)
# print(mcc_Vert_Telo_RNA_drtr_mfe)

mcc_Vert_Telo_RNA_drtr_prob = mcc(tp_Vert_Telo_RNA_drtr_prob_sum, fp_Vert_Telo_RNA_drtr_prob_sum, tn_Vert_Telo_RNA_drtr_prob_sum, fn_Vert_Telo_RNA_drtr_prob_sum)
# print(mcc_Vert_Telo_RNA_drtr_prob)

mcc_Vert_Telo_RNA_drtr_bpp = mcc(tp_Vert_Telo_RNA_drtr_bpp_sum, fp_Vert_Telo_RNA_drtr_bpp_sum, tn_Vert_Telo_RNA_drtr_bpp_sum, fn_Vert_Telo_RNA_drtr_bpp_sum)
# print(mcc_Vert_Telo_RNA_drtr_bpp)


#Viral&Phage
Viral_Phage=family_data.loc[family_data.family=='Viral&Phage']

tp_Viral_Phage_fold_mfe = (Viral_Phage[['RNAfold_mfe_tp']])
tp_Viral_Phage_fold_mfe_sum = boot(tp_Viral_Phage_fold_mfe)

fp_Viral_Phage_fold_mfe = (Viral_Phage[['RNAfold_mfe_fp']])
fp_Viral_Phage_fold_mfe_sum = boot(fp_Viral_Phage_fold_mfe)

tn_Viral_Phage_fold_mfe = (Viral_Phage[['RNAfold_mfe_tn']])
tn_Viral_Phage_fold_mfe_sum = boot(tn_Viral_Phage_fold_mfe)

fn_Viral_Phage_fold_mfe = (Viral_Phage[['RNAfold_mfe_fn']])
fn_Viral_Phage_fold_mfe_sum = boot(fn_Viral_Phage_fold_mfe)


tp_Viral_Phage_fold_bpp = (Viral_Phage[['RNAfold_bpp_tp']])
tp_Viral_Phage_fold_bpp_sum = boot(tp_Viral_Phage_fold_bpp)

fp_Viral_Phage_fold_bpp = (Viral_Phage[['RNAfold_bpp_fp']])
fp_Viral_Phage_fold_bpp_sum = boot(fp_Viral_Phage_fold_bpp)

tn_Viral_Phage_fold_bpp = (Viral_Phage[['RNAfold_bpp_tn']])
tn_Viral_Phage_fold_bpp_sum = boot(tn_Viral_Phage_fold_bpp)

fn_Viral_Phage_fold_bpp = (Viral_Phage[['RNAfold_bpp_fn']])
fn_Viral_Phage_fold_bpp_sum = boot(fn_Viral_Phage_fold_bpp)


tp_Viral_Phage_drtr_mfe = (Viral_Phage[['drtr_tp_mfe']])
tp_Viral_Phage_drtr_mfe_sum = boot(tp_Viral_Phage_drtr_mfe)

fp_Viral_Phage_drtr_mfe = (Viral_Phage[['drtr_fp_mfe']])
fp_Viral_Phage_drtr_mfe_sum = boot(fp_Viral_Phage_drtr_mfe)

tn_Viral_Phage_drtr_mfe = (Viral_Phage[['drtr_tn_mfe']])
tn_Viral_Phage_drtr_mfe_sum = boot(tn_Viral_Phage_drtr_mfe)

fn_Viral_Phage_drtr_mfe = (Viral_Phage[['drtr_fn_mfe']])
fn_Viral_Phage_drtr_mfe_sum = boot(fn_Viral_Phage_drtr_mfe)


tp_Viral_Phage_drtr_prob = (Viral_Phage[['drtr_tp_prob']])
tp_Viral_Phage_drtr_prob_sum = boot(tp_Viral_Phage_drtr_prob)

fp_Viral_Phage_drtr_prob = (Viral_Phage[['drtr_fp_prob']])
fp_Viral_Phage_drtr_prob_sum = boot(fp_Viral_Phage_drtr_prob)

tn_Viral_Phage_drtr_prob = (Viral_Phage[['drtr_tn_prob']])
tn_Viral_Phage_drtr_prob_sum = boot(tn_Viral_Phage_drtr_prob)

fn_Viral_Phage_drtr_prob = (Viral_Phage[['drtr_fn_prob']])
fn_Viral_Phage_drtr_prob_sum = boot(fn_Viral_Phage_drtr_prob)


tp_Viral_Phage_drtr_bpp = (Viral_Phage[['drtr_tp']])
tp_Viral_Phage_drtr_bpp_sum = boot(tp_Viral_Phage_drtr_bpp)

fp_Viral_Phage_drtr_bpp = (Viral_Phage[['drtr_fp']])
fp_Viral_Phage_drtr_bpp_sum = boot(fp_Viral_Phage_drtr_bpp)

tn_Viral_Phage_drtr_bpp = (Viral_Phage[['drtr_tn']])
tn_Viral_Phage_drtr_bpp_sum = boot(tn_Viral_Phage_drtr_bpp)

fn_Viral_Phage_drtr_bpp = (Viral_Phage[['drtr_fn']])
fn_Viral_Phage_drtr_bpp_sum = boot(fn_Viral_Phage_drtr_bpp)

#TPR
tpr_Viral_Phage_fold_mfe = tpr(tp_Viral_Phage_fold_mfe_sum, fn_Viral_Phage_fold_mfe_sum)
# print(tpr_Viral_Phage_fold_mfe)

tpr_Viral_Phage_fold_bpp = tpr(tp_Viral_Phage_fold_bpp_sum, fn_Viral_Phage_fold_bpp_sum)
# print(tpr_Viral_Phage_fold_bpp)

tpr_Viral_Phage_drtr_mfe = tpr(tp_Viral_Phage_drtr_mfe_sum, fn_Viral_Phage_drtr_mfe_sum)
# print(tpr_Viral_Phage_drtr_mfe)

tpr_Viral_Phage_drtr_prob = tpr(tp_Viral_Phage_drtr_prob_sum, fn_Viral_Phage_drtr_prob_sum)
# print(tpr_Viral_Phage_drtr_prob)

tpr_Viral_Phage_drtr_bpp = tpr(tp_Viral_Phage_drtr_bpp_sum, fn_Viral_Phage_drtr_bpp_sum)
# print(tpr_Viral_Phage_drtr_bpp)


#PPV
ppv_Viral_Phage_fold_mfe = ppv(tp_Viral_Phage_fold_mfe_sum, fn_Viral_Phage_fold_mfe_sum)
# print(ppv_Viral_Phage_fold_mfe)

ppv_Viral_Phage_fold_bpp = ppv(tp_Viral_Phage_fold_bpp_sum, fn_Viral_Phage_fold_bpp_sum)
# print(ppv_Viral_Phage_fold_bpp)

ppv_Viral_Phage_drtr_mfe = ppv(tp_Viral_Phage_drtr_mfe_sum, fn_Viral_Phage_drtr_mfe_sum)
# print(ppv_Viral_Phage_drtr_mfe)

ppv_Viral_Phage_drtr_prob = ppv(tp_Viral_Phage_drtr_prob_sum, fn_Viral_Phage_drtr_prob_sum)
# print(ppv_Viral_Phage_drtr_prob)

ppv_Viral_Phage_drtr_bpp = ppv(tp_Viral_Phage_drtr_bpp_sum, fn_Viral_Phage_drtr_bpp_sum)
# print(ppv_Viral_Phage_drtr_bpp)

#F1
f_1_Viral_Phage_fold_mfe = f_1(tp_Viral_Phage_fold_mfe_sum, fp_Viral_Phage_fold_mfe_sum, fn_Viral_Phage_fold_mfe_sum)
# print(f_1_Viral_Phage_fold_mfe)

f_1_Viral_Phage_fold_bpp = f_1(tp_Viral_Phage_fold_bpp_sum, fp_Viral_Phage_fold_bpp_sum, fn_Viral_Phage_fold_bpp_sum)
# print(f_1_Viral_Phage_fold_bpp)

f_1_Viral_Phage_drtr_mfe = f_1(tp_Viral_Phage_drtr_mfe_sum, fp_Viral_Phage_drtr_mfe_sum, fn_Viral_Phage_drtr_mfe_sum)
# print(f_1_Viral_Phage_drtr_mfe)

f_1_Viral_Phage_drtr_prob = f_1(tp_Viral_Phage_drtr_prob_sum, fp_Viral_Phage_drtr_prob_sum, fn_Viral_Phage_drtr_prob_sum)
# print(f_1_Viral_Phage_drtr_prob)

f_1_Viral_Phage_drtr_bpp = f_1(tp_Viral_Phage_drtr_bpp_sum, fp_Viral_Phage_drtr_bpp_sum, fn_Viral_Phage_drtr_bpp_sum)
# print(f_1_Viral_Phage_drtr_bpp)

# #MCC
mcc_Viral_Phage_fold_mfe = mcc(tp_Viral_Phage_fold_mfe_sum, fp_Viral_Phage_fold_mfe_sum, tn_Viral_Phage_fold_mfe_sum, fn_Viral_Phage_fold_mfe_sum)
# print(mcc_Viral_Phage_fold_mfe)

mcc_Viral_Phage_fold_bpp = mcc(tp_Viral_Phage_fold_bpp_sum, fp_Viral_Phage_fold_bpp_sum, tn_Viral_Phage_fold_bpp_sum, fn_Viral_Phage_fold_bpp_sum)
# print(mcc_Viral_Phage_fold_bpp)

mcc_Viral_Phage_drtr_mfe = mcc(tp_Viral_Phage_drtr_mfe_sum, fp_Viral_Phage_drtr_mfe_sum, tn_Viral_Phage_drtr_mfe_sum, fn_Viral_Phage_drtr_mfe_sum)
# print(mcc_Viral_Phage_drtr_mfe)

mcc_Viral_Phage_drtr_prob = mcc(tp_Viral_Phage_drtr_prob_sum, fp_Viral_Phage_drtr_prob_sum, tn_Viral_Phage_drtr_prob_sum, fn_Viral_Phage_drtr_prob_sum)
# print(mcc_Viral_Phage_drtr_prob)

mcc_Viral_Phage_drtr_bpp = mcc(tp_Viral_Phage_drtr_bpp_sum, fp_Viral_Phage_drtr_bpp_sum, tn_Viral_Phage_drtr_bpp_sum, fn_Viral_Phage_drtr_bpp_sum)
# print(mcc_Viral_Phage_drtr_bpp)


#YRNA
YRNA=family_data.loc[family_data.family=='YRNA']

tp_YRNA_fold_mfe = (YRNA[['RNAfold_mfe_tp']])
tp_YRNA_fold_mfe_sum = boot(tp_YRNA_fold_mfe)

fp_YRNA_fold_mfe = (YRNA[['RNAfold_mfe_fp']])
fp_YRNA_fold_mfe_sum = boot(fp_YRNA_fold_mfe)

tn_YRNA_fold_mfe = (YRNA[['RNAfold_mfe_tn']])
tn_YRNA_fold_mfe_sum = boot(tn_YRNA_fold_mfe)

fn_YRNA_fold_mfe = (YRNA[['RNAfold_mfe_fn']])
fn_YRNA_fold_mfe_sum = boot(fn_YRNA_fold_mfe)


tp_YRNA_fold_bpp = (YRNA[['RNAfold_bpp_tp']])
tp_YRNA_fold_bpp_sum = boot(tp_YRNA_fold_bpp)

fp_YRNA_fold_bpp = (YRNA[['RNAfold_bpp_fp']])
fp_YRNA_fold_bpp_sum = boot(fp_YRNA_fold_bpp)

tn_YRNA_fold_bpp = (YRNA[['RNAfold_bpp_tn']])
tn_YRNA_fold_bpp_sum = boot(tn_YRNA_fold_bpp)

fn_YRNA_fold_bpp = (YRNA[['RNAfold_bpp_fn']])
fn_YRNA_fold_bpp_sum = boot(fn_YRNA_fold_bpp)


tp_YRNA_drtr_mfe = (YRNA[['drtr_tp_mfe']])
tp_YRNA_drtr_mfe_sum = boot(tp_YRNA_drtr_mfe)

fp_YRNA_drtr_mfe = (YRNA[['drtr_fp_mfe']])
fp_YRNA_drtr_mfe_sum = boot(fp_YRNA_drtr_mfe)

tn_YRNA_drtr_mfe = (YRNA[['drtr_tn_mfe']])
tn_YRNA_drtr_mfe_sum = boot(tn_YRNA_drtr_mfe)

fn_YRNA_drtr_mfe = (YRNA[['drtr_fn_mfe']])
fn_YRNA_drtr_mfe_sum = boot(fn_YRNA_drtr_mfe)


tp_YRNA_drtr_prob = (YRNA[['drtr_tp_prob']])
tp_YRNA_drtr_prob_sum = boot(tp_YRNA_drtr_prob)

fp_YRNA_drtr_prob = (YRNA[['drtr_fp_prob']])
fp_YRNA_drtr_prob_sum = boot(fp_YRNA_drtr_prob)

tn_YRNA_drtr_prob = (YRNA[['drtr_tn_prob']])
tn_YRNA_drtr_prob_sum = boot(tn_YRNA_drtr_prob)

fn_YRNA_drtr_prob = (YRNA[['drtr_fn_prob']])
fn_YRNA_drtr_prob_sum = boot(fn_YRNA_drtr_prob)


tp_YRNA_drtr_bpp = (YRNA[['drtr_tp']])
tp_YRNA_drtr_bpp_sum = boot(tp_YRNA_drtr_bpp)

fp_YRNA_drtr_bpp = (YRNA[['drtr_fp']])
fp_YRNA_drtr_bpp_sum = boot(fp_YRNA_drtr_bpp)

tn_YRNA_drtr_bpp = (YRNA[['drtr_tn']])
tn_YRNA_drtr_bpp_sum = boot(tn_YRNA_drtr_bpp)

fn_YRNA_drtr_bpp = (YRNA[['drtr_fn']])
fn_YRNA_drtr_bpp_sum = boot(fn_YRNA_drtr_bpp)

#TPR
tpr_YRNA_fold_mfe = tpr(tp_YRNA_fold_mfe_sum, fn_YRNA_fold_mfe_sum)
# print(tpr_YRNA_fold_mfe)

tpr_YRNA_fold_bpp = tpr(tp_YRNA_fold_bpp_sum, fn_YRNA_fold_bpp_sum)
# print(tpr_YRNA_fold_bpp)

tpr_YRNA_drtr_mfe = tpr(tp_YRNA_drtr_mfe_sum, fn_YRNA_drtr_mfe_sum)
# print(tpr_YRNA_drtr_mfe)

tpr_YRNA_drtr_prob = tpr(tp_YRNA_drtr_prob_sum, fn_YRNA_drtr_prob_sum)
# print(tpr_YRNA_drtr_prob)

tpr_YRNA_drtr_bpp = tpr(tp_YRNA_drtr_bpp_sum, fn_YRNA_drtr_bpp_sum)
# print(tpr_YRNA_drtr_bpp)


#PPV
ppv_YRNA_fold_mfe = ppv(tp_YRNA_fold_mfe_sum, fn_YRNA_fold_mfe_sum)
# print(ppv_YRNA_fold_mfe)

ppv_YRNA_fold_bpp = ppv(tp_YRNA_fold_bpp_sum, fn_YRNA_fold_bpp_sum)
# print(ppv_YRNA_fold_bpp)

ppv_YRNA_drtr_mfe = ppv(tp_YRNA_drtr_mfe_sum, fn_YRNA_drtr_mfe_sum)
# print(ppv_YRNA_drtr_mfe)

ppv_YRNA_drtr_prob = ppv(tp_YRNA_drtr_prob_sum, fn_YRNA_drtr_prob_sum)
# print(ppv_YRNA_drtr_prob)

ppv_YRNA_drtr_bpp = ppv(tp_YRNA_drtr_bpp_sum, fn_YRNA_drtr_bpp_sum)
# print(ppv_YRNA_drtr_bpp)

#F1
f_1_YRNA_fold_mfe = f_1(tp_YRNA_fold_mfe_sum, fp_YRNA_fold_mfe_sum, fn_YRNA_fold_mfe_sum)
# print(f_1_YRNA_fold_mfe)

f_1_YRNA_fold_bpp = f_1(tp_YRNA_fold_bpp_sum, fp_YRNA_fold_bpp_sum, fn_YRNA_fold_bpp_sum)
# print(f_1_YRNA_fold_bpp)

f_1_YRNA_drtr_mfe = f_1(tp_YRNA_drtr_mfe_sum, fp_YRNA_drtr_mfe_sum, fn_YRNA_drtr_mfe_sum)
# print(f_1_YRNA_drtr_mfe)

f_1_YRNA_drtr_prob = f_1(tp_YRNA_drtr_prob_sum, fp_YRNA_drtr_prob_sum, fn_YRNA_drtr_prob_sum)
# print(f_1_YRNA_drtr_prob)

f_1_YRNA_drtr_bpp = f_1(tp_YRNA_drtr_bpp_sum, fp_YRNA_drtr_bpp_sum, fn_YRNA_drtr_bpp_sum)
# print(f_1_YRNA_drtr_bpp)

# #MCC
mcc_YRNA_fold_mfe = mcc(tp_YRNA_fold_mfe_sum, fp_YRNA_fold_mfe_sum, tn_YRNA_fold_mfe_sum, fn_YRNA_fold_mfe_sum)
# print(mcc_YRNA_fold_mfe)

mcc_YRNA_fold_bpp = mcc(tp_YRNA_fold_bpp_sum, fp_YRNA_fold_bpp_sum, tn_YRNA_fold_bpp_sum, fn_YRNA_fold_bpp_sum)
# print(mcc_YRNA_fold_bpp)

mcc_YRNA_drtr_mfe = mcc(tp_YRNA_drtr_mfe_sum, fp_YRNA_drtr_mfe_sum, tn_YRNA_drtr_mfe_sum, fn_YRNA_drtr_mfe_sum)
# print(mcc_YRNA_drtr_mfe)

mcc_YRNA_drtr_prob = mcc(tp_YRNA_drtr_prob_sum, fp_YRNA_drtr_prob_sum, tn_YRNA_drtr_prob_sum, fn_YRNA_drtr_prob_sum)
# print(mcc_YRNA_drtr_prob)

mcc_YRNA_drtr_bpp = mcc(tp_YRNA_drtr_bpp_sum, fp_YRNA_drtr_bpp_sum, tn_YRNA_drtr_bpp_sum, fn_YRNA_drtr_bpp_sum)
# print(mcc_YRNA_drtr_bpp)



if args.output == True:

    data_processed = [
    tpr_sixteen_SrRNA_fold_mfe,
    tpr_sixteen_SrRNA_fold_bpp,
    tpr_sixteen_SrRNA_drtr_mfe,
    tpr_sixteen_SrRNA_drtr_prob,
    tpr_sixteen_SrRNA_drtr_bpp,
    ppv_sixteen_SrRNA_fold_mfe,
    ppv_sixteen_SrRNA_fold_bpp,
    ppv_sixteen_SrRNA_drtr_mfe,
    ppv_sixteen_SrRNA_drtr_prob,
    ppv_sixteen_SrRNA_drtr_bpp,
    f_1_sixteen_SrRNA_fold_mfe,
    f_1_sixteen_SrRNA_fold_bpp,
    f_1_sixteen_SrRNA_drtr_mfe,
    f_1_sixteen_SrRNA_drtr_prob,
    f_1_sixteen_SrRNA_drtr_bpp,
    mcc_sixteen_SrRNA_fold_mfe,
    mcc_sixteen_SrRNA_fold_bpp,
    mcc_sixteen_SrRNA_drtr_mfe,
    mcc_sixteen_SrRNA_drtr_prob,
    mcc_sixteen_SrRNA_drtr_bpp,
    tpr_twentythree_SrRNA_fold_mfe,
    tpr_twentythree_SrRNA_fold_bpp,
    tpr_twentythree_SrRNA_drtr_mfe,
    tpr_twentythree_SrRNA_drtr_prob,
    tpr_twentythree_SrRNA_drtr_bpp,
    ppv_twentythree_SrRNA_fold_mfe,
    ppv_twentythree_SrRNA_fold_bpp,
    ppv_twentythree_SrRNA_drtr_mfe,
    ppv_twentythree_SrRNA_drtr_prob,
    ppv_twentythree_SrRNA_drtr_bpp,
    f_1_twentythree_SrRNA_fold_mfe,
    f_1_twentythree_SrRNA_fold_bpp,
    f_1_twentythree_SrRNA_drtr_mfe,
    f_1_twentythree_SrRNA_drtr_prob,
    f_1_twentythree_SrRNA_drtr_bpp,
    mcc_twentythree_SrRNA_fold_mfe,
    mcc_twentythree_SrRNA_fold_bpp,
    mcc_twentythree_SrRNA_drtr_mfe,
    mcc_twentythree_SrRNA_drtr_prob,
    mcc_twentythree_SrRNA_drtr_bpp,
    tpr_five_SrRNA_fold_mfe,
    tpr_five_SrRNA_fold_bpp,
    tpr_five_SrRNA_drtr_mfe,
    tpr_five_SrRNA_drtr_prob,
    tpr_five_SrRNA_drtr_bpp,
    ppv_five_SrRNA_fold_mfe,
    ppv_five_SrRNA_fold_bpp,
    ppv_five_SrRNA_drtr_mfe,
    ppv_five_SrRNA_drtr_prob,
    ppv_five_SrRNA_drtr_bpp,
    f_1_five_SrRNA_fold_mfe,
    f_1_five_SrRNA_fold_bpp,
    f_1_five_SrRNA_drtr_mfe,
    f_1_five_SrRNA_drtr_prob,
    f_1_five_SrRNA_drtr_bpp,
    mcc_five_SrRNA_fold_mfe,
    mcc_five_SrRNA_fold_bpp,
    mcc_five_SrRNA_drtr_mfe,
    mcc_five_SrRNA_drtr_prob,
    mcc_five_SrRNA_drtr_bpp,
    tpr_cili_telo_RNA_fold_mfe,
    tpr_cili_telo_RNA_fold_bpp,
    tpr_cili_telo_RNA_drtr_mfe,
    tpr_cili_telo_RNA_drtr_prob,
    tpr_cili_telo_RNA_drtr_bpp,
    ppv_cili_telo_RNA_fold_mfe,
    ppv_cili_telo_RNA_fold_bpp,
    ppv_cili_telo_RNA_drtr_mfe,
    ppv_cili_telo_RNA_drtr_prob,
    ppv_cili_telo_RNA_drtr_bpp,
    f_1_cili_telo_RNA_fold_mfe,
    f_1_cili_telo_RNA_fold_bpp,
    f_1_cili_telo_RNA_drtr_mfe,
    f_1_cili_telo_RNA_drtr_prob,
    f_1_cili_telo_RNA_drtr_bpp,
    mcc_cili_telo_RNA_fold_mfe,
    mcc_cili_telo_RNA_fold_bpp,
    mcc_cili_telo_RNA_drtr_mfe,
    mcc_cili_telo_RNA_drtr_prob,
    mcc_cili_telo_RNA_drtr_bpp,
    tpr_cis_regulatory_element_fold_mfe,
    tpr_cis_regulatory_element_fold_bpp,
    tpr_cis_regulatory_element_drtr_mfe,
    tpr_cis_regulatory_element_drtr_prob,
    tpr_cis_regulatory_element_drtr_bpp,
    ppv_cis_regulatory_element_fold_mfe,
    ppv_cis_regulatory_element_fold_bpp,
    ppv_cis_regulatory_element_drtr_mfe,
    ppv_cis_regulatory_element_drtr_prob,
    ppv_cis_regulatory_element_drtr_bpp,
    f_1_cis_regulatory_element_fold_mfe,
    f_1_cis_regulatory_element_fold_bpp,
    f_1_cis_regulatory_element_drtr_mfe,
    f_1_cis_regulatory_element_drtr_prob,
    f_1_cis_regulatory_element_drtr_bpp,
    mcc_cis_regulatory_element_fold_mfe,
    mcc_cis_regulatory_element_fold_bpp,
    mcc_cis_regulatory_element_drtr_mfe,
    mcc_cis_regulatory_element_drtr_prob,
    mcc_cis_regulatory_element_drtr_bpp,
    tpr_GIIIntron_fold_mfe,
    tpr_GIIIntron_fold_bpp,
    tpr_GIIIntron_drtr_mfe,
    tpr_GIIIntron_drtr_prob,
    tpr_GIIIntron_drtr_bpp,
    ppv_GIIIntron_fold_mfe,
    ppv_GIIIntron_fold_bpp,
    ppv_GIIIntron_drtr_mfe,
    ppv_GIIIntron_drtr_prob,
    ppv_GIIIntron_drtr_bpp,
    f_1_GIIIntron_fold_mfe,
    f_1_GIIIntron_fold_bpp,
    f_1_GIIIntron_drtr_mfe,
    f_1_GIIIntron_drtr_prob,
    f_1_GIIIntron_drtr_bpp,
    mcc_GIIIntron_fold_mfe,
    mcc_GIIIntron_fold_bpp,
    mcc_GIIIntron_drtr_mfe,
    mcc_GIIIntron_drtr_prob,
    mcc_GIIIntron_drtr_bpp,
    tpr_GIIntron_fold_mfe,
    tpr_GIIntron_fold_bpp,
    tpr_GIIntron_drtr_mfe,
    tpr_GIIntron_drtr_prob,
    tpr_GIIntron_drtr_bpp,
    ppv_GIIntron_fold_mfe,
    ppv_GIIntron_fold_bpp,
    ppv_GIIntron_drtr_mfe,
    ppv_GIIntron_drtr_prob,
    ppv_GIIntron_drtr_bpp,
    f_1_GIIntron_fold_mfe,
    f_1_GIIntron_fold_bpp,
    f_1_GIIntron_drtr_mfe,
    f_1_GIIntron_drtr_prob,
    f_1_GIIntron_drtr_bpp,
    mcc_GIIntron_fold_mfe,
    mcc_GIIntron_fold_bpp,
    mcc_GIIntron_drtr_mfe,
    mcc_GIIntron_drtr_prob,
    mcc_GIIntron_drtr_bpp,
    tpr_ham_ribozyme_fold_mfe,
    tpr_ham_ribozyme_fold_bpp,
    tpr_ham_ribozyme_drtr_mfe,
    tpr_ham_ribozyme_drtr_prob,
    tpr_ham_ribozyme_drtr_bpp,
    ppv_ham_ribozyme_fold_mfe,
    ppv_ham_ribozyme_fold_bpp,
    ppv_ham_ribozyme_drtr_mfe,
    ppv_ham_ribozyme_drtr_prob,
    ppv_ham_ribozyme_drtr_bpp,
    f_1_ham_ribozyme_fold_mfe,
    f_1_ham_ribozyme_fold_bpp,
    f_1_ham_ribozyme_drtr_mfe,
    f_1_ham_ribozyme_drtr_prob,
    f_1_ham_ribozyme_drtr_bpp,
    mcc_ham_ribozyme_fold_mfe,
    mcc_ham_ribozyme_fold_bpp,
    mcc_ham_ribozyme_drtr_mfe,
    mcc_ham_ribozyme_drtr_prob,
    mcc_ham_ribozyme_drtr_bpp,
    tpr_hdvr_ribozyme_fold_mfe,
    tpr_hdvr_ribozyme_fold_bpp,
    tpr_hdvr_ribozyme_drtr_mfe,
    tpr_hdvr_ribozyme_drtr_prob,
    tpr_hdvr_ribozyme_drtr_bpp,
    ppv_hdvr_ribozyme_fold_mfe,
    ppv_hdvr_ribozyme_fold_bpp,
    ppv_hdvr_ribozyme_drtr_mfe,
    ppv_hdvr_ribozyme_drtr_prob,
    ppv_hdvr_ribozyme_drtr_bpp,
    f_1_hdvr_ribozyme_fold_mfe,
    f_1_hdvr_ribozyme_fold_bpp,
    f_1_hdvr_ribozyme_drtr_mfe,
    f_1_hdvr_ribozyme_drtr_prob,
    f_1_hdvr_ribozyme_drtr_bpp,
    mcc_hdvr_ribozyme_fold_mfe,
    mcc_hdvr_ribozyme_fold_bpp,
    mcc_hdvr_ribozyme_drtr_mfe,
    mcc_hdvr_ribozyme_drtr_prob,
    mcc_hdvr_ribozyme_drtr_bpp,
    tpr_ires_fold_mfe,
    tpr_ires_fold_bpp,
    tpr_ires_drtr_mfe,
    tpr_ires_drtr_prob,
    tpr_ires_drtr_bpp,
    ppv_ires_fold_mfe,
    ppv_ires_fold_bpp,
    ppv_ires_drtr_mfe,
    ppv_ires_drtr_prob,
    ppv_ires_drtr_bpp,
    f_1_ires_fold_mfe,
    f_1_ires_fold_bpp,
    f_1_ires_drtr_mfe,
    f_1_ires_drtr_prob,
    f_1_ires_drtr_bpp,
    mcc_ires_fold_mfe,
    mcc_ires_fold_bpp,
    mcc_ires_drtr_mfe,
    mcc_ires_drtr_prob,
    mcc_ires_drtr_bpp,
    tpr_other_ribozyme_fold_mfe,
    tpr_other_ribozyme_fold_bpp,
    tpr_other_ribozyme_drtr_mfe,
    tpr_other_ribozyme_drtr_prob,
    tpr_other_ribozyme_drtr_bpp,
    ppv_other_ribozyme_fold_mfe,
    ppv_other_ribozyme_fold_bpp,
    ppv_other_ribozyme_drtr_mfe,
    ppv_other_ribozyme_drtr_prob,
    ppv_other_ribozyme_drtr_bpp,
    f_1_other_ribozyme_fold_mfe,
    f_1_other_ribozyme_fold_bpp,
    f_1_other_ribozyme_drtr_mfe,
    f_1_other_ribozyme_drtr_prob,
    f_1_other_ribozyme_drtr_bpp,
    mcc_other_ribozyme_fold_mfe,
    mcc_other_ribozyme_fold_bpp,
    mcc_other_ribozyme_drtr_mfe,
    mcc_other_ribozyme_drtr_prob,
    mcc_other_ribozyme_drtr_bpp,
    tpr_other_RNA_fold_mfe,
    tpr_other_RNA_fold_bpp,
    tpr_other_RNA_drtr_mfe,
    tpr_other_RNA_drtr_prob,
    tpr_other_RNA_drtr_bpp,
    ppv_other_RNA_fold_mfe,
    ppv_other_RNA_fold_bpp,
    ppv_other_RNA_drtr_mfe,
    ppv_other_RNA_drtr_prob,
    ppv_other_RNA_drtr_bpp,
    f_1_other_RNA_fold_mfe,
    f_1_other_RNA_fold_bpp,
    f_1_other_RNA_drtr_mfe,
    f_1_other_RNA_drtr_prob,
    f_1_other_RNA_drtr_bpp,
    mcc_other_RNA_fold_mfe,
    mcc_other_RNA_fold_bpp,
    mcc_other_RNA_drtr_mfe,
    mcc_other_RNA_drtr_prob,
    mcc_other_RNA_drtr_bpp,
    tpr_other_rRNA_fold_mfe,
    tpr_other_rRNA_fold_bpp,
    tpr_other_rRNA_drtr_mfe,
    tpr_other_rRNA_drtr_prob,
    tpr_other_rRNA_drtr_bpp,
    ppv_other_rRNA_fold_mfe,
    ppv_other_rRNA_fold_bpp,
    ppv_other_rRNA_drtr_mfe,
    ppv_other_rRNA_drtr_prob,
    ppv_other_rRNA_drtr_bpp,
    f_1_other_rRNA_fold_mfe,
    f_1_other_rRNA_fold_bpp,
    f_1_other_rRNA_drtr_mfe,
    f_1_other_rRNA_drtr_prob,
    f_1_other_rRNA_drtr_bpp,
    mcc_other_rRNA_fold_mfe,
    mcc_other_rRNA_fold_bpp,
    mcc_other_rRNA_drtr_mfe,
    mcc_other_rRNA_drtr_prob,
    mcc_other_rRNA_drtr_bpp,
    tpr_RNAIII_fold_mfe,
    tpr_RNAIII_fold_bpp,
    tpr_RNAIII_drtr_mfe,
    tpr_RNAIII_drtr_prob,
    tpr_RNAIII_drtr_bpp,
    ppv_RNAIII_fold_mfe,
    ppv_RNAIII_fold_bpp,
    ppv_RNAIII_drtr_mfe,
    ppv_RNAIII_drtr_prob,
    ppv_RNAIII_drtr_bpp,
    f_1_RNAIII_fold_mfe,
    f_1_RNAIII_fold_bpp,
    f_1_RNAIII_drtr_mfe,
    f_1_RNAIII_drtr_prob,
    f_1_RNAIII_drtr_bpp,
    mcc_RNAIII_fold_mfe,
    mcc_RNAIII_fold_bpp,
    mcc_RNAIII_drtr_mfe,
    mcc_RNAIII_drtr_prob,
    mcc_RNAIII_drtr_bpp,
    tpr_RNaseE5UTR_fold_mfe,
    tpr_RNaseE5UTR_fold_bpp,
    tpr_RNaseE5UTR_drtr_mfe,
    tpr_RNaseE5UTR_drtr_prob,
    tpr_RNaseE5UTR_drtr_bpp,
    ppv_RNaseE5UTR_fold_mfe,
    ppv_RNaseE5UTR_fold_bpp,
    ppv_RNaseE5UTR_drtr_mfe,
    ppv_RNaseE5UTR_drtr_prob,
    ppv_RNaseE5UTR_drtr_bpp,
    f_1_RNaseE5UTR_fold_mfe,
    f_1_RNaseE5UTR_fold_bpp,
    f_1_RNaseE5UTR_drtr_mfe,
    f_1_RNaseE5UTR_drtr_prob,
    f_1_RNaseE5UTR_drtr_bpp,
    mcc_RNaseE5UTR_fold_mfe,
    mcc_RNaseE5UTR_fold_bpp,
    mcc_RNaseE5UTR_drtr_mfe,
    mcc_RNaseE5UTR_drtr_prob,
    mcc_RNaseE5UTR_drtr_bpp,
    tpr_RNaseMRPRNA_fold_mfe,
    tpr_RNaseMRPRNA_fold_bpp,
    tpr_RNaseMRPRNA_drtr_mfe,
    tpr_RNaseMRPRNA_drtr_prob,
    tpr_RNaseMRPRNA_drtr_bpp,
    ppv_RNaseMRPRNA_fold_mfe,
    ppv_RNaseMRPRNA_fold_bpp,
    ppv_RNaseMRPRNA_drtr_mfe,
    ppv_RNaseMRPRNA_drtr_prob,
    ppv_RNaseMRPRNA_drtr_bpp,
    f_1_RNaseMRPRNA_fold_mfe,
    f_1_RNaseMRPRNA_fold_bpp,
    f_1_RNaseMRPRNA_drtr_mfe,
    f_1_RNaseMRPRNA_drtr_prob,
    f_1_RNaseMRPRNA_drtr_bpp,
    mcc_RNaseMRPRNA_fold_mfe,
    mcc_RNaseMRPRNA_fold_bpp,
    mcc_RNaseMRPRNA_drtr_mfe,
    mcc_RNaseMRPRNA_drtr_prob,
    mcc_RNaseMRPRNA_drtr_bpp,
    tpr_RNasePRNA_fold_mfe,
    tpr_RNasePRNA_fold_bpp,
    tpr_RNasePRNA_drtr_mfe,
    tpr_RNasePRNA_drtr_prob,
    tpr_RNasePRNA_drtr_bpp,
    ppv_RNasePRNA_fold_mfe,
    ppv_RNasePRNA_fold_bpp,
    ppv_RNasePRNA_drtr_mfe,
    ppv_RNasePRNA_drtr_prob,
    ppv_RNasePRNA_drtr_bpp,
    f_1_RNasePRNA_fold_mfe,
    f_1_RNasePRNA_fold_bpp,
    f_1_RNasePRNA_drtr_mfe,
    f_1_RNasePRNA_drtr_prob,
    f_1_RNasePRNA_drtr_bpp,
    mcc_RNasePRNA_fold_mfe,
    mcc_RNasePRNA_fold_bpp,
    mcc_RNasePRNA_drtr_mfe,
    mcc_RNasePRNA_drtr_prob,
    mcc_RNasePRNA_drtr_bpp,
    tpr_snRNA_fold_mfe,
    tpr_snRNA_fold_bpp,
    tpr_snRNA_drtr_mfe,
    tpr_snRNA_drtr_prob,
    tpr_snRNA_drtr_bpp,
    ppv_snRNA_fold_mfe,
    ppv_snRNA_fold_bpp,
    ppv_snRNA_drtr_mfe,
    ppv_snRNA_drtr_prob,
    ppv_snRNA_drtr_bpp,
    f_1_snRNA_fold_mfe,
    f_1_snRNA_fold_bpp,
    f_1_snRNA_drtr_mfe,
    f_1_snRNA_drtr_prob,
    f_1_snRNA_drtr_bpp,
    mcc_snRNA_fold_mfe,
    mcc_snRNA_fold_bpp,
    mcc_snRNA_drtr_mfe,
    mcc_snRNA_drtr_prob,
    mcc_snRNA_drtr_bpp,
    tpr_SRPRNA_fold_mfe,
    tpr_SRPRNA_fold_bpp,
    tpr_SRPRNA_drtr_mfe,
    tpr_SRPRNA_drtr_prob,
    tpr_SRPRNA_drtr_bpp,
    ppv_SRPRNA_fold_mfe,
    ppv_SRPRNA_fold_bpp,
    ppv_SRPRNA_drtr_mfe,
    ppv_SRPRNA_drtr_prob,
    ppv_SRPRNA_drtr_bpp,
    f_1_SRPRNA_fold_mfe,
    f_1_SRPRNA_fold_bpp,
    f_1_SRPRNA_drtr_mfe,
    f_1_SRPRNA_drtr_prob,
    f_1_SRPRNA_drtr_bpp,
    mcc_SRPRNA_fold_mfe,
    mcc_SRPRNA_fold_bpp,
    mcc_SRPRNA_drtr_mfe,
    mcc_SRPRNA_drtr_prob,
    mcc_SRPRNA_drtr_bpp,
    tpr_SyntheticRNA_fold_mfe,
    tpr_SyntheticRNA_fold_bpp,
    tpr_SyntheticRNA_drtr_mfe,
    tpr_SyntheticRNA_drtr_prob,
    tpr_SyntheticRNA_drtr_bpp,
    ppv_SyntheticRNA_fold_mfe,
    ppv_SyntheticRNA_fold_bpp,
    ppv_SyntheticRNA_drtr_mfe,
    ppv_SyntheticRNA_drtr_prob,
    ppv_SyntheticRNA_drtr_bpp,
    f_1_SyntheticRNA_fold_mfe,
    f_1_SyntheticRNA_fold_bpp,
    f_1_SyntheticRNA_drtr_mfe,
    f_1_SyntheticRNA_drtr_prob,
    f_1_SyntheticRNA_drtr_bpp,
    mcc_SyntheticRNA_fold_mfe,
    mcc_SyntheticRNA_fold_bpp,
    mcc_SyntheticRNA_drtr_mfe,
    mcc_SyntheticRNA_drtr_prob,
    mcc_SyntheticRNA_drtr_bpp,
    tpr_tmRNA_fold_mfe,
    tpr_tmRNA_fold_bpp,
    tpr_tmRNA_drtr_mfe,
    tpr_tmRNA_drtr_prob,
    tpr_tmRNA_drtr_bpp,
    ppv_tmRNA_fold_mfe,
    ppv_tmRNA_fold_bpp,
    ppv_tmRNA_drtr_mfe,
    ppv_tmRNA_drtr_prob,
    ppv_tmRNA_drtr_bpp,
    f_1_tmRNA_fold_mfe,
    f_1_tmRNA_fold_bpp,
    f_1_tmRNA_drtr_mfe,
    f_1_tmRNA_drtr_prob,
    f_1_tmRNA_drtr_bpp,
    mcc_tmRNA_fold_mfe,
    mcc_tmRNA_fold_bpp,
    mcc_tmRNA_drtr_mfe,
    mcc_tmRNA_drtr_prob,
    mcc_tmRNA_drtr_bpp,
    tpr_tRNA_fold_mfe,
    tpr_tRNA_fold_bpp,
    tpr_tRNA_drtr_mfe,
    tpr_tRNA_drtr_prob,
    tpr_tRNA_drtr_bpp,
    ppv_tRNA_fold_mfe,
    ppv_tRNA_fold_bpp,
    ppv_tRNA_drtr_mfe,
    ppv_tRNA_drtr_prob,
    ppv_tRNA_drtr_bpp,
    f_1_tRNA_fold_mfe,
    f_1_tRNA_fold_bpp,
    f_1_tRNA_drtr_mfe,
    f_1_tRNA_drtr_prob,
    f_1_tRNA_drtr_bpp,
    mcc_tRNA_fold_mfe,
    mcc_tRNA_fold_bpp,
    mcc_tRNA_drtr_mfe,
    mcc_tRNA_drtr_prob,
    mcc_tRNA_drtr_bpp,
    tpr_Vert_Telo_RNA_fold_mfe,
    tpr_Vert_Telo_RNA_fold_bpp,
    tpr_Vert_Telo_RNA_drtr_mfe,
    tpr_Vert_Telo_RNA_drtr_prob,
    tpr_Vert_Telo_RNA_drtr_bpp,
    ppv_Vert_Telo_RNA_fold_mfe,
    ppv_Vert_Telo_RNA_fold_bpp,
    ppv_Vert_Telo_RNA_drtr_mfe,
    ppv_Vert_Telo_RNA_drtr_prob,
    ppv_Vert_Telo_RNA_drtr_bpp,
    f_1_Vert_Telo_RNA_fold_mfe,
    f_1_Vert_Telo_RNA_fold_bpp,
    f_1_Vert_Telo_RNA_drtr_mfe,
    f_1_Vert_Telo_RNA_drtr_prob,
    f_1_Vert_Telo_RNA_drtr_bpp,
    mcc_Vert_Telo_RNA_fold_mfe,
    mcc_Vert_Telo_RNA_fold_bpp,
    mcc_Vert_Telo_RNA_drtr_mfe,
    mcc_Vert_Telo_RNA_drtr_prob,
    mcc_Vert_Telo_RNA_drtr_bpp,
    tpr_Viral_Phage_fold_mfe,
    tpr_Viral_Phage_fold_bpp,
    tpr_Viral_Phage_drtr_mfe,
    tpr_Viral_Phage_drtr_prob,
    tpr_Viral_Phage_drtr_bpp,
    ppv_Viral_Phage_fold_mfe,
    ppv_Viral_Phage_fold_bpp,
    ppv_Viral_Phage_drtr_mfe,
    ppv_Viral_Phage_drtr_prob,
    ppv_Viral_Phage_drtr_bpp,
    f_1_Viral_Phage_fold_mfe,
    f_1_Viral_Phage_fold_bpp,
    f_1_Viral_Phage_drtr_mfe,
    f_1_Viral_Phage_drtr_prob,
    f_1_Viral_Phage_drtr_bpp,
    mcc_Viral_Phage_fold_mfe,
    mcc_Viral_Phage_fold_bpp,
    mcc_Viral_Phage_drtr_mfe,
    mcc_Viral_Phage_drtr_prob,
    mcc_Viral_Phage_drtr_bpp,
    tpr_YRNA_fold_mfe,
    tpr_YRNA_fold_bpp,
    tpr_YRNA_drtr_mfe,
    tpr_YRNA_drtr_prob,
    tpr_YRNA_drtr_bpp,
    ppv_YRNA_fold_mfe,
    ppv_YRNA_fold_bpp,
    ppv_YRNA_drtr_mfe,
    ppv_YRNA_drtr_prob,
    ppv_YRNA_drtr_bpp,
    f_1_YRNA_fold_mfe,
    f_1_YRNA_fold_bpp,
    f_1_YRNA_drtr_mfe,
    f_1_YRNA_drtr_prob,
    f_1_YRNA_drtr_bpp,
    mcc_YRNA_fold_mfe,
    mcc_YRNA_fold_bpp,
    mcc_YRNA_drtr_mfe,
    mcc_YRNA_drtr_prob,
    mcc_YRNA_drtr_bpp
    ]


    indexer = [
    "tpr_sixteen_SrRNA_fold_mfe",
    "tpr_sixteen_SrRNA_fold_bpp",
    "tpr_sixteen_SrRNA_drtr_mfe",
    "tpr_sixteen_SrRNA_drtr_prob",
    "tpr_sixteen_SrRNA_drtr_bpp",
    "ppv_sixteen_SrRNA_fold_mfe",
    "ppv_sixteen_SrRNA_fold_bpp",
    "ppv_sixteen_SrRNA_drtr_mfe",
    "ppv_sixteen_SrRNA_drtr_prob",
    "ppv_sixteen_SrRNA_drtr_bpp",
    "f_1_sixteen_SrRNA_fold_mfe",
    "f_1_sixteen_SrRNA_fold_bpp",
    "f_1_sixteen_SrRNA_drtr_mfe",
    "f_1_sixteen_SrRNA_drtr_prob",
    "f_1_sixteen_SrRNA_drtr_bpp",
    "mcc_sixteen_SrRNA_fold_mfe",
    "mcc_sixteen_SrRNA_fold_bpp",
    "mcc_sixteen_SrRNA_drtr_mfe",
    "mcc_sixteen_SrRNA_drtr_prob",
    "mcc_sixteen_SrRNA_drtr_bpp",
    "tpr_twentythree_SrRNA_fold_mfe",
    "tpr_twentythree_SrRNA_fold_bpp",
    "tpr_twentythree_SrRNA_drtr_mfe",
    "tpr_twentythree_SrRNA_drtr_prob",
    "tpr_twentythree_SrRNA_drtr_bpp",
    "ppv_twentythree_SrRNA_fold_mfe",
    "ppv_twentythree_SrRNA_fold_bpp",
    "ppv_twentythree_SrRNA_drtr_mfe",
    "ppv_twentythree_SrRNA_drtr_prob",
    "ppv_twentythree_SrRNA_drtr_bpp",
    "f_1_twentythree_SrRNA_fold_mfe",
    "f_1_twentythree_SrRNA_fold_bpp",
    "f_1_twentythree_SrRNA_drtr_mfe",
    "f_1_twentythree_SrRNA_drtr_prob",
    "f_1_twentythree_SrRNA_drtr_bpp",
    "mcc_twentythree_SrRNA_fold_mfe",
    "mcc_twentythree_SrRNA_fold_bpp",
    "mcc_twentythree_SrRNA_drtr_mfe",
    "mcc_twentythree_SrRNA_drtr_prob",
    "mcc_twentythree_SrRNA_drtr_bpp",
    "tpr_five_SrRNA_fold_mfe",
    "tpr_five_SrRNA_fold_bpp",
    "tpr_five_SrRNA_drtr_mfe",
    "tpr_five_SrRNA_drtr_prob",
    "tpr_five_SrRNA_drtr_bpp",
    "ppv_five_SrRNA_fold_mfe",
    "ppv_five_SrRNA_fold_bpp",
    "ppv_five_SrRNA_drtr_mfe",
    "ppv_five_SrRNA_drtr_prob",
    "ppv_five_SrRNA_drtr_bpp",
    "f_1_five_SrRNA_fold_mfe",
    "f_1_five_SrRNA_fold_bpp",
    "f_1_five_SrRNA_drtr_mfe",
    "f_1_five_SrRNA_drtr_prob",
    "f_1_five_SrRNA_drtr_bpp",
    "mcc_five_SrRNA_fold_mfe",
    "mcc_five_SrRNA_fold_bpp",
    "mcc_five_SrRNA_drtr_mfe",
    "mcc_five_SrRNA_drtr_prob",
    "mcc_five_SrRNA_drtr_bpp",
    "tpr_cili_telo_RNA_fold_mfe",
    "tpr_cili_telo_RNA_fold_bpp",
    "tpr_cili_telo_RNA_drtr_mfe",
    "tpr_cili_telo_RNA_drtr_prob",
    "tpr_cili_telo_RNA_drtr_bpp",
    "ppv_cili_telo_RNA_fold_mfe",
    "ppv_cili_telo_RNA_fold_bpp",
    "ppv_cili_telo_RNA_drtr_mfe",
    "ppv_cili_telo_RNA_drtr_prob",
    "ppv_cili_telo_RNA_drtr_bpp",
    "f_1_cili_telo_RNA_fold_mfe",
    "f_1_cili_telo_RNA_fold_bpp",
    "f_1_cili_telo_RNA_drtr_mfe",
    "f_1_cili_telo_RNA_drtr_prob",
    "f_1_cili_telo_RNA_drtr_bpp",
    "mcc_cili_telo_RNA_fold_mfe",
    "mcc_cili_telo_RNA_fold_bpp",
    "mcc_cili_telo_RNA_drtr_mfe",
    "mcc_cili_telo_RNA_drtr_prob",
    "mcc_cili_telo_RNA_drtr_bpp",
    "tpr_cis_regulatory_element_fold_mfe",
    "tpr_cis_regulatory_element_fold_bpp",
    "tpr_cis_regulatory_element_drtr_mfe",
    "tpr_cis_regulatory_element_drtr_prob",
    "tpr_cis_regulatory_element_drtr_bpp",
    "ppv_cis_regulatory_element_fold_mfe",
    "ppv_cis_regulatory_element_fold_bpp",
    "ppv_cis_regulatory_element_drtr_mfe",
    "ppv_cis_regulatory_element_drtr_prob",
    "ppv_cis_regulatory_element_drtr_bpp",
    "f_1_cis_regulatory_element_fold_mfe",
    "f_1_cis_regulatory_element_fold_bpp",
    "f_1_cis_regulatory_element_drtr_mfe",
    "f_1_cis_regulatory_element_drtr_prob",
    "f_1_cis_regulatory_element_drtr_bpp",
    "mcc_cis_regulatory_element_fold_mfe",
    "mcc_cis_regulatory_element_fold_bpp",
    "mcc_cis_regulatory_element_drtr_mfe",
    "mcc_cis_regulatory_element_drtr_prob",
    "mcc_cis_regulatory_element_drtr_bpp",
    "tpr_GIIIntron_fold_mfe",
    "tpr_GIIIntron_fold_bpp",
    "tpr_GIIIntron_drtr_mfe",
    "tpr_GIIIntron_drtr_prob",
    "tpr_GIIIntron_drtr_bpp",
    "ppv_GIIIntron_fold_mfe",
    "ppv_GIIIntron_fold_bpp",
    "ppv_GIIIntron_drtr_mfe",
    "ppv_GIIIntron_drtr_prob",
    "ppv_GIIIntron_drtr_bpp",
    "f_1_GIIIntron_fold_mfe",
    "f_1_GIIIntron_fold_bpp",
    "f_1_GIIIntron_drtr_mfe",
    "f_1_GIIIntron_drtr_prob",
    "f_1_GIIIntron_drtr_bpp",
    "mcc_GIIIntron_fold_mfe",
    "mcc_GIIIntron_fold_bpp",
    "mcc_GIIIntron_drtr_mfe",
    "mcc_GIIIntron_drtr_prob",
    "mcc_GIIIntron_drtr_bpp",
    "tpr_GIIntron_fold_mfe",
    "tpr_GIIntron_fold_bpp",
    "tpr_GIIntron_drtr_mfe",
    "tpr_GIIntron_drtr_prob",
    "tpr_GIIntron_drtr_bpp",
    "ppv_GIIntron_fold_mfe",
    "ppv_GIIntron_fold_bpp",
    "ppv_GIIntron_drtr_mfe",
    "ppv_GIIntron_drtr_prob",
    "ppv_GIIntron_drtr_bpp",
    "f_1_GIIntron_fold_mfe",
    "f_1_GIIntron_fold_bpp",
    "f_1_GIIntron_drtr_mfe",
    "f_1_GIIntron_drtr_prob",
    "f_1_GIIntron_drtr_bpp",
    "mcc_GIIntron_fold_mfe",
    "mcc_GIIntron_fold_bpp",
    "mcc_GIIntron_drtr_mfe",
    "mcc_GIIntron_drtr_prob",
    "mcc_GIIntron_drtr_bpp",
    "tpr_ham_ribozyme_fold_mfe",
    "tpr_ham_ribozyme_fold_bpp",
    "tpr_ham_ribozyme_drtr_mfe",
    "tpr_ham_ribozyme_drtr_prob",
    "tpr_ham_ribozyme_drtr_bpp",
    "ppv_ham_ribozyme_fold_mfe",
    "ppv_ham_ribozyme_fold_bpp",
    "ppv_ham_ribozyme_drtr_mfe",
    "ppv_ham_ribozyme_drtr_prob",
    "ppv_ham_ribozyme_drtr_bpp",
    "f_1_ham_ribozyme_fold_mfe",
    "f_1_ham_ribozyme_fold_bpp",
    "f_1_ham_ribozyme_drtr_mfe",
    "f_1_ham_ribozyme_drtr_prob",
    "f_1_ham_ribozyme_drtr_bpp",
    "mcc_ham_ribozyme_fold_mfe",
    "mcc_ham_ribozyme_fold_bpp",
    "mcc_ham_ribozyme_drtr_mfe",
    "mcc_ham_ribozyme_drtr_prob",
    "mcc_ham_ribozyme_drtr_bpp",
    "tpr_hdvr_ribozyme_fold_mfe",
    "tpr_hdvr_ribozyme_fold_bpp",
    "tpr_hdvr_ribozyme_drtr_mfe",
    "tpr_hdvr_ribozyme_drtr_prob",
    "tpr_hdvr_ribozyme_drtr_bpp",
    "ppv_hdvr_ribozyme_fold_mfe",
    "ppv_hdvr_ribozyme_fold_bpp",
    "ppv_hdvr_ribozyme_drtr_mfe",
    "ppv_hdvr_ribozyme_drtr_prob",
    "ppv_hdvr_ribozyme_drtr_bpp",
    "f_1_hdvr_ribozyme_fold_mfe",
    "f_1_hdvr_ribozyme_fold_bpp",
    "f_1_hdvr_ribozyme_drtr_mfe",
    "f_1_hdvr_ribozyme_drtr_prob",
    "f_1_hdvr_ribozyme_drtr_bpp",
    "mcc_hdvr_ribozyme_fold_mfe",
    "mcc_hdvr_ribozyme_fold_bpp",
    "mcc_hdvr_ribozyme_drtr_mfe",
    "mcc_hdvr_ribozyme_drtr_prob",
    "mcc_hdvr_ribozyme_drtr_bpp",
    "tpr_ires_fold_mfe",
    "tpr_ires_fold_bpp",
    "tpr_ires_drtr_mfe",
    "tpr_ires_drtr_prob",
    "tpr_ires_drtr_bpp",
    "ppv_ires_fold_mfe",
    "ppv_ires_fold_bpp",
    "ppv_ires_drtr_mfe",
    "ppv_ires_drtr_prob",
    "ppv_ires_drtr_bpp",
    "f_1_ires_fold_mfe",
    "f_1_ires_fold_bpp",
    "f_1_ires_drtr_mfe",
    "f_1_ires_drtr_prob",
    "f_1_ires_drtr_bpp",
    "mcc_ires_fold_mfe",
    "mcc_ires_fold_bpp",
    "mcc_ires_drtr_mfe",
    "mcc_ires_drtr_prob",
    "mcc_ires_drtr_bpp",
    "tpr_other_ribozyme_fold_mfe",
    "tpr_other_ribozyme_fold_bpp",
    "tpr_other_ribozyme_drtr_mfe",
    "tpr_other_ribozyme_drtr_prob",
    "tpr_other_ribozyme_drtr_bpp",
    "ppv_other_ribozyme_fold_mfe",
    "ppv_other_ribozyme_fold_bpp",
    "ppv_other_ribozyme_drtr_mfe",
    "ppv_other_ribozyme_drtr_prob",
    "ppv_other_ribozyme_drtr_bpp",
    "f_1_other_ribozyme_fold_mfe",
    "f_1_other_ribozyme_fold_bpp",
    "f_1_other_ribozyme_drtr_mfe",
    "f_1_other_ribozyme_drtr_prob",
    "f_1_other_ribozyme_drtr_bpp",
    "mcc_other_ribozyme_fold_mfe",
    "mcc_other_ribozyme_fold_bpp",
    "mcc_other_ribozyme_drtr_mfe",
    "mcc_other_ribozyme_drtr_prob",
    "mcc_other_ribozyme_drtr_bpp",
    "tpr_other_RNA_fold_mfe",
    "tpr_other_RNA_fold_bpp",
    "tpr_other_RNA_drtr_mfe",
    "tpr_other_RNA_drtr_prob",
    "tpr_other_RNA_drtr_bpp",
    "ppv_other_RNA_fold_mfe",
    "ppv_other_RNA_fold_bpp",
    "ppv_other_RNA_drtr_mfe",
    "ppv_other_RNA_drtr_prob",
    "ppv_other_RNA_drtr_bpp",
    "f_1_other_RNA_fold_mfe",
    "f_1_other_RNA_fold_bpp",
    "f_1_other_RNA_drtr_mfe",
    "f_1_other_RNA_drtr_prob",
    "f_1_other_RNA_drtr_bpp",
    "mcc_other_RNA_fold_mfe",
    "mcc_other_RNA_fold_bpp",
    "mcc_other_RNA_drtr_mfe",
    "mcc_other_RNA_drtr_prob",
    "mcc_other_RNA_drtr_bpp",
    "tpr_other_rRNA_fold_mfe",
    "tpr_other_rRNA_fold_bpp",
    "tpr_other_rRNA_drtr_mfe",
    "tpr_other_rRNA_drtr_prob",
    "tpr_other_rRNA_drtr_bpp",
    "ppv_other_rRNA_fold_mfe",
    "ppv_other_rRNA_fold_bpp",
    "ppv_other_rRNA_drtr_mfe",
    "ppv_other_rRNA_drtr_prob",
    "ppv_other_rRNA_drtr_bpp",
    "f_1_other_rRNA_fold_mfe",
    "f_1_other_rRNA_fold_bpp",
    "f_1_other_rRNA_drtr_mfe",
    "f_1_other_rRNA_drtr_prob",
    "f_1_other_rRNA_drtr_bpp",
    "mcc_other_rRNA_fold_mfe",
    "mcc_other_rRNA_fold_bpp",
    "mcc_other_rRNA_drtr_mfe",
    "mcc_other_rRNA_drtr_prob",
    "mcc_other_rRNA_drtr_bpp",
    "tpr_RNAIII_fold_mfe",
    "tpr_RNAIII_fold_bpp",
    "tpr_RNAIII_drtr_mfe",
    "tpr_RNAIII_drtr_prob",
    "tpr_RNAIII_drtr_bpp",
    "ppv_RNAIII_fold_mfe",
    "ppv_RNAIII_fold_bpp",
    "ppv_RNAIII_drtr_mfe",
    "ppv_RNAIII_drtr_prob",
    "ppv_RNAIII_drtr_bpp",
    "f_1_RNAIII_fold_mfe",
    "f_1_RNAIII_fold_bpp",
    "f_1_RNAIII_drtr_mfe",
    "f_1_RNAIII_drtr_prob",
    "f_1_RNAIII_drtr_bpp",
    "mcc_RNAIII_fold_mfe",
    "mcc_RNAIII_fold_bpp",
    "mcc_RNAIII_drtr_mfe",
    "mcc_RNAIII_drtr_prob",
    "mcc_RNAIII_drtr_bpp",
    "tpr_RNaseE5UTR_fold_mfe",
    "tpr_RNaseE5UTR_fold_bpp",
    "tpr_RNaseE5UTR_drtr_mfe",
    "tpr_RNaseE5UTR_drtr_prob",
    "tpr_RNaseE5UTR_drtr_bpp",
    "ppv_RNaseE5UTR_fold_mfe",
    "ppv_RNaseE5UTR_fold_bpp",
    "ppv_RNaseE5UTR_drtr_mfe",
    "ppv_RNaseE5UTR_drtr_prob",
    "ppv_RNaseE5UTR_drtr_bpp",
    "f_1_RNaseE5UTR_fold_mfe",
    "f_1_RNaseE5UTR_fold_bpp",
    "f_1_RNaseE5UTR_drtr_mfe",
    "f_1_RNaseE5UTR_drtr_prob",
    "f_1_RNaseE5UTR_drtr_bpp",
    "mcc_RNaseE5UTR_fold_mfe",
    "mcc_RNaseE5UTR_fold_bpp",
    "mcc_RNaseE5UTR_drtr_mfe",
    "mcc_RNaseE5UTR_drtr_prob",
    "mcc_RNaseE5UTR_drtr_bpp",
    "tpr_RNaseMRPRNA_fold_mfe",
    "tpr_RNaseMRPRNA_fold_bpp",
    "tpr_RNaseMRPRNA_drtr_mfe",
    "tpr_RNaseMRPRNA_drtr_prob",
    "tpr_RNaseMRPRNA_drtr_bpp",
    "ppv_RNaseMRPRNA_fold_mfe",
    "ppv_RNaseMRPRNA_fold_bpp",
    "ppv_RNaseMRPRNA_drtr_mfe",
    "ppv_RNaseMRPRNA_drtr_prob",
    "ppv_RNaseMRPRNA_drtr_bpp",
    "f_1_RNaseMRPRNA_fold_mfe",
    "f_1_RNaseMRPRNA_fold_bpp",
    "f_1_RNaseMRPRNA_drtr_mfe",
    "f_1_RNaseMRPRNA_drtr_prob",
    "f_1_RNaseMRPRNA_drtr_bpp",
    "mcc_RNaseMRPRNA_fold_mfe",
    "mcc_RNaseMRPRNA_fold_bpp",
    "mcc_RNaseMRPRNA_drtr_mfe",
    "mcc_RNaseMRPRNA_drtr_prob",
    "mcc_RNaseMRPRNA_drtr_bpp",
    "tpr_RNasePRNA_fold_mfe",
    "tpr_RNasePRNA_fold_bpp",
    "tpr_RNasePRNA_drtr_mfe",
    "tpr_RNasePRNA_drtr_prob",
    "tpr_RNasePRNA_drtr_bpp",
    "ppv_RNasePRNA_fold_mfe",
    "ppv_RNasePRNA_fold_bpp",
    "ppv_RNasePRNA_drtr_mfe",
    "ppv_RNasePRNA_drtr_prob",
    "ppv_RNasePRNA_drtr_bpp",
    "f_1_RNasePRNA_fold_mfe",
    "f_1_RNasePRNA_fold_bpp",
    "f_1_RNasePRNA_drtr_mfe",
    "f_1_RNasePRNA_drtr_prob",
    "f_1_RNasePRNA_drtr_bpp",
    "mcc_RNasePRNA_fold_mfe",
    "mcc_RNasePRNA_fold_bpp",
    "mcc_RNasePRNA_drtr_mfe",
    "mcc_RNasePRNA_drtr_prob",
    "mcc_RNasePRNA_drtr_bpp",
    "tpr_snRNA_fold_mfe",
    "tpr_snRNA_fold_bpp",
    "tpr_snRNA_drtr_mfe",
    "tpr_snRNA_drtr_prob",
    "tpr_snRNA_drtr_bpp",
    "ppv_snRNA_fold_mfe",
    "ppv_snRNA_fold_bpp",
    "ppv_snRNA_drtr_mfe",
    "ppv_snRNA_drtr_prob",
    "ppv_snRNA_drtr_bpp",
    "f_1_snRNA_fold_mfe",
    "f_1_snRNA_fold_bpp",
    "f_1_snRNA_drtr_mfe",
    "f_1_snRNA_drtr_prob",
    "f_1_snRNA_drtr_bpp",
    "mcc_snRNA_fold_mfe",
    "mcc_snRNA_fold_bpp",
    "mcc_snRNA_drtr_mfe",
    "mcc_snRNA_drtr_prob",
    "mcc_snRNA_drtr_bpp",
    "tpr_SRPRNA_fold_mfe",
    "tpr_SRPRNA_fold_bpp",
    "tpr_SRPRNA_drtr_mfe",
    "tpr_SRPRNA_drtr_prob",
    "tpr_SRPRNA_drtr_bpp",
    "ppv_SRPRNA_fold_mfe",
    "ppv_SRPRNA_fold_bpp",
    "ppv_SRPRNA_drtr_mfe",
    "ppv_SRPRNA_drtr_prob",
    "ppv_SRPRNA_drtr_bpp",
    "f_1_SRPRNA_fold_mfe",
    "f_1_SRPRNA_fold_bpp",
    "f_1_SRPRNA_drtr_mfe",
    "f_1_SRPRNA_drtr_prob",
    "f_1_SRPRNA_drtr_bpp",
    "mcc_SRPRNA_fold_mfe",
    "mcc_SRPRNA_fold_bpp",
    "mcc_SRPRNA_drtr_mfe",
    "mcc_SRPRNA_drtr_prob",
    "mcc_SRPRNA_drtr_bpp",
    "tpr_SyntheticRNA_fold_mfe",
    "tpr_SyntheticRNA_fold_bpp",
    "tpr_SyntheticRNA_drtr_mfe",
    "tpr_SyntheticRNA_drtr_prob",
    "tpr_SyntheticRNA_drtr_bpp",
    "ppv_SyntheticRNA_fold_mfe",
    "ppv_SyntheticRNA_fold_bpp",
    "ppv_SyntheticRNA_drtr_mfe",
    "ppv_SyntheticRNA_drtr_prob",
    "ppv_SyntheticRNA_drtr_bpp",
    "f_1_SyntheticRNA_fold_mfe",
    "f_1_SyntheticRNA_fold_bpp",
    "f_1_SyntheticRNA_drtr_mfe",
    "f_1_SyntheticRNA_drtr_prob",
    "f_1_SyntheticRNA_drtr_bpp",
    "mcc_SyntheticRNA_fold_mfe",
    "mcc_SyntheticRNA_fold_bpp",
    "mcc_SyntheticRNA_drtr_mfe",
    "mcc_SyntheticRNA_drtr_prob",
    "mcc_SyntheticRNA_drtr_bpp",
    "tpr_tmRNA_fold_mfe",
    "tpr_tmRNA_fold_bpp",
    "tpr_tmRNA_drtr_mfe",
    "tpr_tmRNA_drtr_prob",
    "tpr_tmRNA_drtr_bpp",
    "ppv_tmRNA_fold_mfe",
    "ppv_tmRNA_fold_bpp",
    "ppv_tmRNA_drtr_mfe",
    "ppv_tmRNA_drtr_prob",
    "ppv_tmRNA_drtr_bpp",
    "f_1_tmRNA_fold_mfe",
    "f_1_tmRNA_fold_bpp",
    "f_1_tmRNA_drtr_mfe",
    "f_1_tmRNA_drtr_prob",
    "f_1_tmRNA_drtr_bpp",
    "mcc_tmRNA_fold_mfe",
    "mcc_tmRNA_fold_bpp",
    "mcc_tmRNA_drtr_mfe",
    "mcc_tmRNA_drtr_prob",
    "mcc_tmRNA_drtr_bpp",
    "tpr_tRNA_fold_mfe",
    "tpr_tRNA_fold_bpp",
    "tpr_tRNA_drtr_mfe",
    "tpr_tRNA_drtr_prob",
    "tpr_tRNA_drtr_bpp",
    "ppv_tRNA_fold_mfe",
    "ppv_tRNA_fold_bpp",
    "ppv_tRNA_drtr_mfe",
    "ppv_tRNA_drtr_prob",
    "ppv_tRNA_drtr_bpp",
    "f_1_tRNA_fold_mfe",
    "f_1_tRNA_fold_bpp",
    "f_1_tRNA_drtr_mfe",
    "f_1_tRNA_drtr_prob",
    "f_1_tRNA_drtr_bpp",
    "mcc_tRNA_fold_mfe",
    "mcc_tRNA_fold_bpp",
    "mcc_tRNA_drtr_mfe",
    "mcc_tRNA_drtr_prob",
    "mcc_tRNA_drtr_bpp",
    "tpr_Vert_Telo_RNA_fold_mfe",
    "tpr_Vert_Telo_RNA_fold_bpp",
    "tpr_Vert_Telo_RNA_drtr_mfe",
    "tpr_Vert_Telo_RNA_drtr_prob",
    "tpr_Vert_Telo_RNA_drtr_bpp",
    "ppv_Vert_Telo_RNA_fold_mfe",
    "ppv_Vert_Telo_RNA_fold_bpp",
    "ppv_Vert_Telo_RNA_drtr_mfe",
    "ppv_Vert_Telo_RNA_drtr_prob",
    "ppv_Vert_Telo_RNA_drtr_bpp",
    "f_1_Vert_Telo_RNA_fold_mfe",
    "f_1_Vert_Telo_RNA_fold_bpp",
    "f_1_Vert_Telo_RNA_drtr_mfe",
    "f_1_Vert_Telo_RNA_drtr_prob",
    "f_1_Vert_Telo_RNA_drtr_bpp",
    "mcc_Vert_Telo_RNA_fold_mfe",
    "mcc_Vert_Telo_RNA_fold_bpp",
    "mcc_Vert_Telo_RNA_drtr_mfe",
    "mcc_Vert_Telo_RNA_drtr_prob",
    "mcc_Vert_Telo_RNA_drtr_bpp",
    "tpr_Viral_Phage_fold_mfe",
    "tpr_Viral_Phage_fold_bpp",
    "tpr_Viral_Phage_drtr_mfe",
    "tpr_Viral_Phage_drtr_prob",
    "tpr_Viral_Phage_drtr_bpp",
    "ppv_Viral_Phage_fold_mfe",
    "ppv_Viral_Phage_fold_bpp",
    "ppv_Viral_Phage_drtr_mfe",
    "ppv_Viral_Phage_drtr_prob",
    "ppv_Viral_Phage_drtr_bpp",
    "f_1_Viral_Phage_fold_mfe",
    "f_1_Viral_Phage_fold_bpp",
    "f_1_Viral_Phage_drtr_mfe",
    "f_1_Viral_Phage_drtr_prob",
    "f_1_Viral_Phage_drtr_bpp",
    "mcc_Viral_Phage_fold_mfe",
    "mcc_Viral_Phage_fold_bpp",
    "mcc_Viral_Phage_drtr_mfe",
    "mcc_Viral_Phage_drtr_prob",
    "mcc_Viral_Phage_drtr_bpp",
    "tpr_YRNA_fold_mfe",
    "tpr_YRNA_fold_bpp",
    "tpr_YRNA_drtr_mfe",
    "tpr_YRNA_drtr_prob",
    "tpr_YRNA_drtr_bpp",
    "ppv_YRNA_fold_mfe",
    "ppv_YRNA_fold_bpp",
    "ppv_YRNA_drtr_mfe",
    "ppv_YRNA_drtr_prob",
    "ppv_YRNA_drtr_bpp",
    "f_1_YRNA_fold_mfe",
    "f_1_YRNA_fold_bpp",
    "f_1_YRNA_drtr_mfe",
    "f_1_YRNA_drtr_prob",
    "f_1_YRNA_drtr_bpp",
    "mcc_YRNA_fold_mfe",
    "mcc_YRNA_fold_bpp",
    "mcc_YRNA_drtr_mfe",
    "mcc_YRNA_drtr_prob",
    "mcc_YRNA_drtr_bpp"
    ]


    cols = list()
    for i in range(0,boots):
        cols.append(i)

    data_processed = pd.DataFrame(data_processed, columns = cols, index = indexer)
    print(data_processed)
    export_csv = data_processed.to_csv(output_name)
