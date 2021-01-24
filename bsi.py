#!python


### ======================================= IMPORT ====================================== ###
import numpy as np
import pandas as pd
import math
# import RNA
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
                    default="comp_m_n_mfe_prob.csv",
                    help='parse directory')

parser.add_argument('--output','-o', action='store_true',
                    default=False,
                    help='recalculate boot/csvs')

args=parser.parse_args()


if args.test == False:
    input_name = args.input

if args.test == True:
    input_name = "comp_m_cp_mfe_prob.csv"

output_name = input_name[:-4] + "_s_boot.csv"


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
        f = 2 * (ppv * tpr)/(ppv + tpr)
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

if args.test==True:
    boots = 10
if args.test==False:
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


# seql = raw_data['sequencelength']
# drtr_tp_bpp = raw_data['drtr_tp']
# drtr_fp_bpp = raw_data['drtr_fp']
# drtr_tn_bpp = raw_data['drtr_tn']
# drtr_fn_bpp = raw_data['drtr_fn']
# drtr_tp_mfe = raw_data['drtr_tp_mfe']
# drtr_fp_mfe = raw_data['drtr_fp_mfe']
# drtr_tn_mfe = raw_data['drtr_tn_mfe']
# drtr_fn_mfe = raw_data['drtr_fn_mfe']
# drtr_tp_prob = raw_data['drtr_tp_prob']
# drtr_fp_prob = raw_data['drtr_fp_prob']
# drtr_tn_prob = raw_data['drtr_tn_prob']
# drtr_fn_prob = raw_data['drtr_fn_prob']
# fold_tp_mfe = raw_data['RNAfold_mfe_tp']
# fold_fp_mfe = raw_data['RNAfold_mfe_fp']
# fold_tn_mfe = raw_data['RNAfold_mfe_tn']
# fold_fn_mfe = raw_data['RNAfold_mfe_fn']
# fold_tp_bpp = raw_data['RNAfold_bpp_tp']
# fold_fp_bpp = raw_data['RNAfold_bpp_fp']
# fold_tn_bpp = raw_data['RNAfold_bpp_tn']
# fold_fn_bpp = raw_data['RNAfold_bpp_fn']
# len_data = len(seql)


### ==================================== split by family ============================ ###
size_data = raw_data
size_data.sort_values(by='sequencelength', inplace=True, axis=0)


### =================== differenceplot
size200 = size_data.loc[size_data.sequencelength < 200]

size500 = size_data.loc[size_data.sequencelength < 500]
size500 = size500.loc[size500.sequencelength > 200]

size1000 = size_data.loc[size_data.sequencelength < 1000]
size1000 = size1000.loc[size1000.sequencelength > 500]

size1000_plus = size_data.loc[size_data.sequencelength > 1000]

### ======================== 200

tp_200_fold_mfe = size200[['RNAfold_mfe_tp']]
tp_200_fold_mfe_sum = boot(tp_200_fold_mfe)
tp_200_fold_mfe_sum = pd.Series(tp_200_fold_mfe_sum)
#
# print(tp_200_fold_mfe.iloc[0])
# print(tp_200_fold_mfe_sum)

fp_200_fold_mfe = size200[['RNAfold_mfe_fp']]
fp_200_fold_mfe_sum = boot(fp_200_fold_mfe)
#
# print(fp_200_fold_mfe.iloc[0])
# print(fp_200_fold_mfe_sum)


tn_200_fold_mfe = size200[['RNAfold_mfe_tn']]
tn_200_fold_mfe_sum = boot(tn_200_fold_mfe)
# print(tn_200_fold_mfe.iloc[0])
# print(tn_200_fold_mfe_sum)


fn_200_fold_mfe = size200[['RNAfold_mfe_fn']]
fn_200_fold_mfe_sum = boot(fn_200_fold_mfe)
# print(fn_200_fold_mfe.iloc[0])
# print(fn_200_fold_mfe_sum)


tp_200_fold_bpp = size200[['RNAfold_bpp_tp']]
tp_200_fold_bpp_sum = boot(tp_200_fold_bpp)

fp_200_fold_bpp = size200[['RNAfold_bpp_fp']]
fp_200_fold_bpp_sum = boot(fp_200_fold_bpp)

tn_200_fold_bpp = size200[['RNAfold_bpp_tn']]
tn_200_fold_bpp_sum = boot(tn_200_fold_bpp)

fn_200_fold_bpp = size200[['RNAfold_bpp_fn']]
fn_200_fold_bpp_sum = boot(fn_200_fold_bpp)


tp_200_drtr_prob = size200[['drtr_tp_prob']]
tp_200_drtr_prob_sum = boot(tp_200_drtr_prob)

fp_200_drtr_prob = size200[['drtr_fp_prob']]
fp_200_drtr_prob_sum = boot(fp_200_drtr_prob)

tn_200_drtr_prob = size200[['drtr_tn_prob']]
tn_200_drtr_prob_sum = boot(tn_200_drtr_prob)

fn_200_drtr_prob = size200[['drtr_fn_prob']]
fn_200_drtr_prob_sum = boot(fn_200_drtr_prob)


tp_200_drtr_mfe = size200[['drtr_tp_mfe']]
tp_200_drtr_mfe_sum = boot(tp_200_drtr_mfe)

fp_200_drtr_mfe = size200[['drtr_fp_mfe']]
fp_200_drtr_mfe_sum = boot(fp_200_drtr_mfe)

tn_200_drtr_mfe = size200[['drtr_tn_mfe']]
tn_200_drtr_mfe_sum = boot(tn_200_drtr_mfe)

fn_200_drtr_mfe = size200[['drtr_fn_mfe']]
fn_200_drtr_mfe_sum = boot(fn_200_drtr_mfe)


tp_200_drtr_bpp = size200[['drtr_tp']]
tp_200_drtr_bpp_sum = boot(tp_200_drtr_bpp)

fp_200_drtr_bpp = size200[['drtr_fp']]
fp_200_drtr_bpp_sum = boot(fp_200_drtr_bpp)

tn_200_drtr_bpp = size200[['drtr_tn']]
tn_200_drtr_bpp_sum = boot(tn_200_drtr_bpp)

fn_200_drtr_bpp = size200[['drtr_fn']]
fn_200_drtr_bpp_sum = boot(fn_200_drtr_bpp)



#TPR
tpr_200_fold_mfe = tpr(tp_200_fold_mfe_sum, fn_200_fold_mfe_sum)

tpr_200_fold_bpp = tpr(tp_200_fold_bpp_sum, fn_200_fold_bpp_sum)

tpr_200_drtr_mfe = tpr(tp_200_drtr_mfe_sum, fn_200_drtr_mfe_sum)

tpr_200_drtr_prob = tpr(tp_200_drtr_prob_sum, fn_200_drtr_prob_sum)

tpr_200_drtr_bpp = tpr(tp_200_drtr_bpp_sum, fn_200_drtr_bpp_sum)


#PPV
ppv_200_fold_mfe = ppv(tp_200_fold_mfe_sum, fp_200_fold_mfe_sum)

ppv_200_fold_bpp = ppv(tp_200_fold_bpp_sum, fp_200_fold_bpp_sum)

ppv_200_drtr_mfe = ppv(tp_200_drtr_mfe_sum, fp_200_drtr_mfe_sum)

ppv_200_drtr_prob = ppv(tp_200_drtr_prob_sum, fp_200_drtr_prob_sum)

ppv_200_drtr_bpp = ppv(tp_200_drtr_bpp_sum, fp_200_drtr_bpp_sum)


#F1
f_1_200_fold_mfe = f_1(tp_200_fold_mfe_sum, fp_200_fold_mfe_sum, fn_200_fold_mfe_sum)

f_1_200_fold_bpp = f_1(tp_200_fold_bpp_sum, fp_200_fold_bpp_sum, fn_200_fold_bpp_sum)

f_1_200_drtr_mfe = f_1(tp_200_drtr_mfe_sum, fp_200_drtr_mfe_sum, fn_200_drtr_mfe_sum)

f_1_200_drtr_prob = f_1(tp_200_drtr_prob_sum, fp_200_drtr_prob_sum, fn_200_drtr_prob_sum)

f_1_200_drtr_bpp = f_1(tp_200_drtr_bpp_sum, fp_200_drtr_bpp_sum, fn_200_drtr_bpp_sum)


#MCC
mcc_200_fold_mfe = mcc(tp_200_fold_mfe_sum, fp_200_fold_mfe_sum, tn_200_fold_mfe_sum, fn_200_fold_mfe_sum)
# print(mcc_200_fold_mfe)

mcc_200_fold_bpp = mcc(tp_200_fold_bpp_sum, fp_200_fold_bpp_sum, tn_200_fold_bpp_sum, fn_200_fold_bpp_sum)

mcc_200_drtr_mfe = mcc(tp_200_drtr_mfe_sum, fp_200_drtr_mfe_sum, tn_200_drtr_mfe_sum, fn_200_drtr_mfe_sum)

mcc_200_drtr_prob = mcc(tp_200_drtr_prob_sum, fp_200_drtr_prob_sum,  tn_200_drtr_prob_sum, fn_200_drtr_prob_sum)

mcc_200_drtr_bpp = mcc(tp_200_drtr_bpp_sum, fp_200_drtr_bpp_sum,  tn_200_drtr_bpp_sum, fn_200_drtr_bpp_sum)

### ======================== 500

tp_500_fold_mfe = size500[['RNAfold_mfe_tp']]
tp_500_fold_mfe_sum = boot(tp_500_fold_mfe)

fp_500_fold_mfe = size500[['RNAfold_mfe_fp']]
fp_500_fold_mfe_sum = boot(fp_500_fold_mfe)

tn_500_fold_mfe = size500[['RNAfold_mfe_tn']]
tn_500_fold_mfe_sum = boot(tn_500_fold_mfe)

fn_500_fold_mfe = size500[['RNAfold_mfe_fn']]
fn_500_fold_mfe_sum = boot(fn_500_fold_mfe)


tp_500_fold_bpp = size500[['RNAfold_bpp_tp']]
tp_500_fold_bpp_sum = boot(tp_500_fold_bpp)

fp_500_fold_bpp = size500[['RNAfold_bpp_fp']]
fp_500_fold_bpp_sum = boot(fp_500_fold_bpp)

tn_500_fold_bpp = size500[['RNAfold_bpp_tn']]
tn_500_fold_bpp_sum = boot(tn_500_fold_bpp)

fn_500_fold_bpp = size500[['RNAfold_bpp_fn']]
fn_500_fold_bpp_sum = boot(fn_500_fold_bpp)


tp_500_drtr_prob = size500[['drtr_tp_prob']]
tp_500_drtr_prob_sum = boot(tp_500_drtr_prob)

fp_500_drtr_prob = size500[['drtr_fp_prob']]
fp_500_drtr_prob_sum = boot(fp_500_drtr_prob)

tn_500_drtr_prob = size500[['drtr_tn_prob']]
tn_500_drtr_prob_sum = boot(tn_500_drtr_prob)

fn_500_drtr_prob = size500[['drtr_fn_prob']]
fn_500_drtr_prob_sum = boot(fn_500_drtr_prob)


tp_500_drtr_mfe = size500[['drtr_tp_mfe']]
tp_500_drtr_mfe_sum = boot(tp_500_drtr_mfe)

fp_500_drtr_mfe = size500[['drtr_fp_mfe']]
fp_500_drtr_mfe_sum = boot(fp_500_drtr_mfe)

tn_500_drtr_mfe = size500[['drtr_tn_mfe']]
tn_500_drtr_mfe_sum = boot(tn_500_drtr_mfe)

fn_500_drtr_mfe = size500[['drtr_fn_mfe']]
fn_500_drtr_mfe_sum = boot(fn_500_drtr_mfe)


tp_500_drtr_bpp = size500[['drtr_tp']]
tp_500_drtr_bpp_sum = boot(tp_500_drtr_bpp)

fp_500_drtr_bpp = size500[['drtr_fp']]
fp_500_drtr_bpp_sum = boot(fp_500_drtr_bpp)

tn_500_drtr_bpp = size500[['drtr_tn']]
tn_500_drtr_bpp_sum = boot(tn_500_drtr_bpp)

fn_500_drtr_bpp = size500[['drtr_fn']]
fn_500_drtr_bpp_sum = boot(fn_500_drtr_bpp)


#TPR
tpr_500_fold_mfe = tpr(tp_500_fold_mfe_sum, fn_500_fold_mfe_sum)

tpr_500_fold_bpp = tpr(tp_500_fold_bpp_sum, fn_500_fold_bpp_sum)

tpr_500_drtr_mfe = tpr(tp_500_drtr_mfe_sum, fn_500_drtr_mfe_sum)

tpr_500_drtr_prob = tpr(tp_500_drtr_prob_sum, fn_500_drtr_prob_sum)

tpr_500_drtr_bpp = tpr(tp_500_drtr_bpp_sum, fn_500_drtr_bpp_sum)


#PPV
ppv_500_fold_mfe = ppv(tp_500_fold_mfe_sum, fp_500_fold_mfe_sum)

ppv_500_fold_bpp = ppv(tp_500_fold_bpp_sum, fp_500_fold_bpp_sum)

ppv_500_drtr_mfe = ppv(tp_500_drtr_mfe_sum, fp_500_drtr_mfe_sum)

ppv_500_drtr_prob = ppv(tp_500_drtr_prob_sum, fp_500_drtr_prob_sum)

ppv_500_drtr_bpp = ppv(tp_500_drtr_bpp_sum, fp_500_drtr_bpp_sum)


#F1
f_1_500_fold_mfe = f_1(tp_500_fold_mfe_sum, fp_500_fold_mfe_sum, fn_500_fold_mfe_sum)

f_1_500_fold_bpp = f_1(tp_500_fold_bpp_sum, fp_500_fold_bpp_sum, fn_500_fold_bpp_sum)

f_1_500_drtr_mfe = f_1(tp_500_drtr_mfe_sum, fp_500_drtr_mfe_sum, fn_500_drtr_mfe_sum)

f_1_500_drtr_prob = f_1(tp_500_drtr_prob_sum, fp_500_drtr_prob_sum, fn_500_drtr_prob_sum)

f_1_500_drtr_bpp = f_1(tp_500_drtr_bpp_sum, fp_500_drtr_bpp_sum, fn_500_drtr_bpp_sum)


#MCC
mcc_500_fold_mfe = mcc(tp_500_fold_mfe_sum, fp_500_fold_mfe_sum, tn_500_fold_mfe_sum, fn_500_fold_mfe_sum)
# print(mcc_500_fold_mfe)

mcc_500_fold_bpp = mcc(tp_500_fold_bpp_sum, fp_500_fold_bpp_sum, tn_500_fold_bpp_sum, fn_500_fold_bpp_sum)

mcc_500_drtr_mfe = mcc(tp_500_drtr_mfe_sum, fp_500_drtr_mfe_sum, tn_500_drtr_mfe_sum, fn_500_drtr_mfe_sum)

mcc_500_drtr_prob = mcc(tp_500_drtr_prob_sum, fp_500_drtr_prob_sum,  tn_500_drtr_prob_sum, fn_500_drtr_prob_sum)

mcc_500_drtr_bpp = mcc(tp_500_drtr_bpp_sum, fp_500_drtr_bpp_sum,  tn_500_drtr_bpp_sum, fn_500_drtr_bpp_sum)

### ======================== 1000

tp_1000_fold_mfe = size1000[['RNAfold_mfe_tp']]
tp_1000_fold_mfe_sum = boot(tp_1000_fold_mfe)

fp_1000_fold_mfe = size1000[['RNAfold_mfe_fp']]
fp_1000_fold_mfe_sum = boot(fp_1000_fold_mfe)

tn_1000_fold_mfe = size1000[['RNAfold_mfe_tn']]
tn_1000_fold_mfe_sum = boot(tn_1000_fold_mfe)

fn_1000_fold_mfe = size1000[['RNAfold_mfe_fn']]
fn_1000_fold_mfe_sum = boot(fn_1000_fold_mfe)


tp_1000_fold_bpp = size1000[['RNAfold_bpp_tp']]
tp_1000_fold_bpp_sum = boot(tp_1000_fold_bpp)

fp_1000_fold_bpp = size1000[['RNAfold_bpp_fp']]
fp_1000_fold_bpp_sum = boot(fp_1000_fold_bpp)

tn_1000_fold_bpp = size1000[['RNAfold_bpp_tn']]
tn_1000_fold_bpp_sum = boot(tn_1000_fold_bpp)

fn_1000_fold_bpp = size1000[['RNAfold_bpp_fn']]
fn_1000_fold_bpp_sum = boot(fn_1000_fold_bpp)


tp_1000_drtr_prob = size1000[['drtr_tp_prob']]
tp_1000_drtr_prob_sum = boot(tp_1000_drtr_prob)

fp_1000_drtr_prob = size1000[['drtr_fp_prob']]
fp_1000_drtr_prob_sum = boot(fp_1000_drtr_prob)

tn_1000_drtr_prob = size1000[['drtr_tn_prob']]
tn_1000_drtr_prob_sum = boot(tn_1000_drtr_prob)

fn_1000_drtr_prob = size1000[['drtr_fn_prob']]
fn_1000_drtr_prob_sum = boot(fn_1000_drtr_prob)


tp_1000_drtr_mfe = size1000[['drtr_tp_mfe']]
tp_1000_drtr_mfe_sum = boot(tp_1000_drtr_mfe)

fp_1000_drtr_mfe = size1000[['drtr_fp_mfe']]
fp_1000_drtr_mfe_sum = boot(fp_1000_drtr_mfe)

tn_1000_drtr_mfe = size1000[['drtr_tn_mfe']]
tn_1000_drtr_mfe_sum = boot(tn_1000_drtr_mfe)

fn_1000_drtr_mfe = size1000[['drtr_fn_mfe']]
fn_1000_drtr_mfe_sum = boot(fn_1000_drtr_mfe)


tp_1000_drtr_bpp = size1000[['drtr_tp']]
tp_1000_drtr_bpp_sum = boot(tp_1000_drtr_bpp)

fp_1000_drtr_bpp = size1000[['drtr_fp']]
fp_1000_drtr_bpp_sum = boot(fp_1000_drtr_bpp)

tn_1000_drtr_bpp = size1000[['drtr_tn']]
tn_1000_drtr_bpp_sum = boot(tn_1000_drtr_bpp)

fn_1000_drtr_bpp = size1000[['drtr_fn']]
fn_1000_drtr_bpp_sum = boot(fn_1000_drtr_bpp)



#TPR
tpr_1000_fold_mfe = tpr(tp_1000_fold_mfe_sum, fn_1000_fold_mfe_sum)

tpr_1000_fold_bpp = tpr(tp_1000_fold_bpp_sum, fn_1000_fold_bpp_sum)

tpr_1000_drtr_mfe = tpr(tp_1000_drtr_mfe_sum, fn_1000_drtr_mfe_sum)

tpr_1000_drtr_prob = tpr(tp_1000_drtr_prob_sum, fn_1000_drtr_prob_sum)

tpr_1000_drtr_bpp = tpr(tp_1000_drtr_bpp_sum, fn_1000_drtr_bpp_sum)


#PPV
ppv_1000_fold_mfe = ppv(tp_1000_fold_mfe_sum, fp_1000_fold_mfe_sum)

ppv_1000_fold_bpp = ppv(tp_1000_fold_bpp_sum, fp_1000_fold_bpp_sum)

ppv_1000_drtr_mfe = ppv(tp_1000_drtr_mfe_sum, fp_1000_drtr_mfe_sum)

ppv_1000_drtr_prob = ppv(tp_1000_drtr_prob_sum, fp_1000_drtr_prob_sum)

ppv_1000_drtr_bpp = ppv(tp_1000_drtr_bpp_sum, fp_1000_drtr_bpp_sum)


#F1
f_1_1000_fold_mfe = f_1(tp_1000_fold_mfe_sum, fp_1000_fold_mfe_sum, fn_1000_fold_mfe_sum)

f_1_1000_fold_bpp = f_1(tp_1000_fold_bpp_sum, fp_1000_fold_bpp_sum, fn_1000_fold_bpp_sum)

f_1_1000_drtr_mfe = f_1(tp_1000_drtr_mfe_sum, fp_1000_drtr_mfe_sum, fn_1000_drtr_mfe_sum)

f_1_1000_drtr_prob = f_1(tp_1000_drtr_prob_sum, fp_1000_drtr_prob_sum, fn_1000_drtr_prob_sum)

f_1_1000_drtr_bpp = f_1(tp_1000_drtr_bpp_sum, fp_1000_drtr_bpp_sum, fn_1000_drtr_bpp_sum)


#MCC
mcc_1000_fold_mfe = mcc(tp_1000_fold_mfe_sum, fp_1000_fold_mfe_sum, tn_1000_fold_mfe_sum, fn_1000_fold_mfe_sum)

mcc_1000_fold_bpp = mcc(tp_1000_fold_bpp_sum, fp_1000_fold_bpp_sum, tn_1000_fold_bpp_sum, fn_1000_fold_bpp_sum)

mcc_1000_drtr_mfe = mcc(tp_1000_drtr_mfe_sum, fp_1000_drtr_mfe_sum, tn_1000_drtr_mfe_sum, fn_1000_drtr_mfe_sum)

mcc_1000_drtr_prob = mcc(tp_1000_drtr_prob_sum, fp_1000_drtr_prob_sum,  tn_1000_drtr_prob_sum, fn_1000_drtr_prob_sum)

mcc_1000_drtr_bpp = mcc(tp_1000_drtr_bpp_sum, fp_1000_drtr_bpp_sum,  tn_1000_drtr_bpp_sum, fn_1000_drtr_bpp_sum)

### ===================== 1000 plus


tp_1000_plus_fold_mfe = size1000_plus[['RNAfold_mfe_tp']]
tp_1000_plus_fold_mfe_sum = boot(tp_1000_plus_fold_mfe)

fp_1000_plus_fold_mfe = size1000_plus[['RNAfold_mfe_fp']]
fp_1000_plus_fold_mfe_sum = boot(fp_1000_plus_fold_mfe)

tn_1000_plus_fold_mfe = size1000_plus[['RNAfold_mfe_tn']]
tn_1000_plus_fold_mfe_sum = boot(tn_1000_plus_fold_mfe)

fn_1000_plus_fold_mfe = size1000_plus[['RNAfold_mfe_fn']]
fn_1000_plus_fold_mfe_sum = boot(fn_1000_plus_fold_mfe)


tp_1000_plus_fold_bpp = size1000_plus[['RNAfold_bpp_tp']]
tp_1000_plus_fold_bpp_sum = boot(tp_1000_plus_fold_bpp)

fp_1000_plus_fold_bpp = size1000_plus[['RNAfold_bpp_fp']]
fp_1000_plus_fold_bpp_sum = boot(fp_1000_plus_fold_bpp)

tn_1000_plus_fold_bpp = size1000_plus[['RNAfold_bpp_tn']]
tn_1000_plus_fold_bpp_sum = boot(tn_1000_plus_fold_bpp)

fn_1000_plus_fold_bpp = size1000_plus[['RNAfold_bpp_fn']]
fn_1000_plus_fold_bpp_sum = boot(fn_1000_plus_fold_bpp)


tp_1000_plus_drtr_prob = size1000_plus[['drtr_tp_prob']]
tp_1000_plus_drtr_prob_sum = boot(tp_1000_plus_drtr_prob)

fp_1000_plus_drtr_prob = size1000_plus[['drtr_fp_prob']]
fp_1000_plus_drtr_prob_sum = boot(fp_1000_plus_drtr_prob)

tn_1000_plus_drtr_prob = size1000_plus[['drtr_tn_prob']]
tn_1000_plus_drtr_prob_sum = boot(tn_1000_plus_drtr_prob)

fn_1000_plus_drtr_prob = size1000_plus[['drtr_fn_prob']]
fn_1000_plus_drtr_prob_sum = boot(fn_1000_plus_drtr_prob)


tp_1000_plus_drtr_mfe = size1000_plus[['drtr_tp_mfe']]
tp_1000_plus_drtr_mfe_sum = boot(tp_1000_plus_drtr_mfe)

fp_1000_plus_drtr_mfe = size1000_plus[['drtr_fp_mfe']]
fp_1000_plus_drtr_mfe_sum = boot(fp_1000_plus_drtr_mfe)

tn_1000_plus_drtr_mfe = size1000_plus[['drtr_tn_mfe']]
tn_1000_plus_drtr_mfe_sum = boot(tn_1000_plus_drtr_mfe)

fn_1000_plus_drtr_mfe = size1000_plus[['drtr_fn_mfe']]
fn_1000_plus_drtr_mfe_sum = boot(fn_1000_plus_drtr_mfe)


tp_1000_plus_drtr_bpp = size1000_plus[['drtr_tp']]
tp_1000_plus_drtr_bpp_sum = boot(tp_1000_plus_drtr_bpp)

fp_1000_plus_drtr_bpp = size1000_plus[['drtr_fp']]
fp_1000_plus_drtr_bpp_sum = boot(fp_1000_plus_drtr_bpp)

tn_1000_plus_drtr_bpp = size1000_plus[['drtr_tn']]
tn_1000_plus_drtr_bpp_sum = boot(tn_1000_plus_drtr_bpp)

fn_1000_plus_drtr_bpp = size1000_plus[['drtr_fn']]
fn_1000_plus_drtr_bpp_sum = boot(fn_1000_plus_drtr_bpp)



#TPR
tpr_1000_plus_fold_mfe = tpr(tp_1000_plus_fold_mfe_sum, fn_1000_plus_fold_mfe_sum)

tpr_1000_plus_fold_bpp = tpr(tp_1000_plus_fold_bpp_sum, fn_1000_plus_fold_bpp_sum)

tpr_1000_plus_drtr_mfe = tpr(tp_1000_plus_drtr_mfe_sum, fn_1000_plus_drtr_mfe_sum)

tpr_1000_plus_drtr_prob = tpr(tp_1000_plus_drtr_prob_sum, fn_1000_plus_drtr_prob_sum)

tpr_1000_plus_drtr_bpp = tpr(tp_1000_plus_drtr_bpp_sum, fn_1000_plus_drtr_bpp_sum)


#PPV
ppv_1000_plus_fold_mfe = ppv(tp_1000_plus_fold_mfe_sum, fp_1000_plus_fold_mfe_sum)

ppv_1000_plus_fold_bpp = ppv(tp_1000_plus_fold_bpp_sum, fp_1000_plus_fold_bpp_sum)

ppv_1000_plus_drtr_mfe = ppv(tp_1000_plus_drtr_mfe_sum, fp_1000_plus_drtr_mfe_sum)

ppv_1000_plus_drtr_prob = ppv(tp_1000_plus_drtr_prob_sum, fp_1000_plus_drtr_prob_sum)

ppv_1000_plus_drtr_bpp = ppv(tp_1000_plus_drtr_bpp_sum, fp_1000_plus_drtr_bpp_sum)


#F1
f_1_1000_plus_fold_mfe = f_1(tp_1000_plus_fold_mfe_sum, fp_1000_plus_fold_mfe_sum, fn_1000_plus_fold_mfe_sum)

f_1_1000_plus_fold_bpp = f_1(tp_1000_plus_fold_bpp_sum, fp_1000_plus_fold_bpp_sum, fn_1000_plus_fold_bpp_sum)

f_1_1000_plus_drtr_mfe = f_1(tp_1000_plus_drtr_mfe_sum, fp_1000_plus_drtr_mfe_sum, fn_1000_plus_drtr_mfe_sum)

f_1_1000_plus_drtr_prob = f_1(tp_1000_plus_drtr_prob_sum, fp_1000_plus_drtr_prob_sum, fn_1000_plus_drtr_prob_sum)

f_1_1000_plus_drtr_bpp = f_1(tp_1000_plus_drtr_bpp_sum, fp_1000_plus_drtr_bpp_sum, fn_1000_plus_drtr_bpp_sum)


#MCC
mcc_1000_plus_fold_mfe = mcc(tp_1000_plus_fold_mfe_sum, fp_1000_plus_fold_mfe_sum, tn_1000_plus_fold_mfe_sum, fn_1000_plus_fold_mfe_sum)

mcc_1000_plus_fold_bpp = mcc(tp_1000_plus_fold_bpp_sum, fp_1000_plus_fold_bpp_sum, tn_1000_plus_fold_bpp_sum, fn_1000_plus_fold_bpp_sum)

mcc_1000_plus_drtr_mfe = mcc(tp_1000_plus_drtr_mfe_sum, fp_1000_plus_drtr_mfe_sum, tn_1000_plus_drtr_mfe_sum, fn_1000_plus_drtr_mfe_sum)

mcc_1000_plus_drtr_prob = mcc(tp_1000_plus_drtr_prob_sum, fp_1000_plus_drtr_prob_sum,  tn_1000_plus_drtr_prob_sum, fn_1000_plus_drtr_prob_sum)

mcc_1000_plus_drtr_bpp = mcc(tp_1000_plus_drtr_bpp_sum, fp_1000_plus_drtr_bpp_sum,  tn_1000_plus_drtr_bpp_sum, fn_1000_plus_drtr_bpp_sum)



if args.output == True:

    data_processed = [
    tpr_200_fold_mfe,
    tpr_200_fold_bpp,
    tpr_200_drtr_mfe,
    tpr_200_drtr_prob,
    tpr_200_drtr_bpp,
    ppv_200_fold_mfe,
    ppv_200_fold_bpp,
    ppv_200_drtr_mfe,
    ppv_200_drtr_prob,
    ppv_200_drtr_bpp,
    f_1_200_fold_mfe,
    f_1_200_fold_bpp,
    f_1_200_drtr_mfe,
    f_1_200_drtr_prob,
    f_1_200_drtr_bpp,
    mcc_200_fold_mfe,
    mcc_200_fold_bpp,
    mcc_200_drtr_mfe,
    mcc_200_drtr_prob,
    mcc_200_drtr_bpp,
    tpr_500_fold_mfe,
    tpr_500_fold_bpp,
    tpr_500_drtr_mfe,
    tpr_500_drtr_prob,
    tpr_500_drtr_bpp,
    ppv_500_fold_mfe,
    ppv_500_fold_bpp,
    ppv_500_drtr_mfe,
    ppv_500_drtr_prob,
    ppv_500_drtr_bpp,
    f_1_500_fold_mfe,
    f_1_500_fold_bpp,
    f_1_500_drtr_mfe,
    f_1_500_drtr_prob,
    f_1_500_drtr_bpp,
    mcc_500_fold_mfe,
    mcc_500_fold_bpp,
    mcc_500_drtr_mfe,
    mcc_500_drtr_prob,
    mcc_500_drtr_bpp,
    tpr_1000_fold_mfe,
    tpr_1000_fold_bpp,
    tpr_1000_drtr_mfe,
    tpr_1000_drtr_prob,
    tpr_1000_drtr_bpp,
    ppv_1000_fold_mfe,
    ppv_1000_fold_bpp,
    ppv_1000_drtr_mfe,
    ppv_1000_drtr_prob,
    ppv_1000_drtr_bpp,
    f_1_1000_fold_mfe,
    f_1_1000_fold_bpp,
    f_1_1000_drtr_mfe,
    f_1_1000_drtr_prob,
    f_1_1000_drtr_bpp,
    mcc_1000_fold_mfe,
    mcc_1000_fold_bpp,
    mcc_1000_drtr_mfe,
    mcc_1000_drtr_prob,
    mcc_1000_drtr_bpp,
    tpr_1000_plus_fold_mfe,
    tpr_1000_plus_fold_bpp,
    tpr_1000_plus_drtr_mfe,
    tpr_1000_plus_drtr_prob,
    tpr_1000_plus_drtr_bpp,
    ppv_1000_plus_fold_mfe,
    ppv_1000_plus_fold_bpp,
    ppv_1000_plus_drtr_mfe,
    ppv_1000_plus_drtr_prob,
    ppv_1000_plus_drtr_bpp,
    f_1_1000_plus_fold_mfe,
    f_1_1000_plus_fold_bpp,
    f_1_1000_plus_drtr_mfe,
    f_1_1000_plus_drtr_prob,
    f_1_1000_plus_drtr_bpp,
    mcc_1000_plus_fold_mfe,
    mcc_1000_plus_fold_bpp,
    mcc_1000_plus_drtr_mfe,
    mcc_1000_plus_drtr_prob,
    mcc_1000_plus_drtr_bpp
    ]


    indexer = [
    "tpr_200_fold_mfe",
    "tpr_200_fold_bpp",
    "tpr_200_drtr_mfe",
    "tpr_200_drtr_prob",
    "tpr_200_drtr_bpp",
    "ppv_200_fold_mfe",
    "ppv_200_fold_bpp",
    "ppv_200_drtr_mfe",
    "ppv_200_drtr_prob",
    "ppv_200_drtr_bpp",
    "f_1_200_fold_mfe",
    "f_1_200_fold_bpp",
    "f_1_200_drtr_mfe",
    "f_1_200_drtr_prob",
    "f_1_200_drtr_bpp",
    "mcc_200_fold_mfe",
    "mcc_200_fold_bpp",
    "mcc_200_drtr_mfe",
    "mcc_200_drtr_prob",
    "mcc_200_drtr_bpp",
    "tpr_500_fold_mfe",
    "tpr_500_fold_bpp",
    "tpr_500_drtr_mfe",
    "tpr_500_drtr_prob",
    "tpr_500_drtr_bpp",
    "ppv_500_fold_mfe",
    "ppv_500_fold_bpp",
    "ppv_500_drtr_mfe",
    "ppv_500_drtr_prob",
    "ppv_500_drtr_bpp",
    "f_1_500_fold_mfe",
    "f_1_500_fold_bpp",
    "f_1_500_drtr_mfe",
    "f_1_500_drtr_prob",
    "f_1_500_drtr_bpp",
    "mcc_500_fold_mfe",
    "mcc_500_fold_bpp",
    "mcc_500_drtr_mfe",
    "mcc_500_drtr_prob",
    "mcc_500_drtr_bpp",
    "tpr_1000_fold_mfe",
    "tpr_1000_fold_bpp",
    "tpr_1000_drtr_mfe",
    "tpr_1000_drtr_prob",
    "tpr_1000_drtr_bpp",
    "ppv_1000_fold_mfe",
    "ppv_1000_fold_bpp",
    "ppv_1000_drtr_mfe",
    "ppv_1000_drtr_prob",
    "ppv_1000_drtr_bpp",
    "f_1_1000_fold_mfe",
    "f_1_1000_fold_bpp",
    "f_1_1000_drtr_mfe",
    "f_1_1000_drtr_prob",
    "f_1_1000_drtr_bpp",
    "mcc_1000_fold_mfe",
    "mcc_1000_fold_bpp",
    "mcc_1000_drtr_mfe",
    "mcc_1000_drtr_prob",
    "mcc_1000_drtr_bpp",
    "tpr_1000_plus_fold_mfe",
    "tpr_1000_plus_fold_bpp",
    "tpr_1000_plus_drtr_mfe",
    "tpr_1000_plus_drtr_prob",
    "tpr_1000_plus_drtr_bpp",
    "ppv_1000_plus_fold_mfe",
    "ppv_1000_plus_fold_bpp",
    "ppv_1000_plus_drtr_mfe",
    "ppv_1000_plus_drtr_prob",
    "ppv_1000_plus_drtr_bpp",
    "f_1_1000_plus_fold_mfe",
    "f_1_1000_plus_fold_bpp",
    "f_1_1000_plus_drtr_mfe",
    "f_1_1000_plus_drtr_prob",
    "f_1_1000_plus_drtr_bpp",
    "mcc_1000_plus_fold_mfe",
    "mcc_1000_plus_fold_bpp",
    "mcc_1000_plus_drtr_mfe",
    "mcc_1000_plus_drtr_prob",
    "mcc_1000_plus_drtr_bpp"
    ]


    cols = list()
    for i in range(0,boots):
        cols.append(i)

    data_processed = pd.DataFrame(data_processed, columns = cols, index = indexer)
    print(data_processed)
    export_csv = data_processed.to_csv(output_name)



# data_processed = [
# fold_tpr_mfe,
# fold_tpr_bpp,
# drtr_tpr_mfe,
# drtr_tpr_prob,
# drtr_tpr_bpp,
# fold_ppv_mfe,
# fold_ppv_bpp,
# drtr_ppv_mfe,
# drtr_ppv_prob,
# drtr_ppv_bpp,
# fold_f_1_mfe,
# fold_f_1_bpp,
# drtr_f_1_mfe,
# drtr_f_1_prob,
# drtr_f_1_bpp,
# fold_mcc_mfe,
# fold_mcc_bpp,
# drtr_mcc_mfe,
# drtr_mcc_prob,
# drtr_mcc_bpp
# ]



# ## bootstraping_section
#
#
#
# ### ==================================== graphics ============================================ ###
#
#
# fold_mcc_mfes = [
# mcc_200_fold_mfe,
# mcc_500_fold_mfe,
# mcc_1000_fold_mfe,
# mcc_1000_plus_fold_mfe
# ]
#
# fold_mcc_bpps = [
# mcc_200_fold_bpp,
# mcc_500_fold_bpp,
# mcc_1000_fold_bpp,
# mcc_1000_plus_fold_bpp
# ]
#
# mcc_drtrs = [
# mcc_200_drtr,
# mcc_500_drtr,
# mcc_1000_drtr,
# mcc_1000_plus_drtr
# ]
#
# mcc_drtr_mfes = [
# mcc_200_drtr_mfe,
# mcc_500_drtr_mfe,
# mcc_1000_drtr_mfe,
# mcc_1000_plus_drtr_mfe
# ]
#
# # print(fold_mcc_bpps)
#
#
# if args.graphics == True:
#     titel = args.label
#     barwidth = 0.22
#
#     r_1 = np.arange(len(fold_mcc_mfes))
#     r_2 = [x + barwidth for x in r_1]
#     r_3 = [x + barwidth for x in r_2]
#     r_4 = [x + barwidth for x in r_3]
#
#
#     plt.bar(r_1, fold_mcc_mfes, width=barwidth, color='#247E85', label='RNAfold: mfe')
#     plt.bar(r_2, mcc_drtr_mfes, width=barwidth, color='#BE3151', label='DrTransformer: mfe')
#     plt.bar(r_3, fold_mcc_bpps, width=barwidth, color='#FA5882', label='RNAfold: bpp')
#     plt.bar(r_4, mcc_drtrs, width=barwidth, color='#0E3151', label='DrTransformer: bpp')
#
#     plt.title(titel, fontsize=20, y=1.03)
#     plt.xlabel('Length of sequence', fontsize=12)
#     plt.ylabel('MCC',  fontsize=16)
#     plt.axis([-0.2, 3.5, 0, 1])
#     plt.xticks([r + barwidth for r in range (len(fold_mcc_mfes))],
#                                                    ["200",
#                                                    "500",
#                                                    "1000",
#                                                    "1000+"],
#
#                                                    rotation=55,
#                                                    ha="right",
#                                                    rotation_mode="anchor",
#                                                    fontsize=12)
#
#     plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
#     plt.gcf().subplots_adjust(bottom=0.22, left=0.13, right=0.91, top=0.91)
#     # plt.grid()
#     plt.legend()
#
#     plt.show()
