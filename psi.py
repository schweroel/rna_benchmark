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
                    help='print plots')

parser.add_argument('--printout','-p', action='store_true',
                    default=False,
                    help='print all values per family')

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
                    help='produce csv')


# parser.add_argument('--output','-o', action='store',
#                     default="c:/Users/johan/Documents/VIENNARNAPRAKT/Johannes/testdir/finaltest/results_final.csv",
#                     help='parse directory')

args=parser.parse_args()


if args.test == False:
    input_name = args.input

if args.test == True:
    input_name = "comp_m_cp_fin.csv"

output_name = input_name[:-4] + "_s_mean.csv"
### =====================================-FUNCTIONS-==================================== ###

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

### ====================================== IMPORT =================================== ###
raw_data = pd.read_csv( open(input_name),float_precision = "high")
data_length = len(raw_data)


seql = raw_data['sequencelength']
pps_fold = raw_data['RNAfold_pps']
pps_drtr = raw_data['drtr_pps']
drtr_tp_mfe = raw_data['drtr_tp_mfe']
drtr_fp_mfe = raw_data['drtr_fp_mfe']
drtr_tn_mfe = raw_data['drtr_tn_mfe']
drtr_fn_mfe = raw_data['drtr_fn_mfe']
drtr_tp_prob = raw_data['drtr_tp_prob']
drtr_fp_prob = raw_data['drtr_fp_prob']
drtr_tn_prob = raw_data['drtr_tn_prob']
drtr_fn_prob = raw_data['drtr_fn_prob']
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



#
# pps_avg_fold = average(pps_fold)
# pps_avg_drtr = average(pps_drtr)
#
# ppv_fold_mfe = ppv(RNAfold_mfe_tp,RNAfold_mfe_fp)
# ppv_fold_bpp = ppv(RNAfold_bpp_tp,RNAfold_bpp_fp)
# ppv_drtr = ppv(drtr_tp,drtr_fp)
#
# tpr_fold_mfe = tpr(RNAfold_mfe_tp,RNAfold_mfe_fn)
# tpr_fold_bpp = tpr(RNAfold_bpp_tp,RNAfold_bpp_fn)
# tpr_drtr = tpr(drtr_tp,drtr_fn)
#
# f_measure_fold_mfe = f_measure(RNAfold_mfe_tp,RNAfold_mfe_fp,RNAfold_mfe_fn)
# f_measure_fold_bpp = f_measure(RNAfold_bpp_tp, RNAfold_bpp_fp, RNAfold_bpp_fn)
# f_measure_drtr = f_measure(drtr_tp, drtr_fp, drtr_fn)

# mcc_fold_mfe = mcc(RNAfold_mfe_tp, RNAfold_mfe_fp, RNAfold_mfe_tn, RNAfold_mfe_fn)
# mcc_fold_bpp = mcc(RNAfold_bpp_tp, RNAfold_bpp_fp, RNAfold_bpp_tn, RNAfold_bpp_fn)
# mcc_drtr = mcc(drtr_tp, drtr_fp, drtr_tn, drtr_fn)


#
# print("%8.5f,"
# "%8.5f,"
# "%8.5f,"
# "%8.5f,"
# "%8.8f,"
# "%8.8f,"
# "%8.8f,"
# "%8.8f,"
# "%8.8f,"
# "%8.8f,"
# "%8.8f,"
# "%8.8f,"
# "%8.8f,"
# "%8.8f" % (pps_avg_fold,
# pps_avg_drtr,
# ppv_fold_mfe,
# ppv_fold_bpp,
# ppv_drtr,
# tpr_fold_mfe,
# tpr_fold_bpp,
# tpr_drtr,
# f_measure_fold_mfe,
# f_measure_fold_bpp,
# f_measure_drtr,
# mcc_fold_mfe,
# mcc_fold_bpp,
# mcc_drtr))


### ==================================== split by family ============================ ###
size_data = raw_data
size_data.sort_values(by='sequencelength', inplace=True, axis=0)
# size_data.set_index(keys=['sequencelength'], drop=False, inplace=True)
### =================== differenceplot

size200= size_data.loc[size_data.sequencelength < 200]
# print(size200)

size500 = size_data.loc[size_data.sequencelength < 500]
size500 = size500.loc[size500.sequencelength > 200]
# print(size500)

size1000 = size_data.loc[size_data.sequencelength < 1000]
size1000 = size1000.loc[size1000.sequencelength > 500]
# print(size1000)

size1000_plus = size_data.loc[size_data.sequencelength > 1000]
# print(size1000_plus)

### ======================== 200

tp_200_fold_mfe = size200[['RNAfold_mfe_tp']]
fp_200_fold_mfe = size200[['RNAfold_mfe_fp']]
tn_200_fold_mfe = size200[['RNAfold_mfe_tn']]
fn_200_fold_mfe = size200[['RNAfold_mfe_fn']]

tp_200_fold_bpp = size200[['RNAfold_bpp_tp']]
fp_200_fold_bpp = size200[['RNAfold_bpp_fp']]
tn_200_fold_bpp = size200[['RNAfold_bpp_tn']]
fn_200_fold_bpp = size200[['RNAfold_bpp_fn']]

tp_200_drtr = size200[['drtr_tp']]
fp_200_drtr = size200[['drtr_fp']]
tn_200_drtr = size200[['drtr_tn']]
fn_200_drtr = size200[['drtr_fn']]

tp_200_drtr_mfe = size200[['drtr_tp_mfe']]
fp_200_drtr_mfe = size200[['drtr_fp_mfe']]
tn_200_drtr_mfe = size200[['drtr_tn_mfe']]
fn_200_drtr_mfe = size200[['drtr_fn_mfe']]

tp_200_drtr_prob = size200[['drtr_tp_prob']]
fp_200_drtr_prob = size200[['drtr_fp_prob']]
tn_200_drtr_prob = size200[['drtr_tn_prob']]
fn_200_drtr_prob = size200[['drtr_fn_prob']]


#TPR
tpr_200_fold_mfe = tpr(tp_200_fold_mfe, fn_200_fold_mfe)
# print(tpr_200_fold_mfe)

tpr_200_fold_bpp = tpr(tp_200_fold_bpp, fn_200_fold_bpp)
#print(tpr_200_fold_bpp)

tpr_200_drtr = tpr(tp_200_drtr, fn_200_drtr)
#print(tpr_200_drtr)

tpr_200_drtr_mfe = tpr(tp_200_drtr_mfe, fn_200_drtr_mfe)

tpr_200_drtr_prob = tpr(tp_200_drtr_prob, fn_200_drtr_prob)


#PPV
ppv_200_fold_mfe = ppv(tp_200_fold_mfe, fp_200_fold_mfe)
#print(ppv_200_fold_mfe)

ppv_200_fold_bpp = ppv(tp_200_fold_bpp, fp_200_fold_bpp)
#print(ppv_200_fold_bpp)

ppv_200_drtr = ppv(tp_200_drtr, fp_200_drtr)
#print(ppv_200_drtr)

ppv_200_drtr_mfe = ppv(tp_200_drtr_mfe, fp_200_drtr_mfe)

ppv_200_drtr_prob = ppv(tp_200_drtr_prob, fp_200_drtr_prob)


#F1
f1_200_fold_mfe = f_measure(tp_200_fold_mfe, fp_200_fold_mfe, fn_200_fold_mfe)
# print(f1_200_fold_mfe)

f1_200_fold_bpp = f_measure(tp_200_fold_bpp, fp_200_fold_bpp, fn_200_fold_bpp)
#print(f1_200_fold_bpp)

f1_200_drtr = f_measure(tp_200_drtr, fp_200_drtr, fn_200_drtr)
#print(f1_200_drtr)

f1_200_drtr_mfe = f_measure(tp_200_drtr_mfe, fp_200_drtr_mfe, fn_200_drtr_mfe)

f1_200_drtr_prob = f_measure(tp_200_drtr_prob, fp_200_drtr_prob, fn_200_drtr_prob)


#MCC
mcc_200_fold_mfe = mcc(tp_200_fold_mfe, fp_200_fold_mfe, tn_200_fold_mfe, fn_200_fold_mfe)
# print(mcc_200_fold_mfe)

mcc_200_fold_bpp = mcc(tp_200_fold_bpp, fp_200_fold_bpp, tn_200_fold_bpp, fn_200_fold_bpp)
#print(mcc_200_fold_bpp)

mcc_200_drtr = mcc(tp_200_drtr, fp_200_drtr, tn_200_drtr, fn_200_drtr)
#print(mcc_200_drtr)

mcc_200_drtr_mfe = mcc(tp_200_drtr_mfe, fp_200_drtr_mfe, tn_200_drtr_mfe, fn_200_drtr_mfe)

mcc_200_drtr_prob = mcc(tp_200_drtr_prob, fp_200_drtr_prob, tn_200_drtr_prob, fn_200_drtr_prob)


### ======================== 500

tp_500_fold_mfe = size500[['RNAfold_mfe_tp']]
fp_500_fold_mfe = size500[['RNAfold_mfe_fp']]
tn_500_fold_mfe = size500[['RNAfold_mfe_tn']]
fn_500_fold_mfe = size500[['RNAfold_mfe_fn']]

tp_500_fold_bpp = size500[['RNAfold_bpp_tp']]
fp_500_fold_bpp = size500[['RNAfold_bpp_fp']]
tn_500_fold_bpp = size500[['RNAfold_bpp_tn']]
fn_500_fold_bpp = size500[['RNAfold_bpp_fn']]

tp_500_drtr = size500[['drtr_tp']]
fp_500_drtr = size500[['drtr_fp']]
tn_500_drtr = size500[['drtr_tn']]
fn_500_drtr = size500[['drtr_fn']]

tp_500_drtr_mfe = size500[['drtr_tp_mfe']]
fp_500_drtr_mfe = size500[['drtr_fp_mfe']]
tn_500_drtr_mfe = size500[['drtr_tn_mfe']]
fn_500_drtr_mfe = size500[['drtr_fn_mfe']]

tp_500_drtr_prob = size500[['drtr_tp_prob']]
fp_500_drtr_prob = size500[['drtr_fp_prob']]
tn_500_drtr_prob = size500[['drtr_tn_prob']]
fn_500_drtr_prob = size500[['drtr_fn_prob']]

#TPR
tpr_500_fold_mfe = tpr(tp_500_fold_mfe, fn_500_fold_mfe)
#print(tpr_500_fold_mfe)

tpr_500_fold_bpp = tpr(tp_500_fold_bpp, fn_500_fold_bpp)
#print(tpr_500_fold_bpp)

tpr_500_drtr = tpr(tp_500_drtr, fn_500_drtr)
#print(tpr_500_drtr)

tpr_500_drtr_mfe = tpr(tp_500_drtr_mfe, fn_500_drtr_mfe)

tpr_500_drtr_prob = tpr(tp_500_drtr_prob, fn_500_drtr_prob)


#PPV
ppv_500_fold_mfe = ppv(tp_500_fold_mfe, fp_500_fold_mfe)
#print(ppv_500_fold_mfe)

ppv_500_fold_bpp = ppv(tp_500_fold_bpp, fp_500_fold_bpp)
#print(ppv_500_fold_bpp)

ppv_500_drtr = ppv(tp_500_drtr, fp_500_drtr)
#print(ppv_500_drtr)

ppv_500_drtr_mfe = ppv(tp_500_drtr_mfe, fp_500_drtr_mfe)

ppv_500_drtr_prob = ppv(tp_500_drtr_prob, fp_500_drtr_prob)


#F1
f1_500_fold_mfe = f_measure(tp_500_fold_mfe, fp_500_fold_mfe, fn_500_fold_mfe)
#print(f1_500_fold_mfe)

f1_500_fold_bpp = f_measure(tp_500_fold_bpp, fp_500_fold_bpp, fn_500_fold_bpp)
#print(f1_500_fold_bpp)

f1_500_drtr = f_measure(tp_500_drtr, fp_500_drtr, fn_500_drtr)
#print(f1_500_drtr)

f1_500_drtr_mfe = f_measure(tp_500_drtr_mfe, fp_500_drtr_mfe, fn_500_drtr_mfe)

f1_500_drtr_prob = f_measure(tp_500_drtr_prob, fp_500_drtr_prob, fn_500_drtr_prob)


#MCC
mcc_500_fold_mfe = mcc (tp_500_fold_mfe, fp_500_fold_mfe, tn_500_fold_mfe, fn_500_fold_mfe)
#print(mcc_500_fold_mfe)

mcc_500_fold_bpp = mcc (tp_500_fold_bpp, fp_500_fold_bpp, tn_500_fold_bpp, fn_500_fold_bpp)
#print(mcc_500_fold_bpp)

mcc_500_drtr = mcc (tp_500_drtr, fp_500_drtr, tn_500_drtr, fn_500_drtr)
#print(mcc_500_drtr)

mcc_500_drtr_mfe = mcc (tp_500_drtr_mfe, fp_500_drtr_mfe, tn_500_drtr_mfe, fn_500_drtr_mfe)

mcc_500_drtr_prob = mcc (tp_500_drtr_prob, fp_500_drtr_prob, tn_500_drtr_prob, fn_500_drtr_prob)


### ======================== 1000

tp_1000_fold_mfe = size1000[['RNAfold_mfe_tp']]
fp_1000_fold_mfe = size1000[['RNAfold_mfe_fp']]
tn_1000_fold_mfe = size1000[['RNAfold_mfe_tn']]
fn_1000_fold_mfe = size1000[['RNAfold_mfe_fn']]

tp_1000_fold_bpp = size1000[['RNAfold_bpp_tp']]
fp_1000_fold_bpp = size1000[['RNAfold_bpp_fp']]
tn_1000_fold_bpp = size1000[['RNAfold_bpp_tn']]
fn_1000_fold_bpp = size1000[['RNAfold_bpp_fn']]

tp_1000_drtr = size1000[['drtr_tp']]
fp_1000_drtr = size1000[['drtr_fp']]
tn_1000_drtr = size1000[['drtr_tn']]
fn_1000_drtr = size1000[['drtr_fn']]

tp_1000_drtr_mfe = size1000[['drtr_tp_mfe']]
fp_1000_drtr_mfe = size1000[['drtr_fp_mfe']]
tn_1000_drtr_mfe = size1000[['drtr_tn_mfe']]
fn_1000_drtr_mfe = size1000[['drtr_fn_mfe']]

tp_1000_drtr_prob = size1000[['drtr_tp_prob']]
fp_1000_drtr_prob = size1000[['drtr_fp_prob']]
tn_1000_drtr_prob = size1000[['drtr_tn_prob']]
fn_1000_drtr_prob = size1000[['drtr_fn_prob']]

#TPR
tpr_1000_fold_mfe = tpr(tp_1000_fold_mfe, fn_1000_fold_mfe)
#print(tpr_1000_fold_mfe)

tpr_1000_fold_bpp = tpr(tp_1000_fold_bpp, fn_1000_fold_bpp)
#print(tpr_1000_fold_bpp)

tpr_1000_drtr = tpr(tp_1000_drtr, fn_1000_drtr)
#print(tpr_1000_drtr)

tpr_1000_drtr_mfe = tpr(tp_1000_drtr_mfe, fn_1000_drtr_mfe)

tpr_1000_drtr_prob = tpr(tp_1000_drtr_prob, fn_1000_drtr_prob)



#PPV
ppv_1000_fold_mfe = ppv(tp_1000_fold_mfe, fp_1000_fold_mfe)
#print(ppv_1000_fold_mfe)

ppv_1000_fold_bpp = ppv(tp_1000_fold_bpp, fp_1000_fold_bpp)
#print(ppv_1000_fold_bpp)

ppv_1000_drtr = ppv(tp_1000_drtr, fp_1000_drtr)
#print(ppv_1000_drtr)

ppv_1000_drtr_mfe = ppv(tp_1000_drtr_mfe, fp_1000_drtr_mfe)

ppv_1000_drtr_prob = ppv(tp_1000_drtr_prob, fp_1000_drtr_prob)


#F1
f1_1000_fold_mfe = f_measure(tp_1000_fold_mfe, fp_1000_fold_mfe, fn_1000_fold_mfe)
#print(f1_1000_fold_mfe)

f1_1000_fold_bpp = f_measure(tp_1000_fold_bpp, fp_1000_fold_bpp, fn_1000_fold_bpp)
#print(f1_1000_fold_bpp)

f1_1000_drtr = f_measure(tp_1000_drtr, fp_1000_drtr, fn_1000_drtr)
#print(f1_1000_drtr)

f1_1000_drtr_mfe = f_measure(tp_1000_drtr_mfe, fp_1000_drtr_mfe, fn_1000_drtr_mfe)

f1_1000_drtr_prob = f_measure(tp_1000_drtr_prob, fp_1000_drtr_prob, fn_1000_drtr_prob)


#MCC
mcc_1000_fold_mfe = mcc (tp_1000_fold_mfe, fp_1000_fold_mfe, tn_1000_fold_mfe, fn_1000_fold_mfe)
#print(mcc_1000_fold_mfe)

mcc_1000_fold_bpp = mcc (tp_1000_fold_bpp, fp_1000_fold_bpp, tn_1000_fold_bpp, fn_1000_fold_bpp)
#print(mcc_1000_fold_bpp)

mcc_1000_drtr = mcc (tp_1000_drtr, fp_1000_drtr, tn_1000_drtr, fn_1000_drtr)
#print(mcc_1000_drtr)

mcc_1000_drtr_mfe = mcc (tp_1000_drtr_mfe, fp_1000_drtr_mfe, tn_1000_drtr_mfe, fn_1000_drtr_mfe)

mcc_1000_drtr_prob = mcc (tp_1000_drtr_prob, fp_1000_drtr_prob, tn_1000_drtr_prob, fn_1000_drtr_prob)


### ===================== 1000 plus

tp_1000_plus_fold_mfe = size1000_plus[['RNAfold_mfe_tp']]
fp_1000_plus_fold_mfe = size1000_plus[['RNAfold_mfe_fp']]
tn_1000_plus_fold_mfe = size1000_plus[['RNAfold_mfe_tn']]
fn_1000_plus_fold_mfe = size1000_plus[['RNAfold_mfe_fn']]

tp_1000_plus_fold_bpp = size1000_plus[['RNAfold_bpp_tp']]
fp_1000_plus_fold_bpp = size1000_plus[['RNAfold_bpp_fp']]
tn_1000_plus_fold_bpp = size1000_plus[['RNAfold_bpp_tn']]
fn_1000_plus_fold_bpp = size1000_plus[['RNAfold_bpp_fn']]

tp_1000_plus_drtr = size1000_plus[['drtr_tp']]
fp_1000_plus_drtr = size1000_plus[['drtr_fp']]
tn_1000_plus_drtr = size1000_plus[['drtr_tn']]
fn_1000_plus_drtr = size1000_plus[['drtr_fn']]

tp_1000_plus_drtr_mfe = size1000_plus[['drtr_tp_mfe']]
fp_1000_plus_drtr_mfe = size1000_plus[['drtr_fp_mfe']]
tn_1000_plus_drtr_mfe = size1000_plus[['drtr_tn_mfe']]
fn_1000_plus_drtr_mfe = size1000_plus[['drtr_fn_mfe']]

tp_1000_plus_drtr_prob = size1000_plus[['drtr_tp_prob']]
fp_1000_plus_drtr_prob = size1000_plus[['drtr_fp_prob']]
tn_1000_plus_drtr_prob = size1000_plus[['drtr_tn_prob']]
fn_1000_plus_drtr_prob = size1000_plus[['drtr_fn_prob']]



#TPR
tpr_1000_plus_fold_mfe = tpr(tp_1000_plus_fold_mfe, fn_1000_plus_fold_mfe)
#print(tpr_1000_plus_fold_mfe)

tpr_1000_plus_fold_bpp = tpr(tp_1000_plus_fold_bpp, fn_1000_plus_fold_bpp)
#print(tpr_1000_plus_fold_bpp)

tpr_1000_plus_drtr = tpr(tp_1000_plus_drtr, fn_1000_plus_drtr)
#print(tpr_1000_plus_drtr)

tpr_1000_plus_drtr_mfe = tpr(tp_1000_plus_drtr_mfe, fn_1000_plus_drtr_mfe)

tpr_1000_plus_drtr_prob = tpr(tp_1000_plus_drtr_prob, fn_1000_plus_drtr_prob)


#PPV
ppv_1000_plus_fold_mfe = ppv(tp_1000_plus_fold_mfe, fp_1000_plus_fold_mfe)
#print(ppv_1000_plus_fold_mfe)

ppv_1000_plus_fold_bpp = ppv(tp_1000_plus_fold_bpp, fp_1000_plus_fold_bpp)
#print(ppv_1000_plus_fold_bpp)

ppv_1000_plus_drtr = ppv(tp_1000_plus_drtr, fp_1000_plus_drtr)
#print(ppv_1000_plus_drtr)

ppv_1000_plus_drtr_mfe = ppv(tp_1000_plus_drtr_mfe, fp_1000_plus_drtr_mfe)

ppv_1000_plus_drtr_prob = ppv(tp_1000_plus_drtr_prob, fp_1000_plus_drtr_prob)


#F1
f1_1000_plus_fold_mfe = f_measure(tp_1000_plus_fold_mfe, fp_1000_plus_fold_mfe, fn_1000_plus_fold_mfe)
#print(f1_1000_plus_fold_mfe)

f1_1000_plus_fold_bpp = f_measure(tp_1000_plus_fold_bpp, fp_1000_plus_fold_bpp, fn_1000_plus_fold_bpp)
#print(f1_1000_plus_fold_bpp)

f1_1000_plus_drtr = f_measure(tp_1000_plus_drtr, fp_1000_plus_drtr, fn_1000_plus_drtr)
#print(f1_1000_plus_drtr)

f1_1000_plus_drtr_mfe = f_measure(tp_1000_plus_drtr_mfe, fp_1000_plus_drtr_mfe, fn_1000_plus_drtr_mfe)

f1_1000_plus_drtr_prob = f_measure(tp_1000_plus_drtr_prob, fp_1000_plus_drtr_prob, fn_1000_plus_drtr_prob)


#MCC
mcc_1000_plus_fold_mfe = mcc (tp_1000_plus_fold_mfe, fp_1000_plus_fold_mfe, tn_1000_plus_fold_mfe, fn_1000_plus_fold_mfe)
#print(mcc_1000_plus_fold_mfe)

mcc_1000_plus_fold_bpp = mcc (tp_1000_plus_fold_bpp, fp_1000_plus_fold_bpp, tn_1000_plus_fold_bpp, fn_1000_plus_fold_bpp)
#print(mcc_1000_plus_fold_bpp)

mcc_1000_plus_drtr = mcc (tp_1000_plus_drtr, fp_1000_plus_drtr, tn_1000_plus_drtr, fn_1000_plus_drtr)
#print(mcc_1000_plus_drtr)

mcc_1000_plus_drtr_mfe = mcc (tp_1000_plus_drtr_mfe, fp_1000_plus_drtr_mfe, tn_1000_plus_drtr_mfe, fn_1000_plus_drtr_mfe)

mcc_1000_plus_drtr_prob = mcc (tp_1000_plus_drtr_prob, fp_1000_plus_drtr_prob, tn_1000_plus_drtr_prob, fn_1000_plus_drtr_prob)



### bootstraping_section

### ==================================== printer ============================================= ###
if args.printout == True:
    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    ppv_200_fold_mfe,
    ppv_200_drtr_mfe,
    ppv_200_drtr_prob,
    ppv_200_fold_bpp,
    ppv_200_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    tpr_200_fold_mfe,
    tpr_200_drtr_mfe,
    tpr_200_drtr_prob,
    tpr_200_fold_bpp,
    tpr_200_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    f1_200_fold_mfe,
    f1_200_drtr_mfe,
    f1_200_drtr_prob,
    f1_200_fold_bpp,
    f1_200_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    mcc_200_fold_mfe,
    mcc_200_drtr_mfe,
    mcc_200_drtr_prob,
    mcc_200_fold_bpp,
    mcc_200_drtr
    ))

    # 500
    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    ppv_500_fold_mfe,
    ppv_500_drtr_mfe,
    ppv_500_drtr_prob,
    ppv_500_fold_bpp,
    ppv_500_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    tpr_500_fold_mfe,
    tpr_500_drtr_mfe,
    tpr_500_drtr_prob,
    tpr_500_fold_bpp,
    tpr_500_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    f1_500_fold_mfe,
    f1_500_drtr_mfe,
    f1_500_drtr_prob,
    f1_500_fold_bpp,
    f1_500_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    mcc_500_fold_mfe,
    mcc_500_drtr_mfe,
    mcc_500_drtr_prob,
    mcc_500_fold_bpp,
    mcc_500_drtr
    ))


    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    ppv_1000_fold_mfe,
    ppv_1000_drtr_mfe,
    ppv_1000_drtr_prob,
    ppv_1000_fold_bpp,
    ppv_1000_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    tpr_1000_fold_mfe,
    tpr_1000_drtr_mfe,
    tpr_1000_drtr_prob,
    tpr_1000_fold_bpp,
    tpr_1000_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    f1_1000_fold_mfe,
    f1_1000_drtr_mfe,
    f1_1000_drtr_prob,
    f1_1000_fold_bpp,
    f1_1000_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    mcc_1000_fold_mfe,
    mcc_1000_drtr_mfe,
    mcc_1000_drtr_prob,
    mcc_1000_fold_bpp,
    mcc_1000_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    ppv_1000_plus_fold_mfe,
    ppv_1000_plus_drtr_mfe,
    ppv_1000_plus_drtr_prob,
    ppv_1000_plus_fold_bpp,
    ppv_1000_plus_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    tpr_1000_plus_fold_mfe,
    tpr_1000_plus_drtr_mfe,
    tpr_1000_plus_drtr_prob,
    tpr_1000_plus_fold_bpp,
    tpr_1000_plus_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    f1_1000_plus_fold_mfe,
    f1_1000_plus_drtr_mfe,
    f1_1000_plus_drtr_prob,
    f1_1000_plus_fold_bpp,
    f1_1000_plus_drtr
    ))

    print(
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    %
    (
    mcc_1000_plus_fold_mfe,
    mcc_1000_plus_drtr_mfe,
    mcc_1000_plus_drtr_prob,
    mcc_1000_plus_fold_bpp,
    mcc_1000_plus_drtr
    ))

### ==================================== graphics ============================================ ###


fold_mcc_mfes = [
mcc_200_fold_mfe,
mcc_500_fold_mfe,
mcc_1000_fold_mfe,
mcc_1000_plus_fold_mfe
]

fold_mcc_bpps = [
mcc_200_fold_bpp,
mcc_500_fold_bpp,
mcc_1000_fold_bpp,
mcc_1000_plus_fold_bpp
]

mcc_drtrs = [
mcc_200_drtr,
mcc_500_drtr,
mcc_1000_drtr,
mcc_1000_plus_drtr
]

mcc_drtr_mfes = [
mcc_200_drtr_mfe,
mcc_500_drtr_mfe,
mcc_1000_drtr_mfe,
mcc_1000_plus_drtr_mfe
]

mcc_drtr_probs = [
mcc_200_drtr_prob,
mcc_500_drtr_prob,
mcc_1000_drtr_prob,
mcc_1000_plus_drtr_prob
]

# print(fold_mcc_bpps)

titel = args.label

if args.graphics == True:

    barwidth = 0.19

    r_1 = np.arange(len(fold_mcc_mfes))
    r_2 = [x + barwidth for x in r_1]
    r_3 = [x + barwidth for x in r_2]
    r_4 = [x + barwidth for x in r_3]
    r_5 = [x + barwidth for x in r_4]



    plt.bar(r_1, fold_mcc_mfes, width=barwidth, color='#247E85', label='RNAfold: mfe')
    plt.bar(r_2, mcc_drtr_mfes, width=barwidth, color='#BE3151', label='DrTransformer: mfe')
    plt.bar(r_3, mcc_drtr_probs, width=barwidth, color='#5A3151', label='DrTransformer: prob')
    plt.bar(r_4, fold_mcc_bpps, width=barwidth, color='#FA5882', label='RNAfold: bpp')
    plt.bar(r_5, mcc_drtrs, width=barwidth, color='#0E3151', label='DrTransformer: bpp')


    plt.title(titel, fontsize=20, y=1.03)
    plt.xlabel('Length of sequence', fontsize=12)
    plt.ylabel('MCC',  fontsize=16)
    plt.axis([-0.2, 3.95, 0, 1])
    plt.xticks([r + barwidth for r in range (len(fold_mcc_mfes))],
                                                   ["200 (n=" + str(len(size200)) + ")",
                                                   "500 (n=" + str(len(size500)) + ")",
                                                   "1000 (n=" + str(len(size1000)) + ")",
                                                   ">1000 (n=" + str(len(size1000_plus)) + ")"],

                                                   rotation=55,
                                                   ha="right",
                                                   rotation_mode="anchor",
                                                   fontsize=12)

    plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    plt.gcf().subplots_adjust(bottom=0.32, left=0.13, right=0.91, top=0.88)
    # plt.grid()
    plt.legend()

    plt.show()


if args.output == True:
    data_200 = [
    ppv_200_fold_mfe,
    ppv_200_fold_bpp,
    ppv_200_drtr_mfe,
    ppv_200_drtr_prob,
    ppv_200_drtr,
    tpr_200_fold_mfe,
    tpr_200_fold_bpp,
    tpr_200_drtr_mfe,
    tpr_200_drtr_prob,
    tpr_200_drtr,
    f1_200_fold_mfe,
    f1_200_fold_bpp,
    f1_200_drtr_mfe,
    f1_200_drtr_prob,
    f1_200_drtr,
    mcc_200_fold_mfe,
    mcc_200_fold_bpp,
    mcc_200_drtr_mfe,
    mcc_200_drtr_prob,
    mcc_200_drtr,
    ]
    data_500 = [
    ppv_500_fold_mfe,
    ppv_500_fold_bpp,
    ppv_500_drtr_mfe,
    ppv_500_drtr_prob,
    ppv_500_drtr,
    tpr_500_fold_mfe,
    tpr_500_fold_bpp,
    tpr_500_drtr_mfe,
    tpr_500_drtr_prob,
    tpr_500_drtr,
    f1_500_fold_mfe,
    f1_500_fold_bpp,
    f1_500_drtr_mfe,
    f1_500_drtr_prob,
    f1_500_drtr,
    mcc_500_fold_mfe,
    mcc_500_fold_bpp,
    mcc_500_drtr_mfe,
    mcc_500_drtr_prob,
    mcc_500_drtr,
    ]
    data_1000 = [
    ppv_1000_fold_mfe,
    ppv_1000_fold_bpp,
    ppv_1000_drtr_mfe,
    ppv_1000_drtr_prob,
    ppv_1000_drtr,
    tpr_1000_fold_mfe,
    tpr_1000_fold_bpp,
    tpr_1000_drtr_mfe,
    tpr_1000_drtr_prob,
    tpr_1000_drtr,
    f1_1000_fold_mfe,
    f1_1000_fold_bpp,
    f1_1000_drtr_mfe,
    f1_1000_drtr_prob,
    f1_1000_drtr,
    mcc_1000_fold_mfe,
    mcc_1000_fold_bpp,
    mcc_1000_drtr_mfe,
    mcc_1000_drtr_prob,
    mcc_1000_drtr,
    ]
    data_1000_plus = [
    ppv_1000_plus_fold_mfe,
    ppv_1000_plus_fold_bpp,
    ppv_1000_plus_drtr_mfe,
    ppv_1000_plus_drtr_prob,
    ppv_1000_plus_drtr,
    tpr_1000_plus_fold_mfe,
    tpr_1000_plus_fold_bpp,
    tpr_1000_plus_drtr_mfe,
    tpr_1000_plus_drtr_prob,
    tpr_1000_plus_drtr,
    f1_1000_plus_fold_mfe,
    f1_1000_plus_fold_bpp,
    f1_1000_plus_drtr_mfe,
    f1_1000_plus_drtr_prob,
    f1_1000_plus_drtr,
    mcc_1000_plus_fold_mfe,
    mcc_1000_plus_fold_bpp,
    mcc_1000_plus_drtr_mfe,
    mcc_1000_plus_drtr_prob,
    mcc_1000_plus_drtr,
    ]
    data_all = [data_200, data_500, data_1000, data_1000_plus]

    print(output_name)
    data_all = pd.DataFrame(data_all)
    print(data_all)
    export_csv = data_all.to_csv(output_name)
