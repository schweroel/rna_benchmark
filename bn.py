#!python

# content: python script that calculates a distribution based on bootstrapping
# for the whole data-set

### ================= import modules ======================================= ###
import numpy as np
import pandas as pd
from random import seed
from random import randint
import pandas as pd
import argparse
import math
import sys

### ================= settings for pandas and numpy ======================== ###
# prevents output from bein truncated when printed
pd.set_option('display.max_rows', 10000)
pd.set_option('display.max_columns', 10000)
pd.set_option('display.width', 10000)
np.set_printoptions(threshold=sys.maxsize)


### ================= argparse arguments =================================== ###
parser = argparse.ArgumentParser(description='Compare ct data with RNAfold')

parser.add_argument('--verbose','-v', action='store_true',
                    default=False,
                    help='be verbose')

parser.add_argument('--test', '-t', action='store_true',
                    default=False,
                    help='run with test csv')

parser.add_argument('--input','-i', action='store',
                    default="comp_m_n_mfe_prob.csv",
                    help='parse directory')

parser.add_argument('--output','-o', action='store_true',
                    default=False,
                    help='recalculate csvs')

args=parser.parse_args()

### ================= file-names =========================================== ###
if args.test==True:
    input_name = "comp_m_cp_mfe_prob.csv"
if args.test==False:
    input_name = args.input
output_name = input_name[:-4]
output_name = output_name + "_n_boot.csv"

### ================= functions ============================================ ###
def tpr(tp_data,fn_data):
    result = 0
    list_r = list()
    for i in range(0, len(tp_data)):
        result = tp_data[i] / (tp_data[i] + fn_data[i])
        list_r.append(result)
    return list_r

def ppv(tp_data,fp_data):
    result = 0
    list_r = list()
    for i in range(0, len(tp_data)):
        result = tp_data[i] / (tp_data[i] + fp_data[i])
        list_r.append(result)
    return list_r

def f_1(tp_data, fp_data, fn_data):
    result = 0
    list_r = list()
    f = 0
    for i in range(0, len(tp_data)):
        ppv = tp_data[i] / (tp_data[i] + fp_data[i])
        tpr = tp_data[i] / (tp_data[i] + fn_data[i])
        f = 2 * (ppv * tpr)/(ppv + tpr)
        result = f
        list_r.append(result)
    return list_r

def mcc(tp_data, fp_data, tn_data, fn_data):
    result = 0
    list_r = list()
    x = 0
    x = np.float64(x)
    root = 0
    for i in range(0, len(tp_data)):
        x = ((np.float64(tp_data[i]) + fp_data[i]) * (tp_data[i] + fn_data[i]) * (tn_data[i] + fp_data[i]) * (tn_data[i] + fn_data[i]))
        root = math.sqrt(x)
        mcc = (tp_data[i] * tn_data[i] - fp_data[i] * fn_data[i]) / root
        list_r.append(mcc)
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
            value = (data.iloc[rng])
            list_v.append(value)
        for k in range(0, len(list_v)):
            sum += list_v[k]
        list_s.append(sum)
        list_v = list()
        sum = 0
    return list_s

### ================= open data/ def vars ================================== ###
raw_data = pd.read_csv( open(input_name),float_precision = "high")
data_length = len(raw_data)

seql = raw_data['sequencelength']
drtr_tp_bpp = raw_data['drtr_tp']
drtr_fp_bpp = raw_data['drtr_fp']
drtr_tn_bpp = raw_data['drtr_tn']
drtr_fn_bpp = raw_data['drtr_fn']
drtr_tp_mfe = raw_data['drtr_tp_mfe']
drtr_fp_mfe = raw_data['drtr_fp_mfe']
drtr_tn_mfe = raw_data['drtr_tn_mfe']
drtr_fn_mfe = raw_data['drtr_fn_mfe']
drtr_tp_prob = raw_data['drtr_tp_prob']
drtr_fp_prob = raw_data['drtr_fp_prob']
drtr_tn_prob = raw_data['drtr_tn_prob']
drtr_fn_prob = raw_data['drtr_fn_prob']
fold_tp_mfe = raw_data['RNAfold_mfe_tp']
fold_fp_mfe = raw_data['RNAfold_mfe_fp']
fold_tn_mfe = raw_data['RNAfold_mfe_tn']
fold_fn_mfe = raw_data['RNAfold_mfe_fn']
fold_tp_bpp = raw_data['RNAfold_bpp_tp']
fold_fp_bpp = raw_data['RNAfold_bpp_fp']
fold_tn_bpp = raw_data['RNAfold_bpp_tn']
fold_fn_bpp = raw_data['RNAfold_bpp_fn']
len_data = len(seql)

### ================= calulate sums ======================================== ###
fold_tp_sum_mfe = boot(fold_tp_mfe)
fold_fp_sum_mfe = boot(fold_fp_mfe)
fold_tn_sum_mfe = boot(fold_tn_mfe)
fold_fn_sum_mfe = boot(fold_fn_mfe)

fold_tp_sum_bpp = boot(fold_tp_bpp)
fold_fp_sum_bpp = boot(fold_fp_bpp)
fold_tn_sum_bpp = boot(fold_tn_bpp)
fold_fn_sum_bpp = boot(fold_fn_bpp)

drtr_tp_sum_prob = boot(drtr_tp_prob)
drtr_fp_sum_prob = boot(drtr_fp_prob)
drtr_tn_sum_prob = boot(drtr_tn_prob)
drtr_fn_sum_prob = boot(drtr_fn_prob)

drtr_tp_sum_mfe = boot(drtr_tp_mfe)
drtr_fp_sum_mfe = boot(drtr_fp_mfe)
drtr_tn_sum_mfe = boot(drtr_tn_mfe)
drtr_fn_sum_mfe = boot(drtr_fn_mfe)

drtr_tp_sum_bpp = boot(drtr_tp_bpp)
drtr_fp_sum_bpp = boot(drtr_fp_bpp)
drtr_tn_sum_bpp = boot(drtr_tn_bpp)
drtr_fn_sum_bpp = boot(drtr_fn_bpp)

fold_tpr_mfe = tpr(fold_tp_sum_mfe, fold_fn_sum_mfe)
fold_tpr_bpp = tpr(fold_tp_sum_bpp, fold_fn_sum_bpp)
drtr_tpr_prob = tpr(drtr_tp_sum_prob, drtr_fn_sum_prob)
drtr_tpr_mfe = tpr(drtr_tp_sum_mfe, drtr_fn_sum_mfe)
drtr_tpr_bpp = tpr(drtr_tp_sum_bpp, drtr_fn_sum_bpp)

fold_ppv_mfe = ppv(fold_tp_sum_mfe, fold_fp_sum_mfe)
fold_ppv_bpp = ppv(fold_tp_sum_bpp, fold_fp_sum_bpp)
drtr_ppv_prob = ppv(drtr_tp_sum_prob, drtr_fp_sum_prob)
drtr_ppv_mfe = ppv(drtr_tp_sum_mfe, drtr_fp_sum_mfe)
drtr_ppv_bpp = ppv(drtr_tp_sum_bpp, drtr_fp_sum_bpp)

fold_f_1_mfe = f_1(fold_tp_sum_mfe, fold_fp_sum_mfe, fold_fn_sum_mfe)
fold_f_1_bpp = f_1(fold_tp_sum_bpp, fold_fp_sum_bpp, fold_fn_sum_bpp)
drtr_f_1_prob = f_1(drtr_tp_sum_prob, drtr_fp_sum_prob, drtr_fn_sum_prob)
drtr_f_1_mfe = f_1(drtr_tp_sum_mfe, drtr_fp_sum_mfe, drtr_fn_sum_mfe)
drtr_f_1_bpp = f_1(drtr_tp_sum_bpp, drtr_fp_sum_bpp, drtr_fn_sum_bpp)

fold_mcc_mfe = mcc(fold_tp_sum_mfe, fold_fp_sum_mfe, fold_tn_sum_mfe, fold_fn_sum_mfe)
fold_mcc_bpp = mcc(fold_tp_sum_bpp, fold_fp_sum_bpp, fold_tn_sum_bpp, fold_fn_sum_bpp)
drtr_mcc_prob = mcc(drtr_tp_sum_prob, drtr_fp_sum_prob, drtr_tn_sum_prob, drtr_fn_sum_prob)
drtr_mcc_mfe = mcc(drtr_tp_sum_mfe, drtr_fp_sum_mfe, drtr_tn_sum_mfe, drtr_fn_sum_mfe)
drtr_mcc_bpp = mcc(drtr_tp_sum_bpp, drtr_fp_sum_bpp, drtr_tn_sum_bpp, drtr_fn_sum_bpp)

names = (
"f_tpr_mfe", "f_tpr_bbp", "d_tpr_prob", "d_tpr_mfe", "d_tpr_bpp",
"f_ppv_mfe", "f_ppv_bbp", "d_ppv_prob", "d_ppv_mfe", "d_ppv_bpp",
"f_f1_mfe", "f_f1_bbp", "d_f1_prob", "d_f1_mfe", "d_f1_bpp",
"f_mcc_mfe", "f_mcc_bbp", "d_mcc_prob", "d_mcc_mfe", "d_mcc_bpp"
)

if args.output == True:
    data_processed = [
    fold_ppv_mfe,
    drtr_ppv_mfe,
    drtr_ppv_prob,
    fold_ppv_bpp,
    drtr_ppv_bpp,
    fold_tpr_mfe,
    drtr_tpr_mfe,
    drtr_tpr_prob,
    fold_tpr_bpp,
    drtr_tpr_bpp,
    fold_f_1_mfe,
    drtr_f_1_mfe,
    drtr_f_1_prob,
    fold_f_1_bpp,
    drtr_f_1_bpp,
    fold_mcc_mfe,
    drtr_mcc_mfe,
    drtr_mcc_prob,
    fold_mcc_bpp,
    drtr_mcc_bpp
    ]

    indexer = [
    "fold_ppv_mfe",
    "drtr_ppv_mfe",
    "drtr_ppv_prob",
    "fold_ppv_bpp",
    "drtr_ppv_bpp",
    "fold_tpr_mfe",
    "drtr_tpr_mfe",
    "drtr_tpr_prob",
    "fold_tpr_bpp",
    "drtr_tpr_bpp",
    "fold_f_1_mfe",
    "drtr_f_1_mfe",
    "drtr_f_1_prob",
    "fold_f_1_bpp",
    "drtr_f_1_bpp",
    "fold_mcc_mfe",
    "drtr_mcc_mfe",
    "drtr_mcc_prob",
    "fold_mcc_bpp",
    "drtr_mcc_bpp"
    ]

    cols = list()
    for i in range(0,boots):
        cols.append(i)

    data_processed = pd.DataFrame(data_processed, index = indexer, columns = cols)
    print(data_processed)
    # print(type(data_processed))
    export_csv = data_processed.to_csv(output_name)
