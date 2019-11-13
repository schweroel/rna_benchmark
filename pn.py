#!python
# content: script that calulates tpr/ppv/f-measure and mcc from raw results;
# calculates values for the whole dataset

### ================= import modules ======================================= ###
import numpy as np
import pandas as pd
import math
import sys
import matplotlib.pyplot as plt
import argparse
from matplotlib.pyplot import figure

### ================= settings for pandas and numpy ======================== ###
# prevents output from bein truncated when printed
pd.set_option('display.max_rows', 10000)
pd.set_option('display.max_columns', 10000)
pd.set_option('display.width', 10000)
np.set_printoptions(threshold=sys.maxsize)


### ================= argparse arguments =================================== ###
parser = argparse.ArgumentParser(description='calculate tpr,ppv,mcc,f_measure for given result csv')

parser.add_argument('--verbose','-v', action='store_true',
                    default=False,
                    help='be verbose')

parser.add_argument('--graphics','-g', action='store_true',
                    default=False,
                    help='print graphics')

parser.add_argument('--label','-l', action='store',
                    default="",
                    help='define label for graphics')

parser.add_argument('--test', '-t', action='store_true',
                    default=False,
                    help='run with test csv')

parser.add_argument('--input','-i', action='store',
                    default="comp_m_n_mfe_prob.csv",
                    help='parse directory')

parser.add_argument('--output','-o', action='store_true',
                    default=False,
                    help='make csv')

parser.add_argument('--printout','-p', action='store_true',
                    default=False,
                    help='print output')

args=parser.parse_args()

if args.test == False:
    input_name = args.input

if args.test == True:
    input_name = "comp_m_n_mfe_prob.csv"

output_name = input_name[:-4] + "_n_mean.csv"

label_name = input_name[:-4]

### ================= functions ============================================ ###

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
        tp = tp_data.iloc[i]
        fp = fp_data.iloc[i]
        # print(tp)
        tp_sum += tp
        fp_sum += fp
        # print(tp_sum)
    result = tp_sum / (tp_sum + fp_sum)
    # print(tp_sum)
    # print(fp_sum)
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
raw_data = pd.read_csv( open(input_name))
data_length = len(raw_data)

seql = raw_data['sequencelength']
pps_fold = raw_data['RNAfold_pps']
pps_drtr = raw_data['drtr_pps']

drtr_tp = raw_data['drtr_tp']
drtr_fp = raw_data['drtr_fp']
drtr_tn = raw_data['drtr_tn']
drtr_fn = raw_data['drtr_fn']

drtr_tp_mfe = raw_data['drtr_tp_mfe']
drtr_fp_mfe = raw_data['drtr_fp_mfe']
drtr_tn_mfe = raw_data['drtr_tn_mfe']
drtr_fn_mfe = raw_data['drtr_fn_mfe']

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

ppv_fold_mfe = ppv(RNAfold_mfe_tp, RNAfold_mfe_fp)
ppv_fold_bpp = ppv(RNAfold_bpp_tp, RNAfold_bpp_fp)
ppv_drtr = ppv(drtr_tp, drtr_fp)
ppv_drtr_mfe = ppv(drtr_tp_mfe, drtr_fp_mfe)
ppv_drtr_prob = ppv(drtr_tp_prob, drtr_fp_prob)

tpr_fold_mfe = tpr(RNAfold_mfe_tp, RNAfold_mfe_fn)
tpr_fold_bpp = tpr(RNAfold_bpp_tp, RNAfold_bpp_fn)
tpr_drtr = tpr(drtr_tp, drtr_fn)
tpr_drtr_mfe = tpr(drtr_tp_mfe, drtr_fn_mfe)
tpr_drtr_prob = tpr(drtr_tp_prob, drtr_fn_prob)


f_measure_fold_mfe = f_measure(RNAfold_mfe_tp,RNAfold_mfe_fp,RNAfold_mfe_fn)
f_measure_fold_bpp = f_measure(RNAfold_bpp_tp, RNAfold_bpp_fp, RNAfold_bpp_fn)
f_measure_drtr = f_measure(drtr_tp, drtr_fp, drtr_fn)
f_measure_drtr_mfe = f_measure(drtr_tp_mfe, drtr_fp_mfe, drtr_fn_mfe)
f_measure_drtr_prob = f_measure(drtr_tp_prob, drtr_fp_prob, drtr_fn_prob)


mcc_fold_mfe = mcc(RNAfold_mfe_tp, RNAfold_mfe_fp, RNAfold_mfe_tn, RNAfold_mfe_fn)
mcc_fold_bpp = mcc(RNAfold_bpp_tp, RNAfold_bpp_fp, RNAfold_bpp_tn, RNAfold_bpp_fn)
mcc_drtr = mcc(drtr_tp, drtr_fp, drtr_tn, drtr_fn)
mcc_drtr_mfe = mcc(drtr_tp_mfe, drtr_fp_mfe, drtr_tn_mfe, drtr_fn_mfe)
mcc_drtr_prob = mcc(drtr_tp_prob, drtr_fp_prob, drtr_tn_prob, drtr_fn_prob)


rows=("ppv_fold_mfe",
      "ppv_drtr_mfe",
      "ppv_drtr_prob",
      "ppv_fold_bpp",
      "ppv_drtr",
      "tpr_fold_mfe",
      "tpr_drtr_mfe",
      "tpr_drtr_prob",
      "tpr_fold_bpp",
      "tpr_drtr",
      "f_measure_fold_mfe",
      "f_measure_drtr_mfe",
      "f_measure_drtr_prob",
      "f_measure_fold_bpp",
      "f_measure_drtr",
      "mcc_fold_mfe",
      "mcc_drtr_mfe",
      "mcc_drtr_prob",
      "mcc_fold_bpp",
      "mcc_drtr")

if args.printout == True:
    print(rows)
    print(
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
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    "%8.5f\t"
    % (
    ppv_fold_mfe,
    ppv_drtr_mfe,
    ppv_drtr_prob,
    ppv_fold_bpp,
    ppv_drtr,
    tpr_fold_mfe,
    tpr_drtr_mfe,
    tpr_drtr_prob,
    tpr_fold_bpp,
    tpr_drtr,
    f_measure_fold_mfe,
    f_measure_drtr_mfe,
    f_measure_drtr_prob,
    f_measure_fold_bpp,
    f_measure_drtr,
    mcc_fold_mfe,
    mcc_drtr_mfe,
    mcc_drtr_mfe,
    mcc_fold_bpp,
    mcc_drtr))

if args.graphics == True:

    fold_mfes = (
    ppv_fold_mfe,
    tpr_fold_mfe,
    f_measure_fold_mfe,
    mcc_fold_mfe
    )

    fold_bpps = (
    ppv_fold_bpp,
    tpr_fold_bpp,
    f_measure_fold_bpp,
    mcc_fold_bpp
    )

    drtr_mfes = (
    ppv_drtr_mfe,
    tpr_fold_mfe,
    f_measure_drtr_mfe,
    mcc_fold_mfe
    )

    drtr_probs = (
    ppv_drtr_prob,
    tpr_drtr_prob,
    f_measure_drtr_prob,
    mcc_drtr_prob
    )

    drtr_bpps = (
    ppv_drtr,
    tpr_drtr,
    f_measure_drtr,
    mcc_drtr
    )


    titel = args.label
    barwidth = 0.1

    r_1 = np.arange(len(drtr_mfes))
    r_2 = [x + barwidth for x in r_1]
    r_3 = [x + barwidth for x in r_2]
    r_4 = [x + barwidth for x in r_3]
    r_5 = [x + barwidth for x in r_4]

    plt.bar(r_1, fold_mfes, width=barwidth, color='#247E85', label='RNAfold: mfe')
    plt.bar(r_2, drtr_mfes, width=barwidth, color='#87f542', label='DrTransformer: mfe')
    plt.bar(r_3, drtr_probs, width=barwidth, color='#FA5882', label='DrTransformer: probs')
    plt.bar(r_4, fold_bpps, width=barwidth, color='#0E3151', label='RNAfold: bpp')
    plt.bar(r_5, fold_bpps, width=barwidth, color='#BA3151', label='DrTransformer: bpp')

    plt.title(titel, fontsize=20, y=1.03)
    plt.xlabel('Family', fontsize=14)
    plt.ylabel('MCC',  fontsize=14)
    plt.xticks([r + barwidth for r in range (len(drtr_mfes))],
                                                   ["ppv", "tpr", "f_measure", "mcc"],

                                                   rotation=55,
                                                   ha="right",
                                                   rotation_mode="anchor",
                                                   fontsize=12)

    plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    plt.gcf().subplots_adjust(bottom=0.28, left=0.10, right=0.94, top=0.91)
    # plt.grid()
    plt.legend()
    plt.show()


if args.output == True:
    data_x = [
    ppv_fold_mfe,
    ppv_drtr_mfe,
    ppv_drtr_prob,
    ppv_fold_bpp,
    ppv_drtr,
    tpr_fold_mfe,
    tpr_drtr_mfe,
    tpr_drtr_prob,
    tpr_fold_bpp,
    tpr_drtr,
    f_measure_fold_mfe,
    f_measure_drtr_mfe,
    f_measure_drtr_prob,
    f_measure_fold_bpp,
    f_measure_drtr,
    mcc_fold_mfe,
    mcc_drtr_mfe,
    mcc_drtr_prob,
    mcc_fold_bpp,
    mcc_drtr]
    print(output_name)
    data_final = pd.DataFrame(data_x)
    # print(data_final)
    export_csv = data_final.to_csv(output_name)
