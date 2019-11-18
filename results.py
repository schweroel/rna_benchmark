#!python

### ================= import modules ======================================= ###
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl


### ================= settings for pandas and numpy ======================== ###
# prevents output from bein truncated when printed
pd.set_option('display.max_rows', 10000) #some default pandas settings will chop data
pd.set_option('display.max_columns', 10000)
pd.set_option('display.width', 10000)
np.set_printoptions(threshold=sys.maxsize) #if printed out, matrix won't be truncated

### ================= argparse arguments =================================== ###
parser = argparse.ArgumentParser(description='Compare ct data with RNAfold')
parser.add_argument('--drtrfile','-d', action='store',
                    default="drtr_ct_home.txt",
                    help='parse fasta ID')

parser.add_argument('--foldfile','-f', action='store',
                    default="fold_ct_canon.txt",
                    help='parse directory')

parser.add_argument('--output','-o', action='store',
                    default="c:/Users/filename",
                    help='parse directory')

args=parser.parse_args()



### ================= import of files ====================================== ###
drtr_name = args.drtrfile
fold_name = args.foldfile
output_name = args.output
output_name = output_name + ".csv"

# result-file of DrTrafo
data_drtr_ct = pd.read_csv( open(drtr_name) ,sep="\t", keep_default_na = False,
                                            skiprows = 0,
                                            names = ['ID',
                                                    'drtr_tp_mfe',
                                                    'drtr_fp_mfe',
                                                    'drtr_tn_mfe',
                                                    'drtr_fn_mfe',
                                                    'drtr_tp_prob',
                                                    'drtr_fp_prob',
                                                    'drtr_tn_prob',
                                                    'drtr_fn_prob',
                                                    'drtr_tp',
                                                    'drtr_fp',
                                                    'drtr_tn',
                                                    'drtr_fn',
                                                    'drtr_pps',
                                                    'drtr_ense_dist',
                                                    'drtr_ense_def',
                                                    'sequencelength'])


# ID.txt file contains IDs and corresponding RNA-family tag
family = pd.read_csv( open("ID.txt") ,keep_default_na = False, sep = ',' ,names = ['ID','family'])
family = pd.DataFrame(family)

# result-file of RNAfold
data_fold_ct = pd.read_csv( open(fold_name) ,sep="\t", keep_default_na = False,
                                            skiprows = 0,
                                            names = ['ID',              #0
                                            'RNAfold_mfe_tp',           #1
                                            'RNAfold_mfe_fp',           #2
                                            'RNAfold_mfe_tn',           #3
                                            'RNAfold_mfe_fn',           #4
                                            'RNAfold_bpp_tp',           #5
                                            'RNAfold_bpp_fp',           #6
                                            'RNAfold_bpp_tn',           #7
                                            'RNAfold_bpp_fn',           #8
                                            'RNAfold_pps',              #9
                                            'RNAfold_ens_dis',          #10
                                            'RNAfold_ense_def',         #11
                                            'sequencelength'])          #12

#define column titles
columns_titles = ['ID',                     #0
                'RNAfold_mfe_tp',           #1
                'RNAfold_mfe_fp',           #2
                'RNAfold_mfe_tn',           #3
                'RNAfold_mfe_fn',           #4
                'RNAfold_bpp_tp',           #5
                'RNAfold_bpp_fp',           #6
                'RNAfold_bpp_tn',           #7
                'RNAfold_bpp_fn',           #8
                'RNAfold_pps',              #9
                'RNAfold_ens_dis',          #10
                'RNAfold_ense_def',         #11
                'drtr_tp_mfe',              #12
                'drtr_fp_mfe',              #13
                'drtr_tn_mfe',              #14
                'drtr_fn_mfe',              #15
                'drtr_tp_prob',             #16
                'drtr_fp_prob',             #17
                'drtr_tn_prob',             #18
                'drtr_fn_prob',             #19
                'drtr_tp',                  #20 bpp
                'drtr_fp',                  #21 bpp
                'drtr_tn',                  #22 bpp
                'drtr_fn',                  #23 bpp
                'drtr_pps',                 #24
                'drtr_ense_dist',           #25
                'drtr_ense_def',            #26
                'sequencelength',           #27
                'family']                   #28

# reindexing of datatable
data_fold_ct = data_fold_ct.reindex( columns = data_fold_ct.columns.tolist() + ['family',           #13
                                                                                'drtr_tp_mfe',      #14
                                                                                'drtr_fp_mfe',      #15
                                                                                'drtr_tn_mfe',      #16
                                                                                'drtr_fn_mfe',      #17
                                                                                'drtr_tp_prob',     #18
                                                                                'drtr_fp_prob',     #19
                                                                                'drtr_tn_prob',     #20
                                                                                'drtr_fn_prob',     #21
                                                                                'drtr_tp',          #22
                                                                                'drtr_fp',          #23
                                                                                'drtr_tn',          #24
                                                                                'drtr_fn',          #25
                                                                                'drtr_pps',         #26
                                                                                'drtr_ense_dist',   #27
                                                                                'drtr_ense_def'])   #28


### ================= find family name from file =========================== ###
k = 0
for i in range(0, len(family)):
    for j in range(0, len(data_fold_ct)):
        if data_fold_ct.iloc[j,0] == family.iloc[i,0]:
            data_fold_ct.iloc[j,13] = family.iloc[i,1]
            j += 1
            break

### ================= concat results based on right ID ===================== ###
for p in range(0, len(data_drtr_ct)):
    for n in range(0, len(data_fold_ct)):
        if data_fold_ct.iloc[n,0] == data_drtr_ct.iloc[p,0]:
            data_fold_ct.iloc[n,14] = data_drtr_ct.iloc[p,1]
            data_fold_ct.iloc[n,15] = data_drtr_ct.iloc[p,2]
            data_fold_ct.iloc[n,16] = data_drtr_ct.iloc[p,3]
            data_fold_ct.iloc[n,17] = data_drtr_ct.iloc[p,4]
            data_fold_ct.iloc[n,18] = data_drtr_ct.iloc[p,5]
            data_fold_ct.iloc[n,19] = data_drtr_ct.iloc[p,6]
            data_fold_ct.iloc[n,20] = data_drtr_ct.iloc[p,7]
            data_fold_ct.iloc[n,21] = data_drtr_ct.iloc[p,8]
            data_fold_ct.iloc[n,22] = data_drtr_ct.iloc[p,9]
            data_fold_ct.iloc[n,23] = data_drtr_ct.iloc[p,10]
            data_fold_ct.iloc[n,24] = data_drtr_ct.iloc[p,11]
            data_fold_ct.iloc[n,25] = data_drtr_ct.iloc[p,12]
            data_fold_ct.iloc[n,26] = data_drtr_ct.iloc[p,13]
            data_fold_ct.iloc[n,27] = data_drtr_ct.iloc[p,14]
            data_fold_ct.iloc[n,28] = data_drtr_ct.iloc[p,15]

data_fold_ct = data_fold_ct.reindex(columns = columns_titles)

# export
export_csv = data_fold_ct.to_csv (output_name, header=True)
