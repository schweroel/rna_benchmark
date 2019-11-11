#!python
# content: script to combine results from RNAfold and DrTransformer into one .csv file
# also adds RNAfamily to the csv, when ID.txt is provided

### ================= import modules ======================================= ###
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl

### ================= settings for pandas and numpy ======================== ###
# prevents output from bein truncated when printed
pd.set_option('display.max_rows', 10000)
pd.set_option('display.max_columns', 10000)
pd.set_option('display.width', 10000)
np.set_printoptions(threshold=sys.maxsize)

### ================= argparse arguments =================================== ###
parser = argparse.ArgumentParser(description='script to combine results of DrTrafo w/ RNAfold')
parser.add_argument('--drtrfile','-d', action='store',
                    default="drtr_ct_home.txt",
                    help='define input-file of DrTrafo')

parser.add_argument('--foldfile','-f', action='store',
                    default="fold_ct_canon.txt",
                    help='define input-file of RNAfold')

parser.add_argument('--output','-o', action='store',
                    default="c:/Users",
                    help='define output-directory')

args=parser.parse_args()


### ================= import of files =======================================###
drtr_name = args.drtrfile
fold_name = args.foldfile
output_name = args.output
output_name = output_name + ".csv"

family = pd.read_csv( open("ID.txt") ,keep_default_na = False, sep = ',' ,names = ['ID','family'])
family = pd.DataFrame(family)

data_drtr_ct = pd.read_csv( open(drtr_name) ,sep="\t", keep_default_na = False,
                                            skiprows = 0,
                                            names = ['ID',
                                                    'drtr_tp',
                                                    'drtr_fp',
                                                    'drtr_tn',
                                                    'drtr_fn',
                                                    'drtr_pps',
                                                    'drtr_ense_dist',
                                                    'drtr_ense_def',
                                                    'sequencelength'])

data_fold_ct = pd.read_csv( open(fold_name) ,sep="\t", keep_default_na = False,
                                            skiprows = 0,
                                            names = ['ID',
                                            'RNAfold_mfe_tp',
                                            'RNAfold_mfe_fp',
                                            'RNAfold_mfe_tn',
                                            'RNAfold_mfe_fn',
                                            'RNAfold_bpp_tp',
                                            'RNAfold_bpp_fp',
                                            'RNAfold_bpp_tn',
                                            'RNAfold_bpp_fn',
                                            'RNAfold_pps',
                                            'RNAfold_ens_dis',
                                            'RNAfold_ense_def',
                                            'sequencelength'])

columns_titles = ['ID',
                'RNAfold_mfe_tp',
                'RNAfold_mfe_fp',
                'RNAfold_mfe_tn',
                'RNAfold_mfe_fn',
                'RNAfold_bpp_tp',
                'RNAfold_bpp_fp',
                'RNAfold_bpp_tn',
                'RNAfold_bpp_fn',
                'RNAfold_pps',
                'RNAfold_ens_dis',
                'RNAfold_ense_def',
                'drtr_tp',
                'drtr_fp',
                'drtr_tn',
                'drtr_fn',
                'drtr_pps',
                'drtr_ense_dist',
                'drtr_ense_def',
                'sequencelength',
                'family']

data_fold_ct = data_fold_ct.reindex( columns = data_fold_ct.columns.tolist() + ['family',
                                                                                'drtr_tp',
                                                                                'drtr_fp',
                                                                                'drtr_tn',
                                                                                'drtr_fn',
                                                                                'drtr_pps',
                                                                                'drtr_ense_dist',
                                                                                'drtr_ense_def'])


### ================= find RNAfamily from ID.txt file ====================== ###
k = 0
for i in range(0, len(family)):
    for j in range(0, len(data_fold_ct)):
        if data_fold_ct.iloc[j,0] == family.iloc[i,0]:
            data_fold_ct.iloc[j,13] = family.iloc[i,1]
            j += 1
            break

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

### ================= Reindex data table and export ======================== ###
data_fold_ct = data_fold_ct.reindex(columns = columns_titles)
export_csv = data_fold_ct.to_csv (output_name, header=True)
