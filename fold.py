#!/bin/env python2

# content: python script, that gets a directory and fasta-file passed;
# for this sequence a secondary structure prediction is performed by
# the python interface of RNAfold
# the ID, the MFE in db-notation, the bpp-matrix, the primary sequence and
# the MFE are saved in a .pickle file, named after the ID of the RNA

### ================= import modules ======================================= ###
import sys
import RNA
import pickle
import re
import argparse

### ================= argparse arguments =================================== ###
parser = argparse.ArgumentParser(description='Compare ct data with RNAfold')

parser.add_argument('--identifier','-i', action='store',
                    default="ASE_00089",
                    help='parse fasta ID')
parser.add_argument('--directory','-d', action='store',
                    default=".",
                    help='parse directory')
args=parser.parse_args()

### ================= import of fasta-file + output-file-name ============== ###
filename_fasta = args.directory + "/" + args.identifier + ".fasta"
name = args.identifier
ID = name

f = open(filename_fasta, "r")
skip = (f.readline())
seq = (f.readline())

pickle_file_name = ID + ".pickle"
pickle_path = "/home/" + pickle_file_name

### ================= RNAfold-section ====================================== ###
fc = RNA.fold_compound(seq)
(ss, mfe) = fc.mfe()
fc.exp_params_rescale(mfe)
(propensity,ensemble_energy) = fc.pf()
matrix_prob = fc.bpp()

### ================= Pickle-export ======================================== ###
data_fold = {
    'name': ID,
    'seq': seq,
    'mfe': mfe,
    'db' : ss,
    'bppmatrix': matrix_prob
    }

with open (pickle_path, 'wb') as e:
    pickle.dump( data_fold, e, pickle.HIGHEST_PROTOCOL)
