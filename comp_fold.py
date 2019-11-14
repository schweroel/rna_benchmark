#!/bin/env python3
# content: open pickle-file that contains results from RNAfold and compare it
# with reference structure;
# can remove non-canonical-bps and pseudoknots when chosen;
# results can be piped into file

### ================= import modules ======================================= ###
import sys
import pandas as pd
import numpy as np
import forgi
import math as mt
import argparse
import pickle
from utilities import stuff as fus

### ================= settings for pandas and numpy ======================== ###
# prevents output from bein truncated when printed
pd.set_option('display.max_rows', 10000)
pd.set_option('display.max_columns', 10000)
pd.set_option('display.width', 10000)
np.set_printoptions(threshold=sys.maxsize)

### ================= argparse arguments =================================== ###
parser = argparse.ArgumentParser(description='Compare ct data with RNAfold')
parser.add_argument('--test','-t', action='store_true',
                    default=False,
                    help='run with testseq')
parser.add_argument('--onlycanonical','-c', action='store_true',
                    default=False,
                    help='corrects input reference structure to canonical bps only, default: OFF')
parser.add_argument('--verbose','-v', action='store_true',
                    default=False,
                    help='show matrices and all extra-information')
parser.add_argument('--use-loops','-l', action='store_true',
                    default=False,
                    help='Use loops instead of functions')
parser.add_argument('--identifier','-i', action='store',
                    default="ASE_00089",
                    help='parse fasta ID')
parser.add_argument('--directory','-d', action='store',
                    default=".",
                    help='parse directory')
parser.add_argument('--recompute','-r', action='store_true',
                    default=False,
                    help='recompute instead of reding pre-computed file')

args=parser.parse_args()

### ============= get file-names from argparse ============================= ###
if args.recompute == True:
    filename_ct = args.directory + "/" + args.identifier + ".ct"
    filename_fasta = args.directory + "/" + args.identifier + ".fasta"

if args.recompute == False:
    filename_fasta = "/home/pickledata/" + args.identifier + ".pickle"
    ID = args.identifier

if args.test == False:
    filename_ct = args.directory + "/" + args.identifier + ".ct"
    name = args.identifier

if args.test == True:
    filename_ct = "/home/ASE_00035.ct"
    filename_fasta = "/home/ASE_00035.fasta"
    name = filename_ct[48:57]

# Flag to indicate whether we use functions (True) or on-the-fly computing (False)
use_functions = True

if args.use_loops:
    use_functions = False


### ================= functions ============================================ ###
def performance_mfe(pt_reference, pt_mfe):
    n  = pt_reference[0]
    tp = 0
    fp = 0
    tn = n * (n - 1) / 2
    fn = 0

    # compute true positives and false negatives
    for i in range(1, n + 1):
        # we have a pair
        if pt_reference[i] > i:
            if pt_mfe[i] == pt_reference[i]:
                tp = tp + 1
            else:
                fn = fn + 1

    # compute false positive
    for i in range(1, n + 1):
        # we have a pair in pt_mfe but not in pt_reference
        if pt_mfe[i] > i and pt_reference[i] != pt_mfe[i]:
            fp = fp + 1

    return (tp, fp, tn, fn)


def performance_bpp(pt_reference, bpp_matrix):
    n  = pt_reference[0]
    tp = 0
    fp = 0
    tn = 0
    fn = 0

    # compute true positives and false negatives
    for i in range(1, n + 1):
        # we have a pair
        if pt_reference[i] > i:
            tp = tp + bpp_matrix[i][pt_reference[i]]
            fn = fn + 1 - bpp_matrix[i][pt_reference[i]]

    # compute false positives and true negatives
    for i in range(1, n +1):
        for j in range(i + 1, n + 1):
            if pt_reference[i] != j:
                fp = fp + bpp_matrix[i][j]
                tn = tn + 1 - bpp_matrix[i][j]

    return (tp, fp, tn, fn)

def pair_prob_score(pt_reference, bpp_matrix):
    n = pt_reference[0]
    score = 0
    num_bp = 0
    for i in range(1, n + 1):
        if pt_reference[i] > i:
            score = score + bpp_matrix[i][pt_reference[i]]
            num_bp += 1

    score = score / num_bp
    return score

def ensemble_distance(pt_reference, bpp_matrix):
    n = pt_reference[0]
    distance = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            if pt_reference[i] == j:
                distance = distance + 1 - bpp_matrix[i][j]
            else:
                distance = distance + bpp_matrix[i][j]
    return distance

def ensemble_defect(pt_reference, bpp_matrix):
    n = pt_reference[0]
    # si = 0
    score = 0
    defect = 0
    for i in range(1, n + 1):
        # we have a pair
        if pt_reference[i] > i:
            score += bpp_matrix[i][pt_reference[i]]
        elif pt_reference[i] == 0:
            q = 0
            for j in range(1, i):
                q += bpp_matrix[j][i]
            for j in range(i + 1, n + 1):
                q += bpp_matrix[i][j]
            if q > 1:
                q = 1
            q = 1 - q
            score += q
        else:
            score += bpp_matrix[pt_reference[i]][i]

    score /= n
    defect = 1 - score
    return defect

### ================= function to alter modified nucleobases =============== ###
def rna_mod(sequence, mask_modified = False, mask_RTblocking = False, uppercase = True):
    adenosines = ('A', 'a', '?A', 'm1A', 'm2A', 'i6A', 'ms2i6A', 'm6A', 't6A', 'm6t6A', 'ms2t6A', 'Am', 'I', 'm1I', 'Ar', 'io6A','H','"','/', '+', '*','=','6','E')
    cytidines = ('C', 'c', '?C', 's2C', 'Cm', 'ac4C', 'm5C', 'm3C', 'k2C', 'f5C', 'f5Cm' )
    guanosines = ('G', 'g', '?G', 'm1G', 'm2G', 'Gm', 'm22G', 'm22Gm', 'm7G', 'fa7d7G', 'Q', 'manQ', 'galQ', 'yW o2yW' )
    uridines = ( 'U', 'u', 'T', 't', '?U', 'mnm5U', 's2U', 'Um', 's4U', 'ncm5U', 'mcm5U',
                    'mnm5s2U', 'mcm5s2U', 'cmo5U', 'mo5U', 'cmnm5U', 'cmnm5s2U', 'acp3U', 'mchm5U', 'cmnm5Um', 'ncm5Um', 'D', 'psi', 'm1psi', 'psim', 'm5U', 'm5s2U', 'm5Um' )
    RTBlocking = ( 'm1A' 'I' 'm1G' 'm22G' 'm3U' 'm3psi' 'acp3U' 'm3C' 'Q')

    mods = {
	'.': 'N', # unknown nucleotide
	'_': 'N', # insertion
    # Modified Adenosines
	'H': '?A', # unknown modified adenosine
	'"': 'm1A', # 1-methyladenosine -> blocks RT
	'/': 'm2A', # 2-methyladenosine
	'+': 'i6A', # N6-isopentenyladenosine
	'*': 'ms2i6A', # 2-methylthio-N6-isopentenyladenosine
	'=': 'm6A', # N6-methyladenosine
	'6': 't6A', # N6-threonylcarbamoyladenosine
	'E': 'm6t6A', # N6-methyl-N6-threonylcarbamoyladenosine
	'[': 'ms2t6A', # 2-methylthio-N6-threonylcarbamoyladenosine
	':': 'Am', # 2'-O-methyladenosine
	'I': 'I', # inosine -> blocks RT
	'O': 'm1I', # 1-methylinosine
	'^': 'Ar', # 2'-O-ribosyladenosine (phosphat)
	'`': 'io6A', # N6-(cis-hydroxyisopentenyl)adenosine
    # Modified Cytidines
	'<': '?C', # unknown modified cytidine
	'%': 's2C', # 2-thiocytidine
	'B': 'Cm', # 2'-O-methylcytidine
	'M': 'ac4C', # N4-acetylcytidine
	'?': 'm5C', # 5-methylcytidine
	'\'': 'm3C', # 3-methylcytidine
	'}': 'k2C', # lysidine
	'>': 'f5C', # 5-formylcytidin
	'°': 'f5Cm', # 2-O-methyl-5-formylcytidin
    # Modified Guanosines
	';': '?G', # unknown modified guanosine
	'K': 'm1G', # 1-methylguanosine
	'L': 'm2G', # N2-methylguanosine
	'#': 'Gm', # 2'-O-methylguanosine
	'R': 'm22G', # N2,N2-dimethylguanosine
	'|': 'm22Gm', # N2,N2,2'-O-trimethylguanosine
	'7': 'm7G', # 7-methylguanosine
	'(': 'fa7d7G', # archaeosine
	'Q': 'Q', # queuosine
	'8': 'manQ', # mannosyl-queuosine
	'9': 'galQ', # galactosyl-queuosine
	'Y': 'yW', # wybutosine
	'W': 'o2yW', # peroxywybutosineadenosines
    # Modified Uridines
	'N': '?U', # unknown modified uridine
	'{': 'mnm5U', # 5-methylaminomethyluridine
	'2': 's2U', # 2-thiouridine
	'J': 'Um', # 2'-O-methyluridine
	'4': 's4U', # 4-thiouridine
	'&': 'ncm5U', # 5-carbamoylmethyluridine
	'1': 'mcm5U', # 5-methoxycarbonylmethyluridine
	'S': 'mnm5s2U', # 5-methylaminomethyl-2-thiouridine
	'3': 'mcm5s2U', # 5-methoxycarbonylmethyl-2-thiouridine
	'V': 'cmo5U', # uridine 5-oxyacetic acid
	'5': 'mo5U', # 5-methoxyuridine
	'!': 'cmnm5U', # 5-carboxymethylaminomethyluridine
	'$': 'cmnm5s2U', # 5-carboxymethylaminomethyl-2-thiouridine
	'X': 'acp3U', # 3-(3-amino-3-carboxypropyl)uridine
	',': 'mchm5U', # 5-(carboxyhydroxymethyl)uridinemethyl ester
	')': 'cmnm5Um', # 5-carboxymethylaminomethyl-2'-O-methyluridine
	'~': 'ncm5Um', # 5-carbamoylmethyl-2'-O-methyluridine
	'D': 'D', # dihydrouridine
	'P': 'psi', # pseudouridine
	']': 'm1psi', # 1-methylpseudouridine
	'Z': 'psim', # 2'-O-methylpseudouridine
	'T': 'm5U', # ribosylthymine
	'F': 'm5s2U', # 5-methyl-2-thiouridine
	'\\': 'm5Um' # 5, 2'-O-dimethyluridine
    }

    unmod_sequence = ""
    for i in range(0, len(sequence)):
        c = sequence[i]

        if c in ("A", "C", "G", "U", "T", "N", "a", "c", "g", "u", "t" , "n"):
            if uppercase:
                unmod_sequence += c.upper()
            else:
                unmod_sequence += c
        elif mask_modified:
            unmod_sequence += "N"
        elif c in mods:
            a = mods[c]

            if mask_RTblocking and a in RTBlocking:
                unmod_sequence += "N"
            elif a in adenosines:
                unmod_sequence += "A"
            elif a in cytidines:
                unmod_sequence += "C"
            elif a in guanosines:
                unmod_sequence += "G"
            elif a in uridines:
                unmod_sequence += "U"
            else:
                unmod_sequence += "N"
        else:
            unmod_sequence += "N"

    return unmod_sequence


### ================= header counter ctfile + import ======================= ###
counter = 0
f=open(filename_ct)
for x in f:
    if x[:1] == "#":
        counter += 1
    else:
        counter = counter + 1
        data_ct = pd.read_csv( open(filename_ct) ,sep="\s+", engine='python', keep_default_na = False, skiprows = counter, names = ['n', 'Base','n-1','n+1','bp','nn'])
        data_corrected_ct = data_ct.drop(['n','n-1','n+1'], axis = 1)
        break


### ================= count sequence length ================================ ###
col_ct = len(data_corrected_ct.index)

### ================= make pairtable from ct =============================== ###
bpdata = [col_ct]
sequence = ""
for i in range(0, col_ct):
    base = (data_corrected_ct.iloc[i,0])
    bp = int(data_corrected_ct.iloc[i,1])
    nn = int(data_corrected_ct.iloc[i,2])
    if bp == 0:
        bpdata.append(0)
    else:
        bpdata.append(int(bp))

    sequence += base

### ================= remove non-canonical base pairs if required ========== ###
if args.onlycanonical == True:
    unmod_sequence = rna_mod(sequence)

    for i in range(1, col_ct):
        if bpdata[i] > i:
            bp = (unmod_sequence[i - 1], unmod_sequence[bpdata[i] - 1])
            if bp not in [("G","C"), ("C","G"), ("A","U"), ("U","A"), ("G","U"), ("U","G")]:
                bpdata[bpdata[i]] = 0
                bpdata[i] = 0

### ================= importer for fasta files ============================= ###
if args.recompute == True:
    f = open(filename_fasta, "r")
    bla = (f.readline())
    seq = (f.readline())
    seq = seq[:-1]
    fc = RNA.fold_compound(seq)
    (ss, mfe) = fc.mfe()
    fc.exp_params_rescale(mfe)
    (propensity,ensemble_energy) = fc.pf()
    matrix_prob2 = fc.bpp()

if args.recompute == False:
    name = ID
    pickle_file_name = ID + ".pickle"
    pickle_file_path = r'/home/pickledata/' + pickle_file_name
    with open(pickle_file_path, 'rb') as e:
        data = pickle.load(e)
    ss = data['db']
    mfe = data['mfe']
    matrix_prob2 = data['bppmatrix']

### ================= RNAfold-section ====================================== ###
matrix_index = list()
for i in range(0, col_ct + 1):
    matrix_index.append(float(i))

### ================= confusion matrix for bpp ============================= ###
bpp_fp = 0
bpp_tn = 0
bpp_fn = 0
bpp_tp = 0

if use_functions:
    (bpp_tp, bpp_fp, bpp_tn, bpp_fn) = performance_bpp(bpdata, matrix_prob2)
else:
    for i in range(1, col_ct+1):
        if bpdata[i] == 0:
            bpp_fp_temp = 0
            bpp_fp_tup = list(zip(*matrix_prob2))[i]
            bpp_fp_temp += sum(x for x in bpp_fp_tup)
            bpp_fp += bpp_fp_tem
            if bpp_fp_temp > 0:
                bpp_tn += (1 - bpp_fp_temp)
            else:
                bpp_tn += 0
        else:
            match = bpdata[i]
            bpp_tp_temp = 0
            bpp_tp_temp += (list(zip(*matrix_prob2))[i][match])
            bpp_tp += bpp_tp_temp
            if bpp_tp_temp > 0:
                bpp_fn += (1 - bpp_tp_temp)
            else:
                bpp_fn += 0

### ================= MFE confusion matrix ================================= ###
###loop version
pair_table_pred = forgi.dotbracket_to_pairtable(ss)
mfe_tn = 0
mfe_tp = 0
mfe_fp = 0
mfe_fn = 0
if use_functions:
    (mfe_tp, mfe_fp, mfe_tn, mfe_fn) = performance_mfe(bpdata, pair_table_pred)
else:
    for i in range(1,col_ct + 1):
        if bpdata[i] == 0:
            if pair_table_pred[i] == 0:
                mfe_tn += 1
            else:
                mfe_fp += 1
        else:
            if pair_table_pred[i] == bpdata[i]:
                mfe_tp += 1
            else:
                mfe_fn += 1

    mfe_tp = mfe_tp / 2
    mfe_fp = mfe_fp / 2
    mfe_fn = mfe_fn / 2
    mfe_tn = (col_ct*(col_ct-1)/2)

### ================= PPS Section ========================================== ###
contributing_bp = 0
score_bpp_pred = 0

if use_functions:
    pps = pair_prob_score(bpdata, matrix_prob2)
else:
    for i in range(1,col_ct + 1):
        if bpdata[i] != 0:
            index_pps = bpdata[i]
            contributing_bp += 1
            score_bpp_pred += (list(zip(*matrix_prob2))[i][index_pps])

    contributing_bp = contributing_bp / 2
    if contributing_bp > 0:
        pps = ( score_bpp_pred/ contributing_bp)
    else:
        pps = 0

### ================= ensemble distance / defect =========================== ###
# The ensemble defect is the average number of incorrectly paired nucleotides at equilibrium over the ensemble of (possible) secondary structures(=shape).
#The Ensemble Distance is the L 1 -distance between the predicted ensemble and the reference structure:

ense_def = 0
ense_dist = 0
pk_bpp_fp = 0
pk_bpp_tn = 0
pk_bpp_fn = 0
pk_bpp_tp = 0

for i in range(1, col_ct+1):
    k = i + 1
    k2= k + 1
    if bpdata_pkrem[i] == 0:
        pk_bpp_fp_temp = 0
        pk_bpp_fp_tup = list(zip(*matrix_prob2))[i]
        pk_bpp_fp_temp += sum(x for x in pk_bpp_fp_tup)
        pk_bpp_tn += (1 - pk_bpp_fp_temp)
        pk_bpp_fp += pk_bpp_fp_temp
    else:
        match = bpdata_pkrem[i]
        pk_bpp_tp_temp = 0
        pk_bpp_tp_temp += (list(zip(*matrix_prob2))[i][match])
        pk_bpp_fn += (1 - pk_bpp_tp_temp)
        pk_bpp_tp += pk_bpp_tp_temp
        ense_dist += (1-pk_bpp_tp) + pk_bpp_fp
        ense_def += pk_bpp_fp

if args.recompute == True:
    ense_def = fc.ensemble_defect(reference_db)
    ense_dist = ensemble_distance(bpdata, matrix_prob2)

if args.recompute == False:
    ense_dist = ensemble_distance(bpdata, matrix_prob2)
    ense_def = ensemble_defect(bpdata, matrix_prob2)

### ================= verbose print - section ============================== ###
if args.verbose == True:
    print(filename_ct)
    print(filename_fasta)
    print("number of nucleotides>")
    print(col_ct)
    print("pair_table pk_free:")
    print(bpdata_pkrem)
    print(len(bpdata_pkrem))
    if args.onlycanonical == True:
        print("pt only canon>")
        print(bpdata)
        print(len(bpdata))
    if args.onlycanonical == False:
        print("pair_table raw>")
        print(bpdata)
        print(len(bpdata))
    print("predicted pair table RNAfold>")
    print(pair_table_pred)
    print("mfe_tp>")
    print(mfe_tp)
    print("mfe_fp>")
    print(mfe_fp)
    print("mfe_tn>")
    print(mfe_tn)
    print("mfe_fn>")
    print(mfe_fn)
    # print("probability matrix")
    # ind = list(range(0,col_ct+1))
    # print(ind)
    # for i in range (0,col_ct):
        # print(*matrix_prob2[i:i+1])
    print("confusion matrix bpp:")
    #true_positives
    print("bpp_tp>")
    print(bpp_tp)
    #false_positives
    print("bpp_fp>")
    print(bpp_fp)
    #true_negatives
    #alt:
    # bpp_tn = ((col_ct*(col_ct-1))/2)
    print("bpp_tn>")
    print(bpp_tn)
    #false_negatives
    print("bpp_fn>")
    print(bpp_fn)
    print("")
    print("contributing_bp>")
    print(contributing_bp)
    print("pps>")
    print(pps)
    print("")
    print("ense_def")
    print(ense_def)
    print("ense_dist")
    print(ense_dist)

### ===================== results into table - section ===================== ###
print("%s\t"
      "%d\t%d\t%d\t%d\t"
      "%8.5f\t%8.5f\t%8.5f\t%8.5f\t"
      "%8.8f\t"
      "%8.8f\t"
      "%8.8f\t"
      "%d" % (name,
            mfe_tp,
            mfe_fp,
            mfe_tn,
            mfe_fn,
            bpp_tp,
            bpp_fp,
            bpp_tn,
            bpp_fn,
            pps,
            ense_dist,
            ense_def,
            col_ct))
