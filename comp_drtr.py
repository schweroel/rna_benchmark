#!/usr/bin python3

# content: takes a .log-file of DrTransformer, performs a readout and
# compares it with a reference structure, that has the same ID
# can remove non-canonical base pairs and pseudoknots from reference
# structure

### ================= import modules ======================================= ###
import sys
import numpy as np
import RNA
import pandas as pd
import argparse
import re

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

parser.add_argument('--canonical', '-c', action='store_true',
                    default=False,
                    help='corrects input ct to canonical only')

parser.add_argument('--verbose','-v', action='store_true',
                    default=False,
                    help='be verbose and show all extra-information, default=OFF')

parser.add_argument('--identifier','-i', action='store',
                    default="ASE_00089",
                    help='parse fasta ID')

parser.add_argument('--directory','-d', action='store',
                    default=".",
                    help='parse directory')

parser.add_argument('--pkrem','-p', action='store_true',
                    default=False,
                    help='activate removal of pseudoknots')


args=parser.parse_args()

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


### ================= get filename from argparse =========================== ###
filename_ct = args.directory + "/" + args.identifier + ".ct"
name = args.identifier

if args.test == True:
    filename_ct = "/home/mescalin/bioinf/Johannes/viennaRNAfoldbpp/ASE_00089.ct"

### ================= import of ct files =================================== ###
counter = 0
g = open(filename_ct)
for x in g:
    if x[:1] == "#":
        counter += 1
    else:
        counter = counter + 1
        data_ct = pd.read_csv( open(filename_ct), sep="\s+", engine='python', keep_default_na = False, skiprows = counter, names = ['n', 'Base','n-1','n+1','bp','nn'])
        data_canonical_ct = data_ct
        data_ct = data_ct.drop(['n','n-1','n+1'], axis = 1)
        data_corrected_ct = data_ct
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
if args.canonical == True:
    unmod_sequence = rna_mod(sequence)

    for i in range(1, col_ct):
        if bpdata[i] > i:
            bp = (unmod_sequence[i - 1], unmod_sequence[bpdata[i] - 1])
            if bp not in [("G","C"), ("C","G"), ("A","U"), ("U","A"), ("G","U"), ("U","G")]:
                bpdata[bpdata[i]] = 0
                bpdata[i] = 0
if args.pkrem == True:
    bpdata = RNA.pt_pk_remove(bpdata)
    reference_db = RNA.db_from_ptable(bpdata)

### ============= open file, read until Distribution section begins, ======= ###
### =================  grab all lines and put it into list of lists ======== ###
data_drtr = list()
f = open("temp_data_mod_4.txt", "r")
for x in f:
    if x[:3] == "# D":
        lines = [line for line in f]
        lines = lines[1:len(lines)]
        for i in range(0, len(lines)):
            string = lines[i:i+1]
            data_drtr.append(string)

### ================= creation of bpp matrix =============================== ###
bpp_matrix_drtr = np.zeros((col_ct + 1, col_ct + 1))
for j in range(0,len(data_drtr)):
    single = data_drtr[j][0]
    splitter = single.split()
    pair_table_drtr = (splitter[2])
    probabilties_drtr = (splitter[4])
    pair_table_drtr = RNA.ptable_from_string(pair_table_drtr)
    matrix_len = len(pair_table_drtr)
    matrix_len = int(matrix_len)
    for i in range(1,matrix_len):
        if pair_table_drtr[i] > i and pair_table_drtr[i] != 0:
            probabilties_drtr = float(probabilties_drtr)
            bpp_matrix_drtr[i, pair_table_drtr[i]] += probabilties_drtr

### ================= db for drtr mfe ====================================== ###
probs_list = list()
for j in range(0,len(data_drtr)):
    single = data_drtr[j][0]
    splitter = single.split()
    probabilties_drtr = (splitter[4])
    probs_list.append(probabilties_drtr)
max_prob = float(max(probs_list))

for j in range(0, len(data_drtr)):
    single = data_drtr[j][0]
    splitter = single.split()
    pair_table_drtr = (splitter[2])
    probabilties_drtr = float((splitter[4]))
    if probabilties_drtr == max_prob:
        max_prob_drtr = pair_table_drtr
max_prob_drtr_t = RNA.ptable_from_string(max_prob_drtr)

mfe_data_drtr = data_drtr[0]
single2 = mfe_data_drtr[0]
splitter2 = single2.split()
drtr_mfe_pt_t = (splitter2[2])
drtr_mfe_pt = RNA.ptable_from_string(drtr_mfe_pt_t)

## ======================= call pps, ense_def,ens_dis ====================== ###
(tp, fp, tn, fn) = performance_bpp(bpdata, bpp_matrix_drtr)
pps = pair_prob_score(bpdata, bpp_matrix_drtr)
ense_dist = ensemble_distance(bpdata, bpp_matrix_drtr)
ense_def = ensemble_defect(bpdata, bpp_matrix_drtr)

(tp_mfe, fp_mfe, tn_mfe, fn_mfe) = performance_mfe(bpdata, drtr_mfe_pt)
(tp_prob, fp_prob, tn_prob, fn_prob) = performance_mfe(bpdata, max_prob_drtr_t)

### ================== verbose ============================================= ###
if args.verbose == True:
    # true_positives
    print(name)
    print(bpdata)
    print(data_drtr)
    # print(filename_ct)
    print(bpp_matrix_drtr)
    print("bpp_tp>")
    print(tp)

    # false_positives
    print("bpp_fp>")
    print(fp)

    print("bpp_tn>")
    print(tn)

    # false_negatives
    print("bpp_fn>")
    print(fn)
    print("")

### =================== results ============================================ ###
print("%s\t"
      "%i\t%i\t%i\t%i\t"
      "%i\t%i\t%i\t%i\t"
      "%8.5f\t%8.5f\t%8.5f\t%8.5f\t"
      "%8.5f\t"
      "%8.5f\t"
      "%8.5f\t"
      "%d" % (name,
            tp_mfe,
            fp_mfe,
            tn_mfe,
            fn_mfe,
            tp_prob,
            fp_prob,
            tn_prob,
            fn_prob,
            tp, #bpp
            fp, #bpp
            tn, #bpp
            fn, #bpp
            pps,
            ense_dist,
            ense_def,
            col_ct))
