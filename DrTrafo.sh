#!/bin/bash

#SBATCH --job-name=dr_tr_fasta_jm
#SBATCH --cpus-per-task=1
#SBATCH --cores-per-socket=1
#SBATCH --mem=8G
#SBATCH --time=0:0:0
#SBATCH --output=/home/output/DrTrafo_%A_%a.log


# temporary directory on cluster node
TMPDIR=home/scratch/${SLURM_JOB_USER}/${SLURM_JOB_ID}

# input and output directories
FILE_LIST="/home/input_fasta_files.txt"
OUTDIR="home/output"

# change into tmpdir
cd ${TMPDIR}

# LINE_NUMBER=${SLURM_ARRAY_TASK_ID}

if test "x$1" = "x"
then
    exit 1
else
    LINE_NUMBER=$1
fi

# task_ID defines which line(=filename) to pick from the input file
INPUT_FILE=`awk -v line_number=${LINE_NUMBER} 'NR==line_number {print}' ${FILE_LIST}`

file_prefix=`basename $INPUT_FILE .fasta`

cp ${INPUT_FILE} .

/scr/titan/bioinf/DrTrafo/scripts/RNAmod_convert.pl \
  -u \ # optional settings for RNAmod_convert.pl
  -m \ #
  ${file_prefix}.fasta \
| DrTransformer.py \
  --name $file_prefix \
  --logfile

gzip ${file_prefix}.log

chmod g+w ${file_prefix}.log.gz

mkdir -p ${OUTDIR}

# move results back home
mv ${file_prefix}.log.gz ${OUTDIR}/
