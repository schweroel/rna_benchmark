#!/bin/bash

# content: bash script that takes every fasta file of a defined directory and
# sends it to fold.py, which in turn calculates mfe secondary structure
# and bpp secondary structure; the results are saved in a .pickle file, that
# is named after the RNA ID

for file in /home/*".fasta"
do
  ID=${file##*/}
  ID=${ID%.fasta}
  direct="/home/"
  echo "$file" | python3 fold.py -i "$ID" -d "$direct"
  echo $ID
done
