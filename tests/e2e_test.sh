#!/bin/bash

# The current implementation of BEstimate.py accesses
# some files that are assumed to be in the same folder.
# This forces Bestimate.py to be called from within the
# BEstimate/ folder and therefore we have to run this
# script as:
#   bash ../tests/e2e_test.sh

# Prepare output folder
mkdir -p tests/output

# Run Bestimate
python3 BEstimate.py -gene SRY \
                     -assembly GRCh38 \
                     -pamseq NGG \
                     -pamwin 21-23 \
                     -actwin 4-8 \
                     -protolen 20 \
                     -edit C \
                     -edit_to T \
                     -o tests/output/ \
                     -ofile SRY_CBE_NGG

# check exit status
status=$?
if [ $status -ne 0 ];
then
    printf 'BEstimate failed with status "%s"' "$status";
    exit -1
fi;

# compare output with expected
file1="tests/output/SRY_CBE_NGG_edit_df.csv"
file2="../tests/data/SRY_CBE_NGG_edit_df_expected.csv"

results_dif=$( diff $file1 $file2 | wc -l )
if [ "$results_dif" -ne 0 ]; then
    printf 'The file "%s" is different from "%s"\n' "$file1" "$file2"
    exit -1
fi

rm -rf tests/output
printf 'DONE!'
