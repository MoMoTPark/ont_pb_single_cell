#!/usr/bin/env python


# Author: Moe
# Modify SQANTI3 filter output to make it compatible as a classfication input file for pigeon make-seurat tool (<wld>_classification.filtered_lite_classification.txt); add two additional columns ('fl_assoc' and 'cell_barcodes') to <wld>_RulesFilter_result_classification.txt and remove 'filter_result' column.
# Input: SQ3 filtered classification output, Pigeon lite classification output, output file path

import sys
import pandas as pd

print(sys.executable)

# SQ3 filtered classification output file
sqanti = sys.argv[1]
# Pigeon full dataset classification output file
pigeon = sys.argv[2]
# Path for output file
outfile = sys.argv[3]

print(sqanti)
print(pigeon)
print(outfile)

with open(sqanti, 'r') as f:
    sq_df = pd.read_csv(f, sep="\t")

print(sq_df.head())

with open(pigeon, 'r') as f:
    pg_df = pd.read_csv(f, sep="\t")

print(pg_df.head())

left = sq_df[sq_df['filter_result'] == 'Isoform'].drop(labels='filter_result', axis=1)
right = pg_df[['isoform', 'fl_assoc', 'cell_barcodes']]
# All values in left should be present in right dataframe
assert left.shape[0] <= right.shape[0], "\n\t*** ERROR: SQ3 dataframe has more rows than Pigeons classification dataframe ***"
out_df = pd.merge(left, right, how='left', on='isoform').reset_index(drop=True)

with open(outfile, 'w') as f:
    out_df.to_csv(f, sep='\t', index=False, na_rep='NA')