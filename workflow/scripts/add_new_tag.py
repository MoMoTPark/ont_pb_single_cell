#!/usr/bin/env python

# Author: Moe
# Add a new BAM tag to an input mapped BAM file

import pysam
import sys

print(sys.executable)
in_bam = sys.argv[1]
out_bam = sys.argv[2]

# @RG header
rg = {
    "ID": "e4927d21",
    "PL": "PACBIO",
    "DS": "READTYPE=TRANSCRIPT",
    "PU": "transcript",
    "SM": "UnnamedSample",
    "PM": "SEQUEL"
}

with pysam.AlignmentFile(in_bam, 'rb') as alignment:
    print(alignment)

    header = alignment.header.to_dict()
    if 'RG' not in header:
        header['RG'] = []
    header['RG'].append(rg)

    with pysam.AlignmentFile(out_bam, 'wb', header=header) as output:
        for read in alignment:
            read.query_qualities = None
            read.set_tag('im', read.query_name, value_type='Z')
            read.set_tag('is', 1, value_type='i')
            read.set_tag('rc', 1, value_type='i')
            read.set_tag('RG', "e4927d21", value_type='Z')
            # read.set_tag('mg', 99.8506, value_type='f')
            read.set_tag('XA', "XM-CB", value_type='Z')
            read.set_tag('XM', read.get_tag("UB"), value_type='Z')
            read.set_tag('ic', 1, value_type='i')
            read.set_tag('NM', 1, value_type='i')
            output.write(read)