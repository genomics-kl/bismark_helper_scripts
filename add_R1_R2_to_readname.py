#!/usr/bin/env python

import pysam
import argparse
import os.path
import sys

parser = argparse.ArgumentParser(description='Add R1 and R2 to BAM \
        query names.')
parser.add_argument("in_bam",
                    help="Input BAM.")
parser.add_argument("out_bam",
                    help="Output BAM.")
args = parser.parse_args()

in_bam = args.in_bam
out_bam = args.out_bam

# file tests
assert os.path.isfile(in_bam), "Input file does not exist."
assert not os.path.isfile(out_bam), "Output file exists."

# read BAM
bamfile = pysam.AlignmentFile(in_bam, "rb")

renamed_bam = pysam.AlignmentFile(out_bam, "wb", template=bamfile)

for read in bamfile.fetch():
    if read.is_read1:
        read.query_name = read.query_name + "_R1"
    elif read.is_read2:
        read.query_name = read.query_name + "_R2"
    else:
        sys.exit("Read is neither R1 nor R2.")

    renamed_bam.write(read)

renamed_bam.close()
bamfile.close()
