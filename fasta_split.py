#!/usr/bin/env python2

import argparse
import random
import sys
from Bio import SeqIO

# define arguments
parser = argparse.ArgumentParser(description='Split FASTA into random chunks')
parser.add_argument('-f', dest='fasta', help='FASTA file to split')
parser.add_argument('-p', dest='prefix', help='Prefix name for output file')
parser.add_argument('-min', dest='min_contig', help='Minimum contig length')
parser.add_argument('-max', dest='max_contig', help='Maximum contig length')

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)
else:
    args = parser.parse_args()

seqs = {}

# retrieve sequences from file
with open(args.fasta, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        seqs[record.id] = record.seq

# set interval and output name
prefix = args.prefix
min_len = int(args.min_contig)
contig_n = 0

# go through each sequence
for s in seqs.keys():
    present_length = 0
    total_length = len(seqs[s])
    max_len = int(args.max_contig)
    # while the sequence added is below the total sequence length
    while present_length < total_length:
        contig_n += 1
        # check if max_len is not above the remaining sequence length
        if max_len > (total_length - present_length):
            max_len = total_length - present_length
        try:
            dist = random.sample(range(min_len,max_len),1)
            print ">%s_%i" % (prefix, contig_n)
            print "%s" % (seqs[s][present_length:present_length+dist[0]])
            present_length += dist[0]
        # if max_length is below the set min_length just ignore that contig
        except:
            break
