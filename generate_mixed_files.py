#!/usr/bin/env python3

import os
import glob
import argparse
import sys
import re
import random
import numpy as np
from Bio import SeqIO

def comb_seq_files(primary_file, secondary_file, total_seqs_num, prop_range, output_dir):  
	for num in np.arange(prop_range[0], sum(prop_range[1:]), prop_range[-1]):
		primary_seqs_num = int(total_seqs_num * num)
		target_records = []
		for record in SeqIO.parse(primary_file, "fasta"):
			if random.randint(0, 1):
				target_records.append(record)
				if len(target_records) == primary_seqs_num:
					break
		for record in SeqIO.parse(secondary_file, "fasta"):
			if random.randint(0, 1):
				target_records.append(record)
				if len(target_records) == total_seqs_num:
					break
		random.shuffle(target_records)
		SeqIO.write(target_records, os.path.join(output_dir, "%s_percent_primary.fna" % int(num * 100)), "fasta")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Generate combined fasta file with user-defined proportions of original fasta files")
	parser.add_argument("-p", dest="primary", help="Relative or absolute path to primary fasta file", required=True)
	parser.add_argument("-s", dest="secondary", help="Relative or absolute path to secondary fasta file", required=True)
	parser.add_argument("-t", dest="total", type=int, help="Total number of seqs for combined fasta file", required=True)
	parser.add_argument("-r", dest="prop", help="Range of proportions for primary fasta file - format: lowest_value,highest_value,increment", required=True)
	parser.add_argument("-o", dest = "outdir", help = "Relative path to directory where you want the output files to be stored (default: cwd)", default = ".")
	if len(sys.argv) == 1 :
		parser.print_help()
		sys.exit(1)
	else:
		args = parser.parse_args()
		primary_prop = [float(x) for x in args.prop.split(",")]
		comb_seq_files(args.primary, args.secondary, args.total, primary_prop, args.outdir)
