#!/usr/bin/env python3.7

import os
import pandas as pd
import glob
import sys
import argparse
import re
from Bio import SeqIO

def contig_clean(blast_file):
	contigs_remove = []
	blast_df = pd.read_csv(blast_file, sep = "\t")
	contig_list = list(blast_df["qseqid"].value_counts().index)
	for contig in contig_list:
		#print(contig)
		filtered_df = blast_df[blast_df["qseqid"] == contig]
		filtered_df = filtered_df.reset_index()
		contig_length = filtered_df["qlen"][0]
		if len(filtered_df) > 1:
			#print("More than one alignment")
			tuple_list = [(filtered_df.qstart[i], filtered_df.qend[i]) for i in range(len(filtered_df))]
			#print(tuple_list)
			sorted_tuple_list = sorted([tuple(sorted(t)) for t in tuple_list])
			#print(sorted_tuple_list)
			reference = list(sorted_tuple_list[0])
			final_tuple_list = []
			for i,j in sorted_tuple_list[1:]:
				if i <= reference[1]:
					reference[1] = max(j, reference[1])
				else:
					final_tuple_list.append(tuple(reference))
					reference[0] = i
					reference[1] = j
			final_tuple_list.append(tuple(reference))
			#print(final_tuple_list)
			total_align = sum([j - i + 1 for i,j in final_tuple_list])
			if (total_align/contig_length) >= 0.6 or (contig_length - total_align) < 200:
				contigs_remove.append(contig)
		else:
			#print("Only one alignment")
			align_coords = sorted((filtered_df.qstart[0], filtered_df.qend[0]))
			#print(align_coords)
			total_align = align_coords[1] - align_coords[0] + 1
			if (total_align/contig_length) >= 0.6 or (contig_length - total_align) < 200:
				contigs_remove.append(contig)
	return contigs_remove


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Remove host-related contigs based on blast output and write filtered contig file")
	parser.add_argument("-f", "--fasta", dest = "input_file", help = "Input fasta file", required = True)
	parser.add_argument("-b", "--blastout", dest = "blast", help = "Blast output file", required = True)
	parser.add_argument("-o", "--outdir", dest = "outdir", help = "Output directory", default = ".")

	if len(sys.argv) == 1:
        	parser.print_help()
	else:
		args = parser.parse_args()
		input_file = args.input_file
		blast_file = args.blast
		output_dir = args.outdir
		list_to_remove = contig_clean(blast_file)
		output_name = os.path.join(output_dir, re.split(r"\.[a-z]+$", os.path.basename(input_file))[0] + "_host_filtered.fasta")
		contig_generator = SeqIO.parse(input_file, "fasta")
		records_to_keep = []
		for record in contig_generator:
			if record.name not in list_to_remove:
				records_to_keep.append(record)
		SeqIO.write(records_to_keep, output_name, "fasta")
