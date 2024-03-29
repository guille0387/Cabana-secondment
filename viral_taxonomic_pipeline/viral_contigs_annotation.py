#!/usr/bin/env python3.7

import os
import re
import sys
import argparse
import operator
import pandas as pd
from Bio import SeqIO


def prot_annot_tbl(protein_file, ratio_evalue_file):
	'''This function takes a fasta file containing the proteins predicted in a set of putative viral contigs and a dataframe that collates the
	   results obtained with hmmscan against the ViPhOG database for the same proteins'''
	annotation_list = []
	ratio_evalue_df = pd.read_csv(ratio_evalue_file, sep = "\t")
	for protein in SeqIO.parse(protein_file, "fasta"):
		contig_id = re.split(r"_\d+$", protein.id)[0]
		protein_prop = protein.description.split(" # ")[:-1]
		if protein_prop[0] in ratio_evalue_df["query"].values:
			filtered_df = ratio_evalue_df[ratio_evalue_df["query"] == protein_prop[0]]
			if len(filtered_df) > 1:
				best_value_index = max(filtered_df["Abs_Evalue_exp"].items(), key = operator.itemgetter(1))[0]
				protein_prop.extend(list(filtered_df.loc[best_value_index, ["ViPhOG", "Abs_Evalue_exp", "Taxon"]]))
			else:
				protein_prop.extend(list(filtered_df.loc[filtered_df.index[0], ["ViPhOG", "Abs_Evalue_exp", "Taxon"]]))
		else:
			protein_prop.extend(["No hit", "NA", ""])
		annotation_list.append([contig_id] + protein_prop)
	protein_annot_df = pd.DataFrame(annotation_list, columns = ["Contig", "CDS_ID", "Start", "End", "Direction", "Best_hit", "Abs_Evalue_exp", "Label"])
	return protein_annot_df

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Generate tabular file with ViPhOG annotation results for proteins predicted in viral contigs")
	parser.add_argument("-p", "--prot", dest = "prot_file", help = "Relative or absolute path to protein file of predicted viral contigs", required = True)
	parser.add_argument("-t", "--table", dest = "ratio_file", help = "Relative or absolute path to ratio_evalue tabular file generated for predicted viral contigs", required = True)
	parser.add_argument("-o", "--outdir", dest = "output_dir", help = "Relative path to directory where you want the output file to be stored (default: cwd)", default = ".")

	if len(sys.argv) == 1:
		parser.print_help()

	else:
		args = parser.parse_args()
		final_df = prot_annot_tbl(args.prot_file, args.ratio_file)
		final_df.to_csv(os.path.join(args.output_dir, re.split(r"\.[a-z]+$", os.path.basename(args.prot_file))[0] + "_ViPhOG_annot.tsv"), sep = "\t", index = False)
