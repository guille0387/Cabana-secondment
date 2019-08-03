#!/usr/bin/env python3.7

import re
import argparse
import sys
import operator
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description = "Generate viral annotation tables for testing Federico's R script")
parser.add_argument("-f", "--fasta", dest = "proteins", help = "Fasta file containing proteins predicted with Prodigal", required = True)
parser.add_argument("-t", "--table", dest = "vphmm", help = "Table containing informative annotations", required = True)
if len(sys.argv) == 1:
	parser.print_help()
else:
	args = parser.parse_args()
	protein_file = args.proteins
	annotation_table = args.vphmm
	vphmm_hit_df = pd.read_table(annotation_table)
	dataset_id = protein_file.split("_")[0]
	annotation_list = []
	for protein in SeqIO.parse(protein_file, "fasta"):
		contig_id = re.search(r"[ER0-9]+_\d+[-Pro]*", protein.id).group(0)
		protein_prop = protein.description.split(" # ")
		del(protein_prop[4])
		if protein_prop[0] in vphmm_hit_df["query"].values:
			filtered_df = vphmm_hit_df[vphmm_hit_df["query"] == protein_prop[0]]
			if len(filtered_df) > 1:
				best_value_index = max(filtered_df["Abs_Evalue_exp"].items(), key = operator.itemgetter(1))[0]
				protein_prop.extend(list(filtered_df.loc[best_value_index, ["ViPhOG", "Abs_Evalue_exp", "Taxon"]]))
			else:
				protein_prop.extend(list(filtered_df.loc[filtered_df.index[0], ["ViPhOG", "Abs_Evalue_exp", "Taxon"]]))
		else:
			protein_prop.extend(["No hit", "NA", ""])
		annotation_list.append([contig_id] + protein_prop)
	final_map_df = pd.DataFrame(annotation_list, columns = ["Contig", "CDS_ID", "Start", "End", "Direction", "Best_hit", "Abs_Evalue_exp", "Label"])
	final_map_df.to_csv("%s_annotations_contig_map.tsv" % dataset_id, sep = "\t", index = False)
