#!/usr/bin/env python3.7

import os
import sys
import re
import glob
import argparse
import pandas as pd

sys.path.insert(0, "/homes/garp/scripts/")

from Host_cleaning_process import contig_clean

parser = argparse.ArgumentParser(description = "Generate table reporting metrics for different VirFinder-based prediction methods, employing an user-defined fdr threshold")
parser.add_argument("-i", "--input", dest = "input_file", help = "Relative or absolute path to input table file that collates VirFinder prediction results for different methods", required = True)
parser.add_argument("-o", "--outdir", dest = "outdir", help = "Relative or absolute path to directory where output is to be stored (default: cwd)", default = ".")
parser.add_argument("-t", "--thres", dest = "thres", help = "Fdr threshold for filtering out non-significant viral predictions (default: 0.1)", type = float, default = 0.1 )
if len(sys.argv) == 1 :
	parser.print_help()
	sys.exit(1)
else:
	args = parser.parse_args()
	input_file = args.input_file
	output_dir = args.outdir
	fdr_thres = args.thres
	dataset_id = os.path.basename(input_file).split("_")[0]
	blast_files = []
	extensions = ["*500bp.tbl", "*1kb.tbl"]
	for item in extensions:
		blast_files.extend(glob.glob("tests/*/megablast/%s%s" % (dataset_id, item)))
	print(blast_files)
	dataframe_rows = []
	viral_contigs_500 = contig_clean(blast_files[0])
	viral_contigs_1000 = contig_clean(blast_files[1])
	virfinder_df = pd.read_csv(input_file, sep = "\t")
	significant_df = virfinder_df[virfinder_df["fdr"] < fdr_thres]
	if len(significant_df) < 1:
		print("None of the processed contigs had an associated fdr value below %s" % fdr_thres)
		sys.exit(0)
	method_list = list(significant_df["Method"].value_counts().index)
	for method in method_list:
		method_df = significant_df[significant_df["Method"] == method]
		viral_predictions = list(method_df["name"].values)
		if "500bp" in method:
			TP = len(list(set(viral_predictions) & set(viral_contigs_500)))
			FP = len([x for x in viral_predictions if x not in viral_contigs_500])
			FN = len([x for x in viral_contigs_500 if x not in viral_predictions])
			dataframe_rows.append((dataset_id, method, TP, FP, FN))
		else:
			TP = len(list(set(viral_predictions) & set(viral_contigs_1000)))
			FP = len([x for x in viral_predictions if x not in viral_contigs_1000])
			FN = len([x for x in viral_contigs_1000 if x not in viral_predictions])
			dataframe_rows.append((dataset_id, method, TP, FP, FN))
	final_df = pd.DataFrame(dataframe_rows, columns = ["Dataset", "Method", "TP", "FP", "FN"])
	final_df.to_csv(os.path.join(output_dir, "%s_VirFinder_method_metrics_%s_fdr.tsv" % (dataset_id, fdr_thres)), sep = "\t", index = False)
