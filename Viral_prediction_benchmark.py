#!/usr/bin/env python3.7

import os
import sys
import re
import glob
import subprocess
import argparse
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, "/homes/garp/scripts")

from Host_cleaning_process import contig_clean

parser = argparse.ArgumentParser(description = "Calculate prediction summary metrics for method based on Sullivan 2019 and append to VirFinder methods table")
parser.add_argument("-i", "--ident", dest = "ident", help = "Accession number of run to be analysed", required = True)
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
else:
	args = parser.parse_args()
	dataset_id = args.ident
	if re.search(r"^ERR", dataset_id):
		project_id = "ERP006614"
	else:
		project_id = "SRP074090"
	blast_file = glob.glob("tests/*/megablast/%s*500bp.tbl" % dataset_id)[0]
	blast_viral_contigs = contig_clean(blast_file)
	virsorter_high = [x for x in glob.glob("viral_annotation_pipeline/VirSorter_500bp_contigs/%s_virsorter/Predicted_viral_sequences/*" % dataset_id) if re.search(r"cat-[12].fasta$", x)]
	predicted_viruses = []
	for item in virsorter_high:
		if os.stat(item).st_size > 0:
			for record in SeqIO.parse(item, "fasta"):
				contig_name_search = re.search(r"(NODE\w+cov_\d+)_(\d+)[a-z-]+(cat_\d)", record.name)
				contig_name_string = "%s.%s VirSorter_%s" % (contig_name_search.group(1), contig_name_search.group(2), contig_name_search.group(3))
				predicted_viruses.append(contig_name_string)
	virfinder_file = glob.glob("viral_annotation_pipeline/VirFinder_test/finder_500bpfilter/%s*" % dataset_id)[0]
	virfinder_df = pd.read_csv(virfinder_file, sep = "\t")
	virfinder_high = list(virfinder_df[(virfinder_df["pvalue"] < 0.05) & (virfinder_df["score"] >= 0.9)]["name"].values)
	if len(virfinder_high) > 0 and len(predicted_viruses) > 0:
		virsorter_high_viruses = [x.split(" ")[0] for x in predicted_viruses]
		for i,j in enumerate(predicted_viruses):
			if j.split(" ")[0] in virfinder_high:
				predicted_viruses[i] = j + "_VirFinder_score_0.9"
		predicted_viruses.extend(["%s VirFinder_score_0.9" % x for x in virfinder_high if x not in virsorter_high_viruses])
	elif len(virfinder_high) > 0:
		predicted_viruses.extend(["%s VirFinder_score_0.9" % x for x in virfinder_high])
	virsorter_low = glob.glob("viral_annotation_pipeline/VirSorter_500bp_contigs/%s_virsorter/Predicted_viral_sequences/*cat-3.fasta" % dataset_id)[0]
	if os.stat(virsorter_low).st_size > 0:
		virsorter_low_viruses = []
		for record in SeqIO.parse(virsorter_low, "fasta"):
			contig_name_search = re.search(r"(NODE\w+cov_\d+)_(\d+)", record.name)
			contig_name_string = "%s.%s" % (contig_name_search.group(1), contig_name_search.group(2))
			virsorter_low_viruses.append(contig_name_string)
		virfinder_low = list(virfinder_df[(virfinder_df["pvalue"] < 0.05) & (virfinder_df["score"] >= 0.7) & (virfinder_df["score"] < 0.9)]["name"].values)
		if len(virfinder_low) > 0:
			predicted_viruses.extend(["%s VirSorter_cat_3_VirFinder_score_0.7" % x for x in list(set(virsorter_low_viruses) & set(virfinder_low))])
	viral_prediction_result = [x.split(" ")[0] for x in predicted_viruses]
	TP = len(list(set(viral_prediction_result) & set(blast_viral_contigs)))
	FP = len([x for x in viral_prediction_result if x not in blast_viral_contigs])
	FN = len([x for x in blast_viral_contigs if x not in viral_prediction_result])
	all_row = pd.DataFrame([(dataset_id, "Combined_all", TP, FP, FN, project_id)], columns = ["Dataset", "Method", "TP", "FP", "FN", "Project"])
	viral_prediction_result = [x.split(" ")[0] for x in predicted_viruses if "VirSorter" in x]
	TP = len(list(set(viral_prediction_result) & set(blast_viral_contigs)))
	FP = len([x for x in viral_prediction_result if x not in blast_viral_contigs])
	FN = len([x for x in blast_viral_contigs if x not in viral_prediction_result])
	no_virfinder_row = pd.DataFrame([(dataset_id, "Combined_no_virfinder", TP, FP, FN, project_id)], columns = ["Dataset", "Method", "TP", "FP", "FN", "Project"])
	result_df = pd.read_csv(glob.glob("viral_annotation_pipeline/VirFinder_test/VirFinder_method_metrics/%s*" % dataset_id)[0], sep = "\t")
	final_df = pd.concat([result_df, all_row, no_virfinder_row])
	final_df.to_csv("viral_annotation_pipeline/VirFinder_test/VirFinder_method_metrics/%s_combined.tsv" % dataset_id, sep = "\t", index = False)
	
	
		
