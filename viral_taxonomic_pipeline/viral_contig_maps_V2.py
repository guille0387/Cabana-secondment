#!/usr/bin/env python3.7

import os
import sys
import re
import argparse
import operator
import subprocess
import pandas as pd
from Bio import SeqIO

sys.path.extend(["/homes/garp/scripts", "/homes/garp/scripts/viral_taxonomic_pipeline"])

from filter_contigs_len import filter_contigs
from write_viral_predictions_V2 import virus_pred
from Generate_vphmm_hmmer_matrix import hmmer_domtbl
from Ratio_Evalue_table import ratio_evalue

def prot_annot_tbl(protein_file, ratio_evalue_df):
	'''This function takes a fasta file containing the proteins predicted in a set of putative viral contigs and a dataframe that collates the
	   results obtained with hmmscan against the ViPhOG database for the same proteins'''  
	annotation_list = []
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


def process_viral_files(viral_pred_file, outdir, num_thread):
	'''This function takes a fasta file containing viral predictions and carries out the taxonomic annotation analysis for each viral contig'''
	#Predict proteins encoded in viral sequences and compare them to the ViPhOG database, done by hmmer_domtbl function
	hmmer_domtbl(viral_pred_file, outdir, num_thread)
	prodigal_output_file = os.path.join(outdir, re.split(r"\.[a-z]+$", os.path.basename(viral_pred_file))[0] + "_CDS.faa")
	hmmer_output_file = os.path.join(outdir, re.split(r"\.[a-z]+$", os.path.basename(viral_pred_file))[0] + "_hmmer_ViPhOG.tbl")
	ViPhOG_overall_df = pd.read_csv(hmmer_output_file, sep = "\t")
	#Generate a new table containing only informative hits, model alingment ratios and absolute values of total Evalue exponents
	informative_df = ratio_evalue(ViPhOG_overall_df)
	if isinstance(informative_df, str):
		print("No informative hits against the ViPhOG database were found for: %s" % os.path.basename(viral_pred_file))
	else:
		informative_output = hmmer_output_file.rstrip("_hmmer_ViPhOG.tbl") + "_ViPhOG_informative.tsv"
		informative_df.to_csv(informative_output, sep = "\t", index = False)
		#Generate annotation table with data for each protein per row and use this table for creating gene maps for each viral sequence
		contig_map_df = prot_annot_tbl(prodigal_output_file, informative_df)
		contig_map_dir = os.path.join(outdir, re.split(r"\.[a-z]+$", os.path.basename(viral_pred_file))[0] + "_maps")
		os.mkdir(contig_map_dir)
		contig_map_table = prodigal_output_file.rstrip(".faa") + "_ViPhOG_annot.tsv"
		contig_map_df.to_csv(contig_map_table, sep = "\t", index = False)
		subprocess.call("Make_viral_contig_map.R -t %s -o %s" % (contig_map_table, contig_map_dir), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)
		print("Gene maps generation completed for %s" % os.path.basename(viral_pred_file))


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description = "Generate gene maps with ViPhOG annotations of putative viral contigs and prophage elements")
	parser.add_argument("-f", "--fasta", dest = "input_file", help = "Relative or absolute path of metagenomic assembly fasta file", required = True)
	parser.add_argument("-i", "--identifier", dest = "ident", help = "Identifier or accession number associated with your dataset, used for naming the pipeline's output", required = True)
	parser.add_argument("-o", "--outdir", dest = "outdir", help = "Relative path to directory where you want the output directory to be stored (default: cwd)", default = ".")
	parser.add_argument("-l", "--length", dest = "length", type = float, help = "Length threshold in kb for filtering contigs stored in the assembly file (default: 0.5)", default = "0.5")
	parser.add_argument("-t", "--thread", dest = "thread", type = int, help = "Number of threads to use for the computation (default: 4)", default = "4")
	parser.add_argument("--virome", dest = "virome", action = "store_true", help = "Use flag to indicate that your dataset was derived from a viral-enriched sample")
	parser.add_argument("--prok", dest = "mode", action = "store_true", help = "Use flag to select viral prediction model trained only with prokaryotic viral sequences")

	if len(sys.argv) == 1:
		parser.print_help()
	else:
		args = parser.parse_args()
		input_file = args.input_file
		dataset_id = args.ident
		length_thres = args.length
		num_thread = args.thread
		if args.outdir == ".":
			output_dir = os.path.join(os.getcwd(), "%s_viral_prediction" % dataset_id)
		else:
			output_dir = os.path.join(args.outdir, "%s_viral_prediction" % dataset_id)
		os.mkdir(output_dir)
		#Filter the metagenome assembly and retain only those contigs that are at least 1000 bp long
		filter_contigs(contig_file = input_file, thres = length_thres, output_dir = os.path.dirname(input_file), run_id = dataset_id)
		filtered_file_name = re.split(r"\.[a-z]+$", input_file)[0] + "_filt%sbp.fasta" % int(length_thres * 1000)
		if os.stat(filtered_file_name).st_size == 0:
			print("None of the assembled contigs is at least %s kb long" % length_thres)
			os.remove(filtered_file_name)
			sys.exit(0)
		#Detect viral contigs and prophages using function virus_pred (VirFinder + VirSorter)
		predicted_viruses = virus_pred(filtered_file_name, output_dir, num_thread, args.virome, args.mode)
		if sum([len(x) for x in predicted_viruses]) == 0:
			print("Overall, no putative viral contigs or prophages were retrieved for the analysed metagenomic assembly")
			sys.exit(0)
		else:
			viral_prediction_files = []
			if len(predicted_viruses[0]) > 0:
				HC_prediction_file = os.path.join(output_dir, "High_confidence_putative_viral_contigs.fna")
				SeqIO.write(predicted_viruses[0], HC_prediction_file, "fasta")
				viral_prediction_files.append(HC_prediction_file)
			if len(predicted_viruses[1]) > 0:
				LC_prediction_file = os.path.join(output_dir, "Low_confidence_putative_viral_contigs.fna")
				SeqIO.write(predicted_viruses[1], LC_prediction_file, "fasta")
				viral_prediction_files.append(LC_prediction_file)
			if len(predicted_viruses[2]) > 0:
				Prophage_prediction_file = os.path.join(output_dir, "Putative_prophages.fna")
				SeqIO.write(predicted_viruses[2], Prophage_prediction_file, "fasta")
				viral_prediction_files.append(Prophage_prediction_file)
			for item in viral_prediction_files:
				process_viral_files(item, output_dir, num_thread)
