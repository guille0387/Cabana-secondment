#!/usr/bin/env python3.7

import os
import subprocess
import re
import argparse
import sys
import glob
import pandas as pd
from Bio import SeqIO

def virus_pred(assembly_file, output_dir, num_thread, virome_dataset, prok_mode):

	"""Run VirFinder and VirSorter on provided metagenomic assembly file and parse the output to generate lists of viral predictions sorted in different significance categories, along with prophage predictions"""

	# VirFinder run using either the default or general viral prediction model according to user's instructions

	if prok_mode:
		print("Running VirFinder only for prokaryotic viruses")
		subprocess.call("VirFinder_analysis_Prok.R -f %s -o %s" % (assembly_file, output_dir), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True) 
	else:
		print("Running VirFinder for prokaryotic and eukaryotic viruses")
		subprocess.call("VirFinder_analysis_Euk.R -f %s -o %s" % (assembly_file, output_dir), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)

	# Path to directory containing VirSorter's databases

	virsorter_data_path = "/hps/nobackup2/production/metagenomics/garp/databases/virsorter-data"

	# Run VirSorter using the RefSeq + Virome database and the decontamination mode for viral enriched samples if virome_dataset = True, otherwise run VirSorter without desontamination mode

	if virome_dataset:
		print("Running VirSorter with virome decontamination mode")
		subprocess.call("wrapper_phage_contigs_sorter_iPlant.pl -f %s --db 2 --wdir %s --ncpu %s --virome --data-dir %s" % (assembly_file, output_dir, num_thread, virsorter_data_path), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)
	else:
		print("Running VirSorter without virome decontamination mode")
		subprocess.call("wrapper_phage_contigs_sorter_iPlant.pl -f %s --db 2 --wdir %s --ncpu %s --data-dir %s" % (assembly_file, output_dir, num_thread, virsorter_data_path), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)

	# Define lists that will store putative viral records assigned to different categories

	HC_viral_predictions = []
	LC_viral_predictions = []
	prophage_predictions = []

	# Store VirSorter prediction files in lists based on pre-defined categories. Retrieve results from VirFinder and sort in high (p < 0.05 and score >= 0.9) and low (p < 0.05 and 0.07 <= score < 0.9) categories

	VirFinder_pred = glob.glob(os.path.join(output_dir, "*VirFinder_table-all.tab"))[0]
	VF_result_df = pd.read_csv(VirFinder_pred, sep = "\t")
	VF_high_ids = list(VF_result_df[(VF_result_df["pvalue"] < 0.05) & (VF_result_df["score"] >= 0.90)]["name"].values)
	if len(VF_high_ids) < 1:
		print("No contigs with p < 0.05 and score >= 0.90 were reported by VirFinder")
	VF_low_ids = list(VF_result_df[(VF_result_df["pvalue"] < 0.05) & (VF_result_df["score"] >= 0.70) & (VF_result_df["score"] < 0.9)]["name"].values)
	if len(VF_low_ids) < 1:
		print("No contigs with p < 0.05 and 0.70 <= score < 0.90 were reported by VirFinder")

	if len(os.listdir(os.path.join(output_dir, "Predicted_viral_sequences"))) > 0:
		VirSorter_viral_high = [x for x in glob.glob(os.path.join(output_dir, "Predicted_viral_sequences", "*.fasta")) if re.search(r"cat-[12]\.fasta$", x)]
		VirSorter_viral_low = glob.glob(os.path.join(output_dir, "Predicted_viral_sequences", "*cat-3.fasta"))[0]
		VirSorter_prophages = [x for x in glob.glob(os.path.join(output_dir, "Predicted_viral_sequences", "*.fasta")) if re.search(r"cat-[45]\.fasta$", x)]

		# Retrieve contig IDs of VirSorter predictions in categories 1 and 2

		print("Beginning the identification of high confidence viral predictions...")

		VS_high_tuples = []
		search_id = re.compile(r"VIRSorter_(\w+)-([a-z-]*cat_\d)")
		for item in VirSorter_viral_high:
			if os.stat(item).st_size != 0:
				with open(item) as input_file:
					for line in input_file:
						if search_id.search(line):
							VS_high_tuples.append((search_id.search(line).group(1), search_id.search(line).group(2)))


		# Store high confidence viral records in list, these correspond to viral predictions in VirSorter categories 1 and 2

		if len(VS_high_tuples) > 0:
			for x,y in VS_high_tuples: 
				for record in SeqIO.parse(assembly_file, "fasta"):
					contig_id_search = re.sub(r"[.,:; ]", "_", record.description)
					if contig_id_search == x and record.description in VF_high_ids:
						if "circular" in y:
							suff = "11_H_%s_circular" % y.split("_")[1]
						else:
							suff = "11_H_%s" % y.split("_")[1]
						record.id = "%s_%s" % (record.id, suff)
						record.description = "_".join(record.description.split()[1:])
						HC_viral_predictions.append(record)
						break
					elif contig_id_search == x and record.description in VF_low_ids:
						if "circular" in y:
							suff = "11_L_%s_circular" % y.split("_")[1]
						else:
							suff = "11_L_%s" % y.split("_")[1]
						record.id = "%s_%s" % (record.id, suff)
						record.description = "_".join(record.description.split()[1:])
						HC_viral_predictions.append(record)
						break
					elif contig_id_search == x:
						if "circular" in y:
							suff = "01_%s_circular" % y.split("_")[1]
						else:
							suff = "01_%s" % y.split("_")[1]
						record.id = "%s_%s" % (record.id, suff)
						record.description = "_".join(record.description.split()[1:])
						HC_viral_predictions.append(record)
						break

    		# Retain VirFinder predictions that were not included in VirSorter categories 1 and 2

			VS_high_ids = [x[0] for x in VS_high_tuples]
			VF_high_ids = [x for x in VF_high_ids if re.sub(r"[.,:; ]", "_", x) not in VS_high_ids]
			VF_low_ids = [x for x in VF_low_ids if re.sub(r"[.,:; ]", "_", x) not in VS_high_ids]
			print("High confidence viral predictions identified")

		else:
			print("No contigs were retrieved from VirSorter categories 1 and 2, therefore no high confidence viral contigs were reported")

		# Retrieve contig IDs of VirSorter predictions category 3

		print("Beginning the identification of low confidence viral predictions...")

		VS_low_tuples = []
		if os.stat(VirSorter_viral_low).st_size != 0:
			with open(VirSorter_viral_low) as input_file:
				for line in input_file:
					if search_id.search(line):
						VS_low_tuples.append((search_id.search(line).group(1), search_id.search(line).group(2)))

		# Store low confidence viral records in list

		if len(VS_low_tuples) > 0:
			id_suffix_list = []
			if len(VF_high_ids) > 0:
				for item in VF_high_ids:
					for x,y in VS_low_tuples:
						if re.sub(r"[.,:; ]", "_", item) == x:
							if "circular" in y:
								id_suffix_list.append((item, "11_H_3_circular"))
							else:
								id_suffix_list.append((item, "11_H_3"))
							break
					if item not in [elem[0] for elem in id_suffix_list]:
						id_suffix_list.append((item, "10_H"))
			if len(VF_low_ids) > 0:
				for item in VF_low_ids:
					for x,y in VS_low_tuples:
						if re.sub(r"[.,:; ]", "_", item) == x:
							if "circular" in y:
								id_suffix_list.append((item, "11_L_3_circular"))
							else:
								id_suffix_list.append((item, "11_L_3"))
							break
			if len(id_suffix_list) > 0:
				for x,y in id_suffix_list:
					for record in SeqIO.parse(assembly_file, "fasta"):
						if record.description == x:
							record.id = "%s_%s" % (record.id, y)
							record.description = "_".join(record.description.split()[1:])
							LC_viral_predictions.append(record)
							break
                
		else:
			print("No contigs were reported as VirSorter category 3")
			if len(VF_high_ids) > 0:
				for item in VF_high_ids:
					for record in SeqIO.parse(assembly_file, "fasta"):
						if record.description == item:
							record.id = record.id + "_10_H"
							record.description = "_".join(record.description.split()[1:])
							LC_viral_predictions.append(record)
							break

		if len(LC_viral_predictions) < 1:
			print("No low confidence viral contigs were reported")

		else:
			print("Low confidence viral predictions identified")

		print("Beginning the identification of prophage predictions")

		# Retrieve records from VirSorter prophage predictions and create unified list

		search_id = re.compile(r"VIRSorter_(\w+)[a-z-]*_(gene_\d+_gene_\d+[0-9-]+cat_\d)")
		for item in VirSorter_prophages:
			if os.stat(item).st_size != 0:
				for prophage in SeqIO.parse(item, "fasta"):
					prophage_description = search_id.search(prophage.description).group(1)
					prophage_suffix = search_id.search(prophage.description).group(2)
					for record in SeqIO.parse(assembly_file, "fasta"):
						if re.sub(r"[.,:; ]", "_", record.description) == prophage_description:
							prophage.id = "%s_%s" % (record.id, prophage_suffix)
							prophage.description = "_".join(record.description.split()[1:])
							prophage_predictions.append(prophage)
							break

		if len(prophage_predictions) < 1:
			print("No putative prophages were reported by VirSorter in categories 4 and 5")

		else:
			print("Prophage predictions identified")


	else:
	
		print("No putative viral sequences were reported by VirSorter")

		if len(VF_high_ids) > 0:

			print("Beginning the identification of low confidence viral predictions by VirFinder...")

			for item in VF_high_ids:
				for record in SeqIO.parse(assembly_file, "fasta"):
					if record.description == item:
						record.id = record.id + "_10_H"
						record.description = "_".join(record.description.split()[1:])
						LC_viral_predictions.append(record)
						break

			print("Low confidence viral predictions by VirFinder identified")

	# Return viral predictions sorted in different confidence categories and prophage predictions

	return (HC_viral_predictions, LC_viral_predictions, prophage_predictions)
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Write fasta files with predicted viral contigs sorted in categories and putative prophages")
	parser.add_argument("-f", "--fasta", dest = "input_file", help = "Assembly fasta file", required = True)
	parser.add_argument("-o", "--outdir", dest = "outdir", help = "output directory", default = ".")
	parser.add_argument("-t", "--thread", dest = "thread", type = int, help = "Number of threads to use for the computation (default: 4)", default = "4")
	parser.add_argument("--virome", dest = "virome", action = "store_true", help = "Indicate whether you are processing a metagenomic or viromic dataset")
	parser.add_argument("--prok", dest = "mode", action = "store_true", help = "Indicate whether you would like the pipeline to focus only on prokaryotic viruses")
	if len(sys.argv) == 1:
		parser.print_help()
	else:
		args = parser.parse_args()
		viral_predictions = virus_pred(args.input_file, args.outdir, args.thread, args.virome, args.mode)
		if sum([len(x) for x in viral_predictions]) > 0:
			if len(viral_predictions[0]) > 0:
				SeqIO.write(viral_predictions[0], os.path.join(args.outdir, "High_confidence_putative_viral_contigs.fna"), "fasta")
			if len(viral_predictions[1]) > 0:
				SeqIO.write(viral_predictions[1], os.path.join(args.outdir, "Low_confidence_putative_viral_contigs.fna"), "fasta")
			if len(viral_predictions[2]) > 0:
				SeqIO.write(viral_predictions[2], os.path.join(args.outdir, "Putative_prophages.fna"), "fasta")
		else:
			print("Overall, no putative viral contigs or prophages were detected in the analysed metagenomic assembly")
