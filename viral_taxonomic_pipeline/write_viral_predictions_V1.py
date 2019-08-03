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
	"""Creates fasta file with viral contigs and putative prophages predicted with VirFinder_Euk_Mod and Virsorter"""
	#VirFinder run using the VF.modEPV_k8.rda prediction model for prokaryotic and eukaryotic viruses
	#The output of the script is a file recording the prediction results for all contigs and another file
	#that only includes those results for which fdr < 0.1
	if prok_mode:
		print("Running VirFinder only for prokaryotic viruses")
		subprocess.call("VirFinder_analysis_Prok.R -f %s -o %s" % (assembly_file, output_dir), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True) 
	else:
		print("Running VirFinder for prokaryotic and eukaryotic viruses")
		subprocess.call("VirFinder_analysis_Euk.R -f %s -o %s" % (assembly_file, output_dir), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)
	#Path to directory containing VirSorter's databases
	virsorter_data_path = "/hps/nobackup2/production/metagenomics/garp/databases/virsorter-data"
	#Run VirSorter using the RefSeq + Virome database and the decontamination mode for viral enriched samples
	#if virome_dataset = True, otherwise run VirSorter without desontamination mode
	if virome_dataset:
		print("Running VirSorter with virome decontamination mode")
		subprocess.call("wrapper_phage_contigs_sorter_iPlant.pl -f %s --db 2 --wdir %s --ncpu %s --virome --data-dir %s" % (assembly_file, output_dir, num_thread, virsorter_data_path), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)
	else:
		print("Running VirSorter without virome decontamination mode")
		subprocess.call("wrapper_phage_contigs_sorter_iPlant.pl -f %s --db 2 --wdir %s --ncpu %s --data-dir %s" % (assembly_file, output_dir, num_thread, virsorter_data_path), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)
	VF_viral_pred = []
	VS_viral_pred = {}
	VS_prophage_pred = []
	viral_record_list = []
	#Save ids and descriptions of contigs predicted as viral by VirFinder (fdr < 0.1)
	VirFinder_file = glob.glob(os.path.join(output_dir, "*VirFinder_table-all.tab"))[0]
	VirFinder_result_df = pd.read_csv(VirFinder_file, sep = "\t")
	VirFinder_sig_df = VirFinder_result_df[VirFinder_result_df["fdr"] < 0.1]
	if len(VirFinder_sig_df) > 0:
		VF_viral_pred.extend(list(VirFinder_sig_df["name"].values))
	else:
		print("No putative viral contigs were reported by VirFinder using fdr < 0.1 as inclusion threshold")
	#Store paths of files containing VirSorter output sorted in categories 1, 2, 4 and 5
	VirSorter_viral_list = [x for x in glob.glob(os.path.join(output_dir, "Predicted_viral_sequences", "*")) if re.search(r"cat-[12]\.fasta", x)]
	VirSorter_prophage_list = [x for x in glob.glob(os.path.join(output_dir, "Predicted_viral_sequences", "*")) if re.search(r"cat-[45]\.fasta", x)]
	#Save ids and details of contigs predicted as viral by VirSorter in categories 1 and 2
	id_search = re.compile(r"VIRSorter_(\w+)-([a-z-]*cat_\d)")
	for item in VirSorter_viral_list:
		if os.stat(item).st_size != 0:
			with open(item) as input_file:
				for line in input_file:
					if id_search.search(line):
						VS_viral_pred[id_search.search(line).group(1)] = id_search.search(line).group(2)
	if len(VS_viral_pred) < 1:
		print("No putative viral contigs were reported by VirSorter in categories 1 and 2")
	#Save category 4 and 5 prophage predictions obtained with VirSorter
	id_search = re.compile(r"VIRSorter_(\w+)_(gene_\d+_gene_\d+[0-9-]+cat_\d)")
	for item in VirSorter_prophage_list:
		if os.stat(item).st_size != 0:
			for prophage in SeqIO.parse(item, "fasta"):
				prophage_description = id_search.search(prophage.description).group(1)
				prophage_suffix = id_search.search(prophage.description).group(2)
				for record in SeqIO.parse(assembly_file, "fasta"):
					if re.sub(r"[.,:; ]", "_", record.description) == prophage_description:
						prophage.id = "%s_%s" % (record.id, prophage_suffix)
						prophage.description = "_".join(record.description.split()[1:])
						VS_prophage_pred.append(prophage)
						break
	if len(VS_prophage_pred) < 1:
		print("No putative prophages were reported by VirSorter in categories 4 and 5")
	#Retrieve SeqRecord objects of viral contigs predicted by VirFinder and VirSorter
	for record in SeqIO.parse(assembly_file, "fasta"):
		if record.description in VF_viral_pred or re.sub(r"[.,:; ]", "_", record.description) in list(VS_viral_pred.keys()):
			if record.description in VF_viral_pred and re.sub(r"[.,:; ]", "_", record.description) in list(VS_viral_pred.keys()):
				if "circular" in VS_viral_pred[re.sub(r"[.,:; ]", "_", record.description)]:
					record_suffix = "11_%s_circular" % VS_viral_pred[re.sub(r"[.,:; ]", "_", record.description)].split("_")[1]
				else:
					record_suffix = "11_%s" % VS_viral_pred[re.sub(r"[.,:; ]", "_", record.description)].split("_")[1]
			elif record.description in VF_viral_pred:
				record_suffix = "10"
			elif re.sub(r"[.,:; ]", "_", record.description) in list(VS_viral_pred.keys()):
				if "circular" in VS_viral_pred[re.sub(r"[.,:; ]", "_", record.description)]:
					record_suffix = "01_%s_circular" % VS_viral_pred[re.sub(r"[.,:; ]", "_", record.description)].split("_")[1]
				else:
					record_suffix = "01_%s" % VS_viral_pred[re.sub(r"[.,:; ]", "_", record.description)].split("_")[1]
			record.id = "%s_%s" % (record.id, record_suffix)
			record.description = "_".join(record.description.split()[1:])
			viral_record_list.append(record)

	return (viral_record_list, VS_prophage_pred)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Write fasta file with predicted viral contigs and putative prophages")
	parser.add_argument("-f", "--fasta", dest = "input_file", help = "Assembly fasta file", required = True)
	parser.add_argument("-o", "--outdir", dest = "outdir", help = "output directory", default = ".")
	parser.add_argument("-t", "--thread", dest = "thread", type = int, help = "Number of threads to use for the computation (default: 4)", default = "4")
	parser.add_argument("--virome", dest = "virome", action = "store_true", help = "Indicate whether you are processing a metagenomic or viromic dataset")
	parser.add_argument("--prok", dest = "mode", action = "store_true", help = "Indicate whether you would like the pipeline to focus only on prokaryotic viruses")
	if len(sys.argv) == 1:
		parser.print_help()
	else:
		args = parser.parse_args()
		viral_sequences = virus_pred(args.input_file, args.outdir, args.thread, args.virome, args.mode)
		if sum([len(x) for x in viral_sequences]) == 0:
			print("Overall, no putative viral contigs or prophages were reported for the analysed metagenomic dataset")
		else:
			if len(viral_sequences[0]) > 0:
				SeqIO.write(viral_sequences[0], os.path.join(args.outdir, "Putative_viral_contigs.fna"), "fasta")
			if len(viral_sequences[1]) > 0:
				SeqIO.write(viral_sequences[1], os.path.join(args.outdir, "Putative_prophages.fna"), "fasta")
