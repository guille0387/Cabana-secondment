#!/usr/bin/env python3.7

import os
import subprocess
import argparse
import re
import pandas as pd
import sys

def hmmer_domtbl(viral_sequence_file, output_dir, threads):
	"""This function takes a fasta file with predicted viral contigs and putative prophage sequences,
	   identifies their putative CDS, compare the corresponding proteins against the ViPhOG HMM database
	   and generates a pandas dataframe that stores the obtained results""" 
	hmm_database_path = "/hps/nobackup2/production/metagenomics/garp/Data/hmmFiles/vpHMM_database"
	output_file_name = os.path.join(output_dir, re.split(r"\.[a-z]+$", os.path.basename(viral_sequence_file))[0])
	prodigal_output = "%s_CDS.faa" % output_file_name
	hmmer_output = "%s_hmmer_ViPhOG.tbl" % output_file_name
	subprocess.call("prodigal -a %s -i %s -p meta" % (prodigal_output, viral_sequence_file), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)
	subprocess.call("hmmscan --domtblout %s --noali -E 0.001 --cpu %s %s %s" % (hmmer_output, threads, hmm_database_path, prodigal_output), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)
	hmmer_output_lines = []
	with open(hmmer_output, "r+") as output_file:
		for line in output_file:
			if not re.search(r"^#", line):
				line = re.sub("\s+", "\t", line.rstrip())
				hmmer_output_lines.append(line)
		output_file.seek(0)
		output_file.write("target name\ttarget accession\ttlen\tquery name\tquery accession\tqlen\tfull sequence E-value\tfull sequence score\tfull sequence bias\t#\tof\tc-Evalue\ti-Evalue\tdomain score\tdomain bias\thmm coord from\thmm coord to\tali coord from\tali coord to\tenv coord from\tenv coord to\tacc\tdescription of target\n")
		output_file.write("\n".join(hmmer_output_lines))
		output_file.truncate()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Generate ViPhOG hit table with 0.001 E-value threshold for viral contigs file")
	parser.add_argument("-f", "--fasta", dest = "input_file", help = "Viral and prophage contigs file", required = True)
	parser.add_argument("-o", "-outdir", dest = "output_dir", help = "Output directory for storing hmmer output file", default = ".")
	parser.add_argument("-t", "--thread", dest = "thread", type = int, help = "Number of threads to use for the computation (default: 4)", default = "4")
	if len(sys.argv) == 1:
		parser.print_help()
	else:
		args = parser.parse_args()
		input_file = args.input_file
		output_dir = args.output_dir
		num_thread = args.thread
		hmmer_domtbl(input_file, output_dir, num_thread)
