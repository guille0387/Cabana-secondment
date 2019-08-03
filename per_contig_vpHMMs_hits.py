#!/usr/bin/env python3

import re
import os


def per_contig_vphmms(hmm_output):
	"""Generate table indicating vphmms hits per contig"""
	with open(hmm_output, "r") as input_file:
		pattern = re.compile(r"[A-Z]{3}\d{6}_\d+")
		contig_list = []
		for line in input_file:
			match = pattern.search(line)
			if match != None and match.group(0) not in contig_list:
				contig_list.append(match.group(0))
		input_file.seek(0)
		output_name = hmm_output.rstrip("_tblout.txt") + "_vphmms_per_contig.tsv"
		with open(output_name, "w") as output_file:
			output_file.write("Contig" + "\t" + "vpHMM_hit" + "\t" + "Total_E-value" + "\t" + "Best_Dom_E-value" + "\t" + "Num_Dom" + "\n")
			for contig in contig_list:
				pattern = contig + "_"
				for i in range(3):
					next(input_file)
				for line in input_file:
					if pattern in line:
						line = line.split()
						output_file.write(contig + "\t" + line[0] + "\t" + line[4] + "\t" + line[7] + "\t" + line[10] + "\n")
				input_file.seek(0)



os.chdir("/hps/nobackup2/production/metagenomics/garp/tests/hmmer_output")

file_list = os.listdir()

search_pattern = re.compile("_tblout")

query_list = [x for x in file_list if search_pattern.search(x) != None]

for item in query_list:
	per_contig_vphmms(item)
