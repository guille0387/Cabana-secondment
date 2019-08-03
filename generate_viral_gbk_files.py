#!/usr/bin/env python3.7

import os
import sys
import re
import argparse
import operator
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import IUPAC	
from Bio.SeqFeature import SeqFeature, FeatureLocation

sys.path.insert(0, "/homes/garp/scripts")

from write_viral_predictions import virus_pred
from Generate_vphmm_hmmer_matrix import hmmer_domtbl
from Ratio_Evalue_table import ratio_evalue

parser = argparse.ArgumentParser(description = "Generate genbank files of viral contigs and putative prophage sequences")
parser.add_argument("-f", "--fasta", dest = "input_file", help = "Metagenomic assembly fasta file", required = True)
if len(sys.argv) == 1:
	parser.print_help()
else:
	args = parser.parse_args()
	input_file = args.input_file
	dataset_id = input_file.split("_")[0]
	contig_length_filter = "/hps/nobackup2/production/metagenomics/aalmeida/scripts/EMBL-EBI/filter_contigs_len.py"
	subprocess.call("%s -f %s -l 1" % (contig_length_filter, input_file), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)
	new_name_search = re.search(r"(%s\w+)\.fna(\w+\.fasta)" % dataset_id, ",".join(os.listdir()))
	new_name = new_name_search.group(1) + new_name_search.group(2)
	os.rename(new_name_search.group(0), new_name)
	input_file = new_name
	predicted_viruses = virus_pred(input_file)
	SeqIO.write(predicted_viruses, "%s_viral_sequences.fna" % dataset_id, "fasta")
	input_file = "%s_viral_sequences.fna" % dataset_id
	hmmer_result = hmmer_domtbl(input_file)
	informative_df = ratio_evalue(hmmer_result)
	os.mkdir("%s_annotated_viral_sequences" % dataset_id)
	for contig in SeqIO.parse(input_file, "fasta", IUPAC.IUPACUnambiguousDNA()):
		for protein in SeqIO.parse("%s_viral_CDS.faa" % dataset_id, "fasta", IUPAC.IUPACProtein()):
			if contig.id in protein.id:
				protein_fields = protein.description.split(" # ")
				CDS_annotation = SeqFeature(FeatureLocation(int(protein_fields[1]), int(protein_fields[2])), type = "CDS", strand = int(protein_fields[3]))
				CDS_annotation.qualifiers["locus_tag"]  = [protein.id]
				CDS_annotation.qualifiers["transl_table"] = [11]
				CDS_annotation.qualifiers["translation"] = [str(protein.seq)]
				if protein.id in list(informative_df["query"].values):
					if list(informative_df["query"].values).count(protein.id) > 1:
						best_hit = max(informative_df[informative_df["query"] == protein.id]["Abs_Evalue_exp"].items(), key = operator.itemgetter(1))
						if best_hit[1] >= 10:
							CDS_annotation.qualifiers["result"] = ["high confidence"]
						else:
							CDS_annotation.qualifiers["result"] = ["low confidence"]
						CDS_annotation.qualifiers["taxon"] = [informative_df.loc[best_hit[0], "Taxon"], informative_df.loc[best_hit[0], "Abs_Evalue_exp"]]
					else:
						query_index = informative_df["query"][informative_df["query"] == protein.id].index[0]
						if informative_df["Abs_Evalue_exp"][query_index] >= 10:
							CDS_annotation.qualifiers["result"] = ["high confidence"]
						else:
							CDS_annotation.qualifiers["result"] = ["low confidence"]
						CDS_annotation.qualifiers["taxon"] = [informative_df.loc[query_index, "Taxon"], informative_df.loc[query_index, "Abs_Evalue_exp"]]
				else:
					CDS_annotation.qualifiers["result"] = ["no hit"]
					CDS_annotation.qualifiers["taxon"] = ["NA"]
				contig.features.append(CDS_annotation)
		SeqIO.write(contig, "%s_annotated_viral_sequences/%s.gbk" % (dataset_id, contig.id), "genbank")
