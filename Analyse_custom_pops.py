#!/usr/bin/env python3
# coding: utf-8

# **Import required modules**

# In[1]:


import os
import re
import sys
import glob
import itertools
import argparse
import numpy as np
import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Analyse taxonomic lineage asignment results for in silico communities")
parser.add_argument("-c", dest="com_path", help="Relative or absolute path to directory containing viral pipeline results for in silico community", required=True)
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

else:
	args = parser.parse_args()
	community_dir = args.com_path


# In[2]:


	ncbi = NCBITaxa()



# **Define function for retrieving seq IDs from fasta files**

# In[6]:


	def retrieve_seq_ids(seq_file):
		seq_ids = []
		with open(seq_file) as input_file:
			search_pattern = re.compile(r"^>(.+)\n$")
			for line in input_file:
				if search_pattern.search(line):
					seq_ids.append(search_pattern.search(line).group(1))
		return seq_ids


# **Retrieve IDs of viral genomic fragments**

# In[7]:


	all_viruses_file = "Data/viral_genome_chunks.fna"


# In[8]:


	all_viruses_ids = retrieve_seq_ids(all_viruses_file)


# **Retrieve viral IDs of custom community**

# In[9]:


	community_file = os.path.join("Data", "%s.fna" % os.path.basename(community_dir))
	# Actuvate the following line to keep IDs of all viral contigs in the community
	#community_ids = retrieve_seq_ids(community_file)
	#community_viruses = [x for x in community_ids if x in all_viruses_ids]
	# Activate the following line to keep IDs of viral contigs >= 1 kb
	community_viruses = []
	for record in SeqIO.parse(community_file, "fasta"):
		if record.id in all_viruses_ids and len(record) >= 1000:
			community_viruses.append(record.id)


# **Generate dictionary of viral accession numbers and corresponding TaxIDs**

# In[10]:


	accession_taxid_dict = {}
	with open("Data/viral_genomes_accession_taxid_pairs.tsv") as input_file:
		for line in input_file:
			line = line.rstrip()
			genome_accession, taxid = line.split("\t")[0], int(line.split("\t")[1])
			accession_taxid_dict[genome_accession] = taxid


# **Obtain generator of tuples that contain different combinations of tax assign parameter values**

# In[11]:


	param_comb_gen = itertools.product(range(10, 100, 10), repeat = 2)


# **Analyse taxonomic annotation results**

# In[12]:


	def check_numeric(tax_series):
		tax_mod = tax_series.apply(lambda x: re.sub("\.", "", x) if pd.notnull(x) else x)
		for item in tax_mod:
			if pd.isnull(item) or item.isnumeric():
				continue
			else:
				return False
		return True


# In[13]:


	final_df_rows = []
	for i,j in param_comb_gen:
		tax_assign_files = glob.glob(os.path.join(community_dir, "viral_pred_decon_parse", "tax_assign_benchmark", "*%s_%s.tsv" % (i, j)))
		tax_assign_df = pd.concat([pd.read_csv(x, sep = "\t") for x in tax_assign_files])
		tax_assign_viral_df = tax_assign_df[tax_assign_df["contig_ID"].apply(lambda x: re.search(r"^(\w+\.\d+_\d+)_[a-z0-9A-Z_-]+$", x).group(1) in community_viruses)]
		tax_assign_viral_df = tax_assign_viral_df.reset_index(drop = True)
		true_positives = 0
		false_positives = 0
		false_negatives = 0
		for contig,results in tax_assign_viral_df.iterrows():
			if results.iloc[1:].isnull().all() or check_numeric(results.iloc[1:]):
				false_negatives += 1
			else:
				genome_accession = re.search(r"^(\w+\.\d+)_[a-z0-9A-Z_-]+$", results["contig_ID"]).group(1)
				viral_lineage = ncbi.get_lineage(accession_taxid_dict[genome_accession])
				viral_lineage_names = ncbi.get_taxid_translator(viral_lineage)
				for elem in results.iloc[1:]:
					if elem in viral_lineage_names.values():
						true_positives += 1
						break
				else:
					false_positives += 1
		final_df_rows.append([i, j, true_positives, false_positives, false_negatives])


# **Create dataframe with analysis results**

# In[15]:


	final_df = pd.DataFrame(final_df_rows, columns=["Param1", "Param2", "TP", "FP", "FN"])
	
	#Activate following line if you evaluated all viral contigs
	#final_df.to_csv(os.path.join(community_dir, "viral_pred_decon_parse", "tax_assign_benchmark", "tax_lineage_assignment_results.tsv"), sep = "\t", index = False)

	#Activate following line if you evaluated viral contigs >= 1 kb
	final_df.to_csv(os.path.join(community_dir, "viral_pred_decon_parse", "tax_assign_benchmark", "tax_lineage_assignment_results_1kb.tsv"), sep = "\t", index = False)
