#!/usr/bin/env python3 

import sys
import numpy as np
import os
from Bio import SeqIO

dataset_id = sys.argv[1].rstrip(".fna").split("_")[0]

length_list = []

record_generator = SeqIO.parse(sys.argv[1], "fasta")

for record in record_generator:
	length_list.append(len(record))

length_array = np.asarray(length_list)

print("Dataset: %s\nNumber of contigs: %s\nMin length: %s\nq25: %s\nq50: %s\nq75: %s\nMax length: %s\nMean: %s\n" % (dataset_id, str(len(length_list)), str(np.min(length_array)), str(np.quantile(length_array, 0.25)), str(np.quantile(length_array, 0.5)), str(np.quantile(length_array, 0.75)), str(np.max(length_array)), str(np.mean(length_array))))
