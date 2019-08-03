#!/usr/bin/env python3

import os
import re
import glob
import sys
import operator
import pandas as pd
import matplotlib.pyplot as plt


hmmer_result_files = glob.glob("ViPhOG_database_analysis/hmm*domtbl/*")

viphog_unique = []

for elem in hmmer_result_files:
	targets = []
	unique_targets = []
	with open(elem) as input_file:
		for line in input_file:
			if not re.search(r"^#", line) and float(line.split()[6]) <= 1e-3:
				targets.append(line.split()[0].split("|")[1])
    
	if len(targets) > 0:
		for item in targets:
			if item not in unique_targets:
				unique_targets.append(item)
        
		viphog_unique.append((os.path.basename(elem).split("_")[0], len(unique_targets)))

viphog_unique_sorted = sorted(viphog_unique, key = operator.itemgetter(1))

def chunks(l, n):
	chunk_size = int(len(l)/n)
	if len(l)%n == 0:
		for i in range(0, len(l), chunk_size):
			yield l[i:i + chunk_size]
	else:
		for i in range(0, len(l), chunk_size):
			if len(l[i+chunk_size:]) >= chunk_size:
				yield l[i:i + chunk_size]
			else:
				yield l[i:]
				break

viphog_gen = chunks(viphog_unique_sorted, 10)

for i,j in enumerate(viphog_gen):
	with open("ViPhOG_database_analysis/unique_hits_chunks/chunk_%s" % (i + 1), "w") as output_file:
		output_file.write("profile\tunique_hits\n")
		for x,y in j:
			output_file.write("%s\t%s\n" % (x, y))
