#!/usr/bin/env python3

import os
from ete3 import NCBITaxa

ncbi = NCBITaxa()

with open("/hps/nobackup2/production/metagenomics/garp/Data/ViPhOG_tax_id.txt") as input_file:
	ViPhOG_lineage_dict = {}
	rank_order = ["order", "family", "subfamily", "genus"]
	for line in input_file:
		line = int(line.rstrip())
		id_name_dict = ncbi.get_taxid_translator([line])
		lineage = ncbi.get_lineage(line)
		rank_dict = ncbi.get_rank(lineage)
		name_dict = ncbi.get_taxid_translator(lineage)
		lineage_list = []
		for rank in rank_order:
			for k,v in rank_dict.items():
				if v == rank:
					lineage_list.append(name_dict[k])
					break
		ViPhOG_lineage_dict[id_name_dict[line]] = lineage_list
	ViPhOG_lineage_dict["Arterivirus"] = ViPhOG_lineage_dict["Dipartevirus"]
	ViPhOG_lineage_dict["Arterivirus"][2] = "Arterivirus"
	del(ViPhOG_lineage_dict["Dipartevirus"])
