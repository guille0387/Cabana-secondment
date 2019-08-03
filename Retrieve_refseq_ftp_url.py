#!/usr/bin/env python3
# coding: utf-8

# In[1]:


import os
import re
import sys
import pandas as pd
from Bio import SeqIO, Entrez


# In[2]:


Entrez.email = "guille.rangel87@gmail.com"


# In[4]:


#os.chdir("/hps/nobackup2/production/metagenomics/garp/")


# In[64]:


input_file = open("Data/Bacteria_NCBI_accessions.txt")


# In[67]:


output_file = open("Data/Bacteria_ftp_refseq_urls.txt", "w")
for line in input_file:
    acc_num = line.rstrip()
    search_handle = Entrez.esearch(db = "assembly", term = "%s[Accn]" % acc_num, idtype = "acc")
    record = Entrez.read(search_handle)
    search_handle.close()
    print("%s\t%s" % (acc_num, record["Count"]))
    if len(record["IdList"]) > 0:
        fetch_handle = Entrez.efetch(db = "assembly", id = record["IdList"][0], rettype = "docsum")
        result = Entrez.read(fetch_handle)
        fetch_handle.close()
        ftp_url = result["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"] + "/%s_genomic.fna.gz" % re.split("/", result["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"])[-1]
        output_file.write("%s\n" % ftp_url)
output_file.close()
input_file.close()


# In[ ]:




