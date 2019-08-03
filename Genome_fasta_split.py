#!/usr/bin/env python3
# coding: utf-8

# In[1]:


import random
from Bio import SeqIO
import os


# In[13]:


min_len = 500


# In[18]:


with open("Data/bacterial_genome_chunks.fna", "w") as output_file:
    for record in SeqIO.parse("Data/Bacterial_genomes.fna", "fasta"):
        contig_n = 0
        present_len = 0
        max_len = 50000
        total_len = len(record)
        while present_len < total_len:
                contig_n += 1
                if max_len > total_len - present_len:
                    max_len = total_len - present_len
                try:
                    dist = random.sample(range(min_len, max_len), 1)
                    output_file.write(">%s_%s\n%s\n" % (record.id, contig_n, record.seq[present_len:present_len + dist[0]]))
                    present_len += dist[0]
                except:
                    break

