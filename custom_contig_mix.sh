#!/usr/bin/env bash

total_seqs=20000

for item in $(seq 2 2 8); do
	viral_seqs=$(( ${total_seqs} * ${item} / 10 ))
	bacterial_seqs=$(( ${total_seqs} - ${viral_seqs} ))
	grep '>' viral_genome_chunks.fna | shuf | head -n ${viral_seqs} >> $(( ${item} * 10 ))_percent_viral.txt
	grep '>' bacterial_genome_chunks.fna | shuf | head -n ${bacterial_seqs} >> $(( ${item} * 10 ))_percent_viral.txt
done
