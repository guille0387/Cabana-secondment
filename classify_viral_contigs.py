#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import operator

def run_BLAST(in_fasta, fasta_db, evalue):
	contig_name = os.path.basename(in_fasta).split(".fasta")[0]
        out_file = in_fasta.split(".fasta")[0]+"_blastp.out"
        hits_per_species = {}
        if not os.path.exists(out_file):
                print "Performing blastp of %s against %s" % (contig_name, os.path.basename(fasta_db))
                cmd =   ['blastp',
                        '-query', in_fasta,
                        '-subject', fasta_db,
                        '-out', out_file,
                        '-outfmt', '6 qseqid sallseqid qcovs pident evalue',
                        '-evalue', evalue,
                        '-max_target_seqs', '1',
                        '-max_hsps', '1']
                subprocess.check_call(cmd)

def store_ref(ref_table):
	vir_table = open(ref_table, "r")
	virus_family = {}
        virus_host = {}
	line_n = 0
	for line in vir_table:
		line_n += 1
		if line_n > 1:
        		virus = line.split("\t")[0]
        		family = line.split("\t")[1]
        		host = line.split("\t")[2].strip("\n")
        		virus_family[virus] = family
        		virus_host[virus] = host
	return (virus_family,virus_host)


def classify_mVCs(in_fasta, hmm_prsd, vir_groups):
	contig_name = os.path.basename(in_fasta).split(".fasta")[0]
	out_file = in_fasta.split(".fasta")[0]+"_blastp.out"
	hits_per_species = {}
	virus_family = vir_groups[0]
	virus_host = vir_groups[1]
	best_family = {}
	best_host = {}
	if os.path.getsize(out_file) > 0:
		print "Analysing and processing blast results"
		for hit in open(out_file, "r"):
			fields = hit.split("\t")
			species = fields[1].split("[")[-1].split("]")[0]
			if species not in hits_per_species:
				hits_per_species[species] = 1
			else:
				hits_per_species[species] += 1
		# convert names to families and select the most frequent
		for name,count in hits_per_species.items():
			if name in virus_family.keys():
				group = virus_family[name]
			else:
				group = "NA"
			if group not in best_family.keys():
				best_family[group] = count
			else:
				best_family[group] += count
			highest = max(best_family.values())
			selected = [k for k,v in best_family.items() if v == highest]
			if len(selected) == 1:
				if selected[0] == "NA":
					end_family = "NA"
				else:
					end_family = "%s (%i CDS)" \
					% (selected[0], highest)
			else:
				end_str = "NA"
		# convert names to hosts and select the most frequent
		for name,count in hits_per_species.items():
                        if name in virus_host.keys():
                                group = virus_host[name]
                        else:
                                group = "NA"
			if "," in group:
				group = group.split(",")
			else:
				group = [group]
			for ele in group:
                        	if ele not in best_host.keys():
                                	best_host[ele] = count
                        	else:
                                	best_host[ele] += count
                        highest = max(best_host.values())
                        selected = [k for k,v in best_host.items() if v == highest]
                        if len(selected) == 1:
                                if selected[0] == "NA":
                                        end_host = "NA"
                                else:
                                        end_host = "%s" \
                                        % (selected[0])
                        else:
                                end_host = "NA"		
		# print final results to output file based on previous selections
		hmm_in = open(hmm_prsd, "r")
		hmm_out_name = hmm_prsd.split(".tab")[0]+"_taxa.tab"
		if os.path.exists(hmm_out_name):
			hmm_out = open(hmm_out_name, "a")
		else:
			hmm_out = open(hmm_out_name, "w")
			print >> hmm_out, \
			("Metagenome\tContig\t"\
			"Virus_CDS\tTotal_CDS\t"\
			"Hit_percentage\tResult\t"\
			"Mapped_reads\tCoverage\t"\
			"Cov_perc\tRPKM\tBLASTP_hits\t"\
			"Family\tHost") 
		sorted_hits = sorted(hits_per_species.items(), key=operator.itemgetter(1), reverse=True)
		for line in hmm_in:
			fields = line.split()
			if fields[1] == contig_name:
				to_append = '; '.join('{}={}'.format(*el) for el in sorted_hits)
				print >> hmm_out, "%s\t%s\t%s\t%s"\
				 % ("\t".join(fields), to_append, end_family, end_host)
		hmm_in.close()
		hmm_out.close()
	else:
		print "No matches in NCBI viral database"
                hmm_in = open(hmm_prsd, "r")
                hmm_out_name = hmm_prsd.split(".tab")[0]+"_taxa.tab"
                if os.path.exists(hmm_out_name):
                        hmm_out = open(hmm_out_name, "a")
                else:
                        hmm_out = open(hmm_out_name, "w")
                        print >> hmm_out, \
			("Metagenome\tContig\t"\
			"Virus_CDS\tTotal_CDS\t"\
			"Hit_percentage\tResult\t"\
			"Mapped_reads\tCoverage\t"\
			"Cov_perc\tRPKM\tBLASTP_hits\t"\
			"Family\tHost")
                for line in hmm_in:
                        fields = line.split()
                        if fields[1] == contig_name:
                                print >> hmm_out, "%s\tNA\tNA\tNA" % ("\t".join(fields))
                hmm_in.close()
                hmm_out.close()

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Annotate viral contigs with known species')
        parser.add_argument('-c', dest='contig', help='Contig FASTA file', required=True)
	parser.add_argument('-p', dest='hmm_prsd', help='Parsed HMM output file', required=True)
	parser.add_argument('-d', dest='db', help='Database with viral CDS (NCBI)', required=True)
	parser.add_argument('-f', dest='family', help='Tab file with virus family info', required=True)
	parser.add_argument('-e', dest='evalue', help='Evalue for blastp (default: 0.001)'\
				, default="0.001")
        if len(sys.argv) == 1 :
                parser.print_help()
                sys.exit(1)
        else:
		# preparing arguments
		args = parser.parse_args()
                contig = args.contig
		hmm_in = args.hmm_prsd
		database = args.db
		vir_group = args.family
		evalue = args.evalue
		# running functions
		run_BLAST(contig, database, evalue)
		ref_dicts = store_ref(vir_group)
		classify_mVCs(contig, hmm_in, ref_dicts)
		print "Finished classifying contigs and added results to hmm_parsed_taxa file"
