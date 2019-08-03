#!/usr/bin/env python

import os
import sys
import argparse
import glob
import subprocess
from Bio import SeqIO

def CDS_annot(assemb):
	CDS_count = {}
	out_file = os.path.basename(assemb).split(".fasta")[0]+"prod.out"
	out_pep = os.path.basename(assemb).split(".fasta")[0]+"_pep.fasta"
	out_nuc = os.path.basename(assemb).split(".fasta")[0]+"_nuc.fasta"
	if not os.path.exists(out_file):
		cmd = ["prodigal",
			"-p", "meta",
			"-i", assemb,
			"-o", out_file,
			"-a", out_pep,
			"-d", out_nuc]
		subprocess.check_call(cmd)
	sample = out_pep.split("_pep")[0]
	cds_in = open(out_pep, "r")
	for record in SeqIO.parse(cds_in, "fasta"):
		contig = "_".join(record.id.split("_")[:-1])
		if contig not in CDS_count.keys():
			CDS_count[contig] = 1
		else:
			CDS_count[contig] += 1
	print "Counting seqs in %s" % (sample)
	return CDS_count
	cds_in.close()


def hmmer_srch(assemb, hmm_file, evalue):
	in_pep = os.path.basename(assemb).split(".fasta")[0]+"_pep.fasta"
	sample = in_pep.split("_pep")[0]
	hmm_out = "%s_hmmsearch.out" % (sample)
	if not os.path.exists(hmm_out):
		cmd =   ["hmmsearch",
			"--tblout", "%s_hmmsearch.out" % (sample),
			"-o", "%s_hmmsearch.out" % (sample),
			"-E", evalue,
			hmm_file, in_pep]
		print "Running hmmsearch with %s on %s" %\
		(os.path.basename(hmm_file), sample)
		subprocess.check_call(cmd) 


def parse_hmmout(assemb, count_CDS):
        sample = os.path.basename(assemb).split(".fasta")[0]
        hmm_out = open("%s_hmmsearch.out" % (sample), "r")
	hits_cds = {}
	hits_contig = {}
	hits_thresh = []
	off_counts = 0
        for entry in hmm_out:
                if entry[0] == "#":
                        off_counts += 1
                elif off_counts == 3:
			fields = entry.split()
			cds = fields[0]
			if cds not in hits_cds.keys():
				hits_cds[cds] = 1
			else:
				hits_cds[cds] += 1
                elif off_counts > 3:
                       break 
        hmm_out.close()
	for prot in hits_cds.keys():
		contig = "_".join(prot.split("_")[:-1])
		if contig not in hits_contig.keys():
			hits_contig[contig] = 1
		else:
			hits_contig[contig] += 1
	out_file = open("%s_hmmsearch_parsed.tab" % (sample), "w")
	print >> out_file, ("Metagenome\tContig\tVirus_CDS\t"\
			    "Total_CDS\tHit_percentage\tResult")
	for contig in count_CDS.keys():
		Metagenome = sample.split("_virf")[0]
		TotalCDS = count_CDS[contig]
		if contig in hits_contig.keys():
			ViralM = hits_contig[contig]
			HitPerc = float(ViralM)/float(TotalCDS)*100
			if HitPerc >= 60:
				hits_thresh.append(contig)
				print >> out_file,\
				 "%s\t%s\t%s\t%s\t%.1f\tmVC" % (Metagenome, contig, ViralM, TotalCDS, HitPerc)
			else:
				print >> out_file,\
				 "%s\t%s\t%s\t%s\t%.1f\t" % (Metagenome, contig, ViralM, TotalCDS, HitPerc)
		else:
			print >> out_file, "%s\t%s\t0\t%s\t0\t" % (Metagenome, contig, TotalCDS)
	out_file.close()
	return hits_thresh


def select_mVCs(assemb, hits):
	in_pep = os.path.basename(assemb).split(".fasta")[0]+"_pep.fasta"
        in_file = list(SeqIO.parse(open(in_pep, "r"), "fasta"))
	sample = in_pep.split("_pep")[0]
	hmm_prsd = "%s_hmmsearch_parsed.tab" % (sample)
	direct = "%s_mVCs" % (sample)
	if not os.path.exists(direct):
		os.makedirs(direct)
	print "Selecting mVCs from hmmer results"
	for hit in hits:
		out_file = open(direct+"/"+hit+".fasta", "w")
		newRecs = []
		for record in in_file:
			contig = "_".join(record.id.split("_")[:-1])
			if contig == hit:
				newRecs.append(record)
		SeqIO.write(newRecs, out_file, "fasta")
		out_file.close()


if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Recover contigs with viral sequences')
        parser.add_argument('-a', dest='assembly', help='Assembly FASTA file', required=True)
	parser.add_argument('-p', dest='hmm_file', help='HMM file to use as query', required=True)
        parser.add_argument('-e', dest='evalue', help='Evalue for phmmer (default: 0.01)'\
				, default="0.01")
	if len(sys.argv) == 1 :
                parser.print_help()
                sys.exit(1)
        else:
		# preparing arguments
		args = parser.parse_args()
                assemb = args.assembly
		hmm_file = args.hmm_file
		evalue = args.evalue
		# running functions
                count_CDS = CDS_annot(assemb)
		hmmer_srch(assemb, hmm_file, evalue)
		hits = parse_hmmout(assemb, count_CDS)
		select_mVCs(assemb, hits)
		print "Finished detection analysis"
