#!/usr/bin/env python

import os
import sys
import argparse

def covStats(infile, readlength):
	filein = open(infile, "r")
	contig_info = {}
	totalMapped = 0
	totalUnmapped = 0
	totalLen = 0
	totalCov = 0.0
	for line in filein:
		contig = line.split()[0]
		length = int(line.split()[1])
		mapped = int(line.split()[2])
		unmapped = int(line.split()[3])
		if contig != "*":
			cov = float(mapped*readlength)/length
			contig_info[contig] = [length, mapped, unmapped, cov]
			totalCov += cov
		totalMapped += mapped
		totalUnmapped += unmapped
		totalLen += length
	totalReads = totalMapped+totalUnmapped
        scale_factor = float(totalReads/1000000)
	filein.close()
	return (contig_info, totalCov, scale_factor)	


def parseCov(hmm_parsed, covInfo):
	contig_info = covInfo[0]
	totalCov = covInfo[1]
	scale_factor = covInfo[2]
	file_name = os.path.basename(hmm_parsed)
	out_name = file_name.split(".tab")[0]+"_abund.tab"
	out_file = open(out_name, "w")
	hmm = open(hmm_parsed, "r")
	line_n = 0
	for entry in hmm:
		line_n += 1
		if line_n == 1:
			print >> out_file, entry.strip("\n")+"\tMapped_reads\tCoverage\tCov_perc\tRPKM"
		else:
			contig = entry.split()[1]
			propCov = contig_info[contig][3]/totalCov*100
			RPM = float(contig_info[contig][1])/scale_factor
			RPKM = RPM/(float(contig_info[contig][0])/1000)
			print >> out_file, "%s\t%i\t%.2f\t%.4f\t%.2f" %\
			(entry.strip("\n"), contig_info[contig][1], contig_info[contig][3], propCov, RPKM)
	hmm.close()
	out_file.close()

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Estimate relative abundance of each contig')
        parser.add_argument('-c', dest='coverage', help='samtools idxstats output file (.tab)', required=True)
        parser.add_argument('-p', dest='hmm_parsed',\
                                  help='HMM output tab file', required=True)
        parser.add_argument('-r', dest='read_length', help='read length', required=True)
	if len(sys.argv) == 1 :
                parser.print_help()
                sys.exit(1)
        else:
                # preparing arguments
                args = parser.parse_args()
                infile = args.coverage
                rlength = float(args.read_length)
                hmm_file = args.hmm_parsed
                # running functions
		print "Calculating coverage metrics"
                outCov = covStats(infile, rlength)
		print "Annotating hmm file"
                parseCov(hmm_file, outCov)
                print "Finished abundance analysis"
