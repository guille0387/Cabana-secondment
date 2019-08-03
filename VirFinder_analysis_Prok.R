#!/usr/bin/env Rscript

# load libraries
library(VirFinder)
library(optparse)
library(purrr)

# prepare arguments
option_list <- list(
  make_option(c("-f", "--fasta"), type="character", default=NULL, 
              help="fasta file", metavar="fasta"),
  make_option(c("-o", "--outdir"), type="character", default=".",
	      help="output dir", metavar="outdir"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$fasta)){
  print_help(opt_parser)
  stop("Provide input FASTA file with -f")
}

# prepare input files
path <- normalizePath(opt$fasta)
#thres <- 0.1

# run VirFinder
cat("Running VirFinder on:",opt$fasta,"\n")
predResult <- VF.pred(path)
qvalue_calc <- possibly(VF.qvalue, NA)
predResult$qvalue <- qvalue_calc(predResult$pvalue)

# correct for multiple testing
predResult$fdr <- p.adjust(predResult$pvalue, method="BH")

# subselect significant contigs
#sign <- predResult[which(predResult$fdr<thres),]
# write to output file
cat("Saving output file","\n")
if (opt$outdir == ".") {
        #sig_table <- paste(sub("\\.f.*", "", opt$fasta), "_VirFinder_table-sig.tab", sep="")
        all_table <- paste(sub("\\.f.*", "", opt$fasta), "_VirFinder_table-all.tab", sep="")
} else {
        #sig_table <- paste(opt$outdir, paste(sub("\\.f.*", "", basename(opt$fasta)), "_VirFinder_table-sig.tab", sep=""), sep="/")
        all_table <- paste(opt$outdir, paste(sub("\\.f.*", "", basename(opt$fasta)), "_VirFinder_table-all.tab", sep=""), sep="/")
}

#write.table(sign, file=sig_table, sep="\t", quote=FALSE, row.names=FALSE)
write.table(predResult, file=all_table, sep="\t", quote=FALSE, row.names=FALSE)
