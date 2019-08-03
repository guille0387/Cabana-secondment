#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

This script detects and characterizes metagenome viral contigs (mVCs) in a metagenome assembly

OPTIONS:
   -h      Show help message
   -a      Assembly file (.fasta) [REQUIRED]
   -d      HMM file containing virus families (.hmms) [REQUIRED]
   -b	   BLAST db to classify contigs (.fasta) [REQUIRED]
   -r	   Forward or single fastq file for abundance estimation (.fastq or .fastq.gz) [REQUIRED]
   -n	   Reverse fastq file for abundance estimation (.fastq or .fastq.gz) [OPTIONAL]
   -l	   Read length for coverage estimation [REQUIRED]
   -t	   Tab file with virus to family match [REQUIRED]
   -j      Add -j if running script on a job submission system [OPTIONAL]
EOF
}

ASSEMB_FILE=
HMM_FILE=
BLAST_DB=
JOB_SYS=
READS=
READS_2=
READ_LENGTH=
VIR_FAMILY=


while getopts “ha:d:b:r:n:l:t:j” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         a)
             ASSEMB_FILE=$OPTARG
             ;;
         d)
             HMM_FILE=$OPTARG
             ;;
	 b)
	     BLAST_DB=$OPTARG
	     ;;
	 r)
	     READS=$OPTARG
	     ;;
	 n)
	     READS_2=$OPTARG
	     ;;
	 l)
	     READ_LENGTH=$OPTARG
	     ;;
	 t)
	     VIR_FAMILY=$OPTARG
	     ;;
         j)
	     JOB_SYS=1
	     ;;
         ?)
             usage
             exit
             ;;
     esac
done

timestamp() {
  date +"%H:%M:%S"
}


echo "$(timestamp) [mVCs script] Parsing command-line"
echo "mVCs_from_assembly.sh -a ${ASSEMB_FILE} -d ${HMM_FILE} -b ${BLAST_DB} -r ${READS} -n ${READS_2} -l ${READ_LENGTH} -j ${JOB_SYS}"
if [[ -z ${ASSEMB_FILE} ]] || [[ -z ${HMM_FILE} ]] || [[ -z ${BLAST_DB} ]] || [[ -z ${READS} ]] || [[ -z ${READ_LENGTH} ]]
then
     echo "ERROR : Please supply the required arguments"
     usage
     exit 1
fi

# set output names
ori_name=$(basename ${ASSEMB_FILE%%.fasta})
vir_name=${ori_name%%.fasta}"_virf"
folder=${vir_name}"_mVCs"
hmm_prsd=${vir_name}"_hmmsearch_parsed.tab"
hmm_abund=${vir_name}"_hmmsearch_parsed_abund.tab"

if [[ ! -e ${vir_name}".fasta" ]]
then
	echo "$(timestamp) [mVCs script] Running VirFinder"
	VirFinder_analysis.R -f ${ASSEMB_FILE}
	select_seqs_by_ID.py -i ${ASSEMB_FILE} -d ${ori_name}"_VirFinder_contigs.txt" -o ${vir_name}".fasta"
fi

if [[ -s ${vir_name}".fasta" ]]
then
	echo "$(timestamp) [mVCs script] Running detection python script"
	detect_viral_contigs.py -a ${vir_name}".fasta" -p ${HMM_FILE}

	if [[ ! -e ${vir_name}".fasta.bwt" ]]
	then
		echo "$(timestamp) [mVCs script] Indexing reference file"
		bwa index ${vir_name}".fasta"
	fi

	if [[ ! -e ${hmm_abund} ]]
	then
		echo "$(timestamp) [mVCs script] Running BWA for abundance estimation"
		if [[ -z ${READS_2} ]]
		then
			echo "bwa mem ${vir_name}".fasta" ${READS} | samtools view -uS - -o ${vir_name}_unsorted.bam"
			bwa mem ${vir_name}".fasta" ${READS} | samtools view -uS - -o ${vir_name}"_unsorted.bam"
		else
			echo "bwa mem ${vir_name}.fasta ${READS} ${READS_2} | samtools view -uS - -o ${vir_name}_unsorted.bam"
			bwa mem ${vir_name}".fasta" ${READS} ${READS_2} | samtools view -uS - -o ${vir_name}"_unsorted.bam"
		fi

		samtools sort ${vir_name}"_unsorted.bam" -o ${vir_name}"_sorted.bam"

		echo "$(timestamp) [mVCs script] Indexing and analysing mapping file"
		samtools index ${vir_name}"_sorted.bam"
		samtools idxstats ${vir_name}"_sorted.bam" > ${vir_name}"_coverage.tab"

		echo "$(timestamp) [mVCs script] Parsing coverage file and outputting results"
		cov_viral_contigs.py -c ${vir_name}"_coverage.tab" -p ${hmm_prsd} -r ${READ_LENGTH}
	fi

	echo "Checking ${folder} and ${hmm_abund}"
	if [[ -d ${folder} ]] && [[ -e ${hmm_abund} ]]
	then
        	echo "$(timestamp) [mVCs script] Running classification python script"
        	for contig in ${folder}"/"*".fasta"
        	do
                	if [[ ! -z ${JOB_SYS} ]]
                	then
				sleep 5
                        	bsub -q production-rh7 -M 50000 -n 8 -o ${vir_name}_classmVCs.log -J ${vir_name}_classmVCs\
                        	"classify_viral_contigs.py -c ${contig} -p ${hmm_abund} -d ${BLAST_DB} -f ${VIR_FAMILY}"
                	else
                        	classify_viral_contigs.py -c ${contig} -p ${hmm_abund} -d ${BLAST_DB} -f ${VIR_FAMILY}
                	fi
        	done
	else
        	echo "$(timestamp) [mVCs script] Oops, something went wrong!"

	fi
else
	echo "$(timestamp) [mVCs script] No contigs passed the VirFinder and length filter, exiting"
fi
