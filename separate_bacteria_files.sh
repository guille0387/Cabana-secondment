#!/usr/bin/env bash

counter=1
bacteria=1

for line in $(cat Bacteria_ftp_refseq_urls.txt); do
if [ ${bacteria} -le 2000 ]
then
	wget -P Bacteria_${counter} ${line}
	echo ${bacteria}
	echo ${counter}
	bacteria=$((${bacteria}+1))
else
	counter=$((${counter}+1))
	bacteria=1
	wget -P Bacteria_${counter} ${line}
	echo ${bacteria}
	echo ${counter}
	bacteria=$((${bacteria}+1))
fi;
done 
