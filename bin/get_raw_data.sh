#! /bin/bash

mkdir -p raw_data
cat etc/accessions.txt | xargs -i fastq-dl -o raw_data -a {}
