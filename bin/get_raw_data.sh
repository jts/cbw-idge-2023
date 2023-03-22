#! /bin/bash

mkdir -p raw_data
cat $1 | xargs -i fastq-dl -o raw_data -a {}
