#! /bin/bash

mkdir -p raw_data
cat $1 | xargs -i fastq-dl -o raw_data -a {}

# clean up
find raw_data -type f \! -name "*_*" -exec rm {} +
rename 's/(.*)_/$1_R/' raw_data/*.fastq.gz
