#1 /bin/bash

cat etc/accessions.txt | xargs -i bin/downsample_accession.sh

