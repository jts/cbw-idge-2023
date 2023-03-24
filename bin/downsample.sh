#1 /bin/bash

#ACC=$1
#DATA=$2
#REPO=$2

while getopts "a:r:" option; do
	case "${option}" in
		a) ACC=$OPTARG;;
		r) REPO=$OPTARG;;
	esac
done

while read i; do $(dirname "$0")/downsample_accession.sh -a $i -r $REPO; done < $ACC

#cat etc/accessions.txt | xargs -i bin/downsample_accession.sh

