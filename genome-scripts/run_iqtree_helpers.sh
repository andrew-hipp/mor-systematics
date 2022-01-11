#!/bin/bash
## source: Ryan Fuller, 2021-12-17
parallel --jobs 3 "./runiqt.sh {}.algn.fasta 1" :::: alignments.txt

#!/bin/bash
# this script will extract the TreeFile from each iqtree archive resulting from the GNU parallel run of 'runiqtree.sh'
# you can amend the second input of the 'tar' command here as you see fit to extract other files or remove it to unarchive all zipped files

while read sample;
do
  tar -xf ${sample}.algn.fasta.iqt.tgz ${sample}.algn.fasta.treefile;
done < genes.txt
