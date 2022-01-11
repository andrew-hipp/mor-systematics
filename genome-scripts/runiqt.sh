#!/bin/bash
## Source: Ryan Fuller, 2021-12-17
#usage: ./runiqtree your.algn.fila.fasta

#input a gene alignment
aln=$1 # alignment (file)
tMax=$2 # maximum number of threads
boots=$3 # bootstraps
prefix=$4 # prefix for files out

# AH modified to use iqtree instead of iqtree2
# AH modified to add arguments for threads, boots, prefix
iqtree -s $aln -B $boots -T AUTO --threads-max $tMax -pre $prefix

for suf in ckp.gz mldist bionj model.gz ; do
    if [ -e $aln.$suf ] ; then
        rm $aln.$suf
    fi
done

if [ -e $aln.log ] ; then
    tar cvfz $aln.iqt.tgz \
        --remove-files \
        --ignore-failed-read \
        $aln.log \
        $aln.treefile \
        $aln.iqtree \
        $aln.contree \
        $aln.splits.nex \
        $aln.uniqueseq.phy ;
fi
