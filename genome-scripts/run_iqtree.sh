#!/bin/bash
## Source: Ryan Fuller, 2021-12-17
#usage: ./runiqtree your.algn.fila.fasta

#input a gene alignment
aln=$1

iqtree2 -s $aln -B 1000 -T AUTO --threads-max 4

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
