#!/bin/bash
set -e

# this run the extraction and alignment of the Clermont PCR amplicons.
# its pretty slow; the bottleneck is the extractRegion script; i think its O(n),
# its a tradeoff of ease of the ortholog searching (by combing all the sequences
# into a single blast database) and the extraction process.
# TODO: make it better

DATE=`date +%Y-%m-%d`
# make the output directory
OUTDIR="./${DATE}_clermont_region_extraction/"
mkdir $OUTDIR

# for each of the amplicons in the Clermnt scheme
for GENEF in refs/*.fasta
do
  # get the name of the gene from the file name
  BN=$(basename $GENEF)
  REG=${BN%.fasta}
  echo $REG
  # find the coords of the ortholog
  ~/GitHub/open_utils/orthoML/simpleOrtho.py $GENEF ./combined/ -n -v 1 -o $OUTDIR/${REG}_out -t 4
  # get the sequence of the othologs
  ~/GitHub/open_utils/extractRegion/extractRegion.py -l ./$OUTDIR/${REG}_out/simpleOrtho_regions.txt ./combined.fasta > $OUTDIR/${REG}_regions.fasta
  # align all those
  mafft --thread 4 ./$OUTDIR/${REG}_regions.fasta > ./$OUTDIR/${REG}_regions_aligned.fasta
  # get rid of big temp files
  rm $OUTDIR/${REG}_out -rf
done
