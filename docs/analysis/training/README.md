# Building the training data
In order to determine how much variability there is in the priming sites, we will download a pile of genomes, extract the regions of interest, and make alignments 
Get the genomes from
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA352562/

Also grab the ones from 
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA352198/

and also from PRJNA218110
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA218110/



Move them all into a folder called `combined`, because naming is important. Thats about 500 strains or so.  Now, to use the `extractRegions.py` script that I wrote, its easier if they are all in a single file.  Yes, its clunky, but it works.



```
cat ./ncbi-genomes-2018-04-24/* > ./combined.fasta

```

# Selecting cannonical primer targets

We then get nucleotide fasta files for an instance of each of the loci being targetted.

ArpA and aceK
https://www.ncbi.nlm.nih.gov/gene/944933
ChuA
https://www.ncbi.nlm.nih.gov/gene/7150945
yjaA
https://www.ncbi.nlm.nih.gov/gene/948515
TspE4.C3
https://www.ncbi.nlm.nih.gov/nuccore/7330942
ArpAgpE
https://www.ncbi.nlm.nih.gov/gene/944933
TrpBA
https://www.ncbi.nlm.nih.gov/gene/945768

Now,for each of the gene fasta files, find the primer site, and trim the sequence to the "amplicon" with +5bp on either end
>Note that this makes the coordinates in the fasta header incorrect.  but life is too short to correct those now.


#  Extracting the relevant regions
We then use `simpleOrtho.py` to find the reciprocal best match of that locus in each of the genomes.  For each gene, we run `simpleOrtho.py`, use `extractRegion.py` to get a fasta sequence of that reciprocal match, and once all those have been found, we generate a multiple sequence alignment with mafft using default parameters:
```
for GENEF in refs/*.fasta
do
  BN=$(basename $GENEF)
  GENE=${GENE}.fastay  # double check this line; I forgot to single quote the log
  REG=yjaA
  echo yjaA
  ~/GitHub/open_utils/orthoML/simpleOrtho.py $GENEF ./combined/ -n -v 1 -o ${REG}_out -t 4
  ~/GitHub/open_utils/extractRegion/extractRegion.py -l ./${REG}_out/simpleOrtho_regions.txt ./combined.fasta > ${REG}_regions.fasta
  mafft --thread 4 ./${REG}_regions.fasta > ${REG}_regions_aligned.fasta
  rm ${REG}_out -rf
done

```

# Incorporating mismatches
Then, Inspect the MSAs.  Ignoring the first 5bp at the 3' and 5' ends, note any ambiguities, and incorporate those into the primer matches in the `main()` method of the `run.py` script.