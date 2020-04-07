# Building the training data
In order to determine how much variability there is in the priming sites, we will download a pile of genomes, extract the regions of interest, and make alignments

### Training genomes
The initial analysis used the *E. coli*  genomes from the following bioprojects:

- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA352562/
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA352198/
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA218110/
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA231221/


However, in order to ensure even representation, version 0.6 introduced primer training based on a subset of enterobase where a random isolate from each sequence type (as of mid-2019) was downloaded.

These 1395 strains were moved to  a folder called `combined`.  Because we need to have unique contig names when these are turnned into a blast database, the following command was run to rename the headers with the file prefix:

```
for i in /mnt/shared/projects/Escherichia_coli/201409_environmental_samples/analysis/2019-04-26-enterobase_subset/*fasta; do bb=$(filename $i) ; cat $i | sed "s/>/>${bb}_/g" > combined/${bb}.fasta ; done
```

### Make concatenated genomes file
To use the `extractRegions.py` requires all them to be in a single file, so we concatenate all the genomes:


```
cat combined/*fasta > combined.fasta
```

### Selecting cannonical primer targets

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
ybgD
https://www.ncbi.nlm.nih.gov/assembly/GCF_900499975.1/

Now,for each of the gene fasta files, find the primer site, and trim the sequence to the "amplicon" with +5bp on either end.  Note that this makes the coordinates in the fasta headers technically incorrect.


#  Extracting the relevant regions
We then use `simpleOrtho.py` to find the reciprocal best match of that locus in each of the genomes.  For each gene, we run `simpleOrtho.py`^[https://github.com/nickp60/simpleOrtho], use `extractRegion.py`^[https://github.com/nickp60/open_utils] to get a fasta sequence of that reciprocal match, and once all those have been found, we generate a multiple sequence alignment with mafft using default parameters except for allowing reverse complementation:

```
# see the extract_amplicons.sh script
```

# Incorporating mismatches

Then, inspect the MSAs by identifying the primer site and noting any  noting any ambiguities, and incorporate those into the primer matches in the `main()` method of the `run.py` script.  Ignore differences in the last 5 bases in the 3' end, as those are used to differentiate alleles in some cases.
