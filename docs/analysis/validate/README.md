# Getting the validation data
The Strains listed in Clermont 2015  are mostly from bioproject PRJNA321606 (all the ECOR strains), and a handful of others.  The other strains were identified manually by seaching NCBI; 6 strains from the paper were not availible on Genbank. The strains from sims2011 were were identified similarly, and added to the `validation_metadata.csv` sheet.  The accession from this used with the Entrez Batch tool to download; download the Genbank ones rather than the refeq, and refseq is missing a number of the accessions.


# Running EzClermont
The strains were unzipped and were then added to a `./genome_assemblies_genome_fasta/ncbi-genomes-2020-04-06/` directory, and EzClermont was run:

```
for i in ./genome_assemblies_genome_fasta/ncbi-genomes-2020-04-06/*.gz
do
	nam=$(basename $i)
	echo $nam
	cat $i | gunzip  |  ezclermont  - -e $nam >> 2020-04-06-ezclermont.txt
done
```



# Generating kmer-base phylogentic tree


```
mkdir ./genome_assemblies_genome_fasta/tmp/
gunzip ./genome_assemblies_genome_fasta/ncbi-genomes-2020-04-06/GCA_00*
for i in ./genome_assemblies_genome_fasta/ncbi-genomes-2020-04-06/*.gz ; do bas=$(filename $i) ; cat $i | gunzip > ./genome_assemblies_genome_fasta/tmp/${bas}.fasta ; done

/Applications/Harvest-OSX64-v1.1.2/parsnp -c -d ./genome_assemblies_genome_fasta/tmp/ -r !  -p 4 -o ./alignment/
/Applications/Harvest-OSX64-v1.1.2/harvesttools  -i ./alignment/parsnp.ggr -M  ./alignment/parsnp.msa

cat alignment/parsnp.msa | sed -e 's/^>\(.\{5\}\)\(.\{8\}\).*/\2/'  > alignment/parsnp.clean.msa


python -c "import sys,os; print(sys.argv[1]) ; outname=os.path.basename(os.path.splitext(sys.argv[1])[0]); from Bio import AlignIO; alignments = AlignIO.parse(sys.argv[1], 'fasta'); AlignIO.write(alignments, outname + '.phy', 'phylip')"  ./alignment/parsnp.clean.msa

seqret -sequence ./alignment/parsnp.msa -outseq ./alignment/parsnp.phy -osformat2 phylip
Read and write (return) sequences

```

# Timing and running ClermontTyping
```

time for i in ../genome_assemblies_genome_fasta/tmp/GCA_00*.fasta; do   ~/GitHub/ClermonTyping/clermonTyping.sh --fasta $i ; done

cat analysis_2020-04-07_15*/*phylogroups.txt >2020-04-07-CT-results.txt

#  3m7.051s
```

```
time for i in ../genome_assemblies_genome_fasta/tmp/GCA_00*.fasta; do   ezclermont  $i  >> 2020-04-06-ezclermont-timing.txt ; done

# 3m43.620s
```
# Plotting the results, etc
The resulting files were then processed with the `plot_Clermont_comparison.R` script.
