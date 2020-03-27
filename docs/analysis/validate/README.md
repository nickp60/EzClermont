# Getting the validation data
The Strains listed in Clermont 2015  are mostly from bioproject PRJNA321606 (all the ECOR strains), and a handful of others.  We get the bioproject strains directly ffrom the link on the bioproject page.  The other strains are found in the `other_accessions.txt` doc. 6 strains from the paper were not availible on Genbank. That was used with the Entrez Batch tool to download.

I then manually copied and pasted table 2 from "Guide to the various phylogenetic classification schemes for Escherichia coliand the correspondence among schemes" Olivier Clermont, David Gordon and Erick Denamur


Remember to pad the 0s in the metadata so things match up.


# Running EzClermont
The strains were unzipped and were then added to a `combined` directory, and EzClermont was run:

```
ls ./docs/analysis/validate/combined/*.fa | parallel "ezclermont {} -p 1>> 2018-05-08_p.txt  2>> 2018-05-08_p.log"

# or, without parallel
for i in ./docs/analysis/data/combined/*.fna ;
do
  echo "processing $i"
  clermontpcr $i 1>> 2018-04-30.txt 2>> 2018-04-30.log
done

```

# Generating kmer-base phylogentic tree


```
	# generate the file list
	~/kSNP3.021_Linux_package/kSNP3/MakeKSNP3infile ./combined/ ksnp_manifest
	# make a combined fasta for Kchooser to determinine optimal k
	MakeFasta ./ksnp_manifest ksnp_combined.fasta A
	# run Kchooser
	Kchooser ./ksnp_combined.fasta
	#  found out that k=19 was ideal
	#do a bunch of nonsense to fix perfectly sensible file names
	mkdir renamed_ksnp
	while read path name
	do
	  echo $path
	  echo $namedone
      cp $path ./renamed_ksnp/${name}.fa
	done < ksnp_manifest
	# remake manifest
	MakeKSNP3infile ./renamed_ksnp/ fixed_ksnp_manifest A
	# run kSNP3 whole hog: neighbor joeining, etc.
	~/kSNP3.021_Linux_package/kSNP3/kSNP3 -in fixed_ksnp_manifest -outdir kSNP_output_k19 -k 19 -ML -NJ -vcf -core -min_frac 0.9
	# rerunning it to output to external drive, because this takes an obscene amount of temp space.  Like over 10GB.  for 95 strains.  Why.

```


# Plotting the results, etc
The resulting files were then processed with the `plot_Clermont_comparison.R` script.
