# Getting the validation data
The Strains listed in Clermont 2015  are mostly from bioproject PRJNA321606 (all the ECOR strains), and a handful of others.  The other strains were identified manually by seaching NCBI; 6 strains from the paper were not availible on Genbank. The strains from sims2011 were were identified similarly, and added to the `validation_metadata.csv` sheet.  The accession from this used with the Entrez Batch tool to download; download the Genbank ones rather than the refeq, and refseq is missing a number of the accessions.

Additionally, as the recent publication about the G phylogroup has changed the phylotyping scheme.  112 strains (not our 112 strains listed above, coincidentally)  were included in their annalysis; not all were sequenced, and none included accessions except for the 13 SRA accessions listed here.  These strains in the table below were downlaoded from enterobase and added to the genomes directory.

| Name    | Barcode         | Accession    |
|---------|-----------------|--------------|
| eo1696  | ESC_BA6587AA_AS | clermont2019 |
| eo916   | ESC_CA3101AA_AS | clermont2019 |
| eo669   | ESC_CA2711AA_AS | clermont2019 |
| ECO0122 | ESC_BA1311AA_AS | clermont2019 |
| eo3196  | ESC_BA6371AA_AS | clermont2019 |
| eo2375  | ESC_BA8369AA_AS | clermont2019 |
| eo2659  | ESC_AA8494AA_AS | clermont2019 |
| eo1698  | ESC_AA9527AA_AS | clermont2019 |
| eo464   | ESC_AA8884AA_AS | clermont2019 |
| eo230   | ESC_BA9755AA_AS | clermont2019 |
| eo2839  | ESC_EA6188AA_AS | clermont2019 |
| ECO0183 | ESC_BA0658AA_AS | clermont2019 |
| eo3033  | ESC_BA0319AA_AS | clermont2019 |

# Running EzClermont and ClermonTyping
The strains were unzipped if necessary and were then added to a `./genome_assemblies_genome_fasta/tmp/` directory.  Then EzClermont and ClermonTyping (git version 5ae1a2ba) was run:

```
mkdir ./genome_assemblies_genome_fasta/tmp/
gunzip ./genome_assemblies_genome_fasta/ncbi-genomes-2020-04-06/GCA_00*
for i in ./genome_assemblies_genome_fasta/ncbi-genomes-2020-04-06/*.gz ; do bas=$(filename $i) ; cat $i | gunzip > ./genome_assemblies_genome_fasta/tmp/${bas}.fasta ; done
mv ~/Downloads/ESC_* ./genome_assemblies_genome_fasta/tmp/


time for i in ./genome_assemblies_genome_fasta/tmp/*.fasta ; do nam=$(basename $i) ; cat $i ezclermont  - -e $nam >> 2020-04-07-ezclermont.txt; done
# 3m57.222s
mkdir ct_results ; cd ct_results
time for i in ../genome_assemblies_genome_fasta/tmp/*.fasta; do   ~/GitHub/ClermonTyping/clermonTyping.sh --fasta $i ; done
# real    3m12.872s
cat analysis_2020-04-07_15*/*phylogroups.txt > 2020-04-07-CT-results.txt


```



# Generating kmer-base phylogentic tree


```
/Applications/Harvest-OSX64-v1.1.2/parsnp -c -d ./genome_assemblies_genome_fasta/tmp/ -r !  -p 4 -o ./alignment/
/Applications/Harvest-OSX64-v1.1.2/harvesttools  -i ./alignment/parsnp.ggr -M  ./alignment/parsnp.msa

cat alignment/parsnp.msa | sed -e 's/^>\(.\{4\}\)\(.\{9\}\).*/>\2/'  > alignment/parsnp.clean.msa
# from commit 65129b2
python ~/GitHub/open_utils/frommsa/frommsa.py ~/GitHub/ezclermont/docs/analysis/validate/alignment/parsnp.clean.msa  > ~/GitHub/ezclermont/docs/analysis/validate/alignment/parsnp.clean.phy

# version 3.3.20190909
phyml -i ./alignment/parsnp.clean.phy~
```

# Plotting the results, etc
The resulting files were then processed with the `plot_Clermont_comparison.R` script.


# Addressing inconsistenties
## ECOR-46
As both ClermonTyping and EzClermont typed GCA_002191015.1 as "G", I checked whether the assembly could be at fault.  I found an additional assembly of ECOR46 in enterobase, ESC_HA7443AA_AS, which is the result of three SRAs: SRR3951511, SRR3987677, SRR4099204. This assembly types as a "D" strain.  Further, GCA_002191015.1 appears to have very poor assembly quality, with 2145 contigs falling below the 500bp default threshold.  It was based an another SRA, SRR3886568.

Aligning the two assemblies with Mauve showed these to differ quite strongly.  The poor assembly quality makes it difficult to confirm, but the two assemblies do not appear to represent the same organism.


SRR3886568 was downloaded

```
docker run --memory 10G --rm  -v $PWD:/inputdata/ nickp60/ezblobtools  -r /inputdata/BA000007.2.fasta -F /inputdata/ECOR46/SRR3951511_1.fastq   -d ref_prok_rep_genomes -o /inputdata/ECOR46/SRR3951511_1_blob/ -t 2 -m 10

docker run --memory 10G --rm  -v $PWD:/inputdata/ nickp60/ezblobtools  -r /inputdata/BA000007.2.fasta -F /inputdata/ECOR46/SRR3886568_1.fastq   -d ref_prok_rep_genomes -o /inputdata/ECOR46/SRR3886568_1_blob/ -t 2 -m 10
```

The blobplot revealed a clear band of *E. coli* with coverage differentfrom teh bulk of the genome.  This is indicative of a mixed culture.  Interestingly, running either EzClermont or ClermonTyping  on the re-assembled genome (generated as part of the BlobTools scheme) types the strain according to the phylogeny, as group F.

## Ecor44

We downloaded SRR3885711
```
ECOR="ECOR44"; SRA="SRR3885711"; rm -r ./${ECOR}/${SRA}_blob/; docker run --memory 10G --rm  -v $PWD:/inputdata/ nickp60/ezblobtools  -r /inputdata/BA000007.2.fasta -F /inputdata/${SRA}_1.fastq   -d ref_prok_rep_genomes -o /inputdata/${ECOR}/${SRA}_blob/ -t 2 -m 10
```

Typing the reassembled genome still gave "G", indicating a missmatch in the arpA allele. simpleOrtho was used to pull out that gene from the genome:
```
simpleOrtho.py -d ./tmp/ -i ~/GitHub/ezclermont/docs/analysis/training/refs/aceK_arpA.fasta --nuc -o tmp_so | extractRegion - -f ./tmp/contigs.fasta
```
Alignment with the training data revealed a mismatch in the reverse primer, but in the final 5 bases (3'), a region excluded from degeneracy as this is used to differentiate alleles. 9 other strains had the same:
ESC_FA0109AA
ESC_GA9273AA
ESC_IA1730AA
ESC_EA2533AA
ESC_FA0386AA
ESC_EA2560AA
ESC_IA4116AA
ESC_GA3114AA
ESC_EA4735AA


## ECOR04

SRR3823860
```
ECOR="ECOR04"; SRA="SRR3823860"; rm -r ./${ECOR}/${SRA}_blob/; docker run --memory 10G --rm  -v $PWD:/inputdata/ nickp60/ezblobtools  -r /inputdata/BA000007.2.fasta -F /inputdata/${SRA}_1.fastq   -d ref_prok_rep_genomes -o /inputdata/${ECOR}/${SRA}_blob/ -t 2 -m 10
```
Reassembled genome types out as a A.

## ECOR51
SRR3921609
```
ECOR="ECOR51"; SRA="SRR3921609"; rm -r ./${ECOR}/${SRA}_blob/; docker run --memory 10G --rm  -v $PWD:/inputdata/ nickp60/ezblobtools  -r /inputdata/BA000007.2.fasta -F /inputdata/${SRA}_1.fastq   -d ref_prok_rep_genomes -o /inputdata/${ECOR}/${SRA}_blob/ -t 2 -m 10
```
## ECOR70
SRR3923908
```
ECOR="ECOR70"; SRA="SRR3923908"; rm -r ./${ECOR}/${SRA}_blob/; docker run --memory 10G --rm  -v $PWD:/inputdata/ nickp60/ezblobtools  -r /inputdata/BA000007.2.fasta -F /inputdata/${SRA}_1.fastq   -d ref_prok_rep_genomes -o /inputdata/${ECOR}/${SRA}_blob/ -t 2 -m 10
```
## ECOR72
SRR3923910
```
ECOR="ECOR72"; SRA="SRR3923910"; rm -r ./${ECOR}/${SRA}_blob/; docker run --memory 10G --rm  -v $PWD:/inputdata/ nickp60/ezblobtools  -r /inputdata/BA000007.2.fasta -F /inputdata/${SRA}_1.fastq   -d ref_prok_rep_genomes -o /inputdata/${ECOR}/${SRA}_blob/ -t 2 -m 10
```
