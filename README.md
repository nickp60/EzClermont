# Clermont PCR typing tool

## Description

This is a crude tool for using the Clermont 2013 PCR typing method for *in silico* analysis of *E. coli* whole genomes or assembled contigs. This has no fancy bells, whistles, mismatch handling, fuzzy searches, or any of that nonsense.  This is the work of Sunday night and Monday morning, December 11-12 2016, and whatever modification need to be done to keep it chugging along.

The concept is based on the SMS2 webtool for doing PCR product size prediction.

There are (at this time) two command line options:

```
usage: clermontPCR.py [-p] [-c] [-h] contigs

run a 'PCR' to get clermont types

positional arguments:
  contigs           FASTA formatted genome or set of contigs

optional arguments:
  -p, --partial     If scanning contigs, breaks between contigs could
                    potentially contain your sequence of interest. if
                    --partial, partial matches that could be intereupted by
                    contig breaks are reported
  -c, --no_control  ignore failure of control PCR
  -h, --help        Displays this help message
```


It prints out the presense or absence of the PCR product to stderr.  It doesnt check the length or anything fancy.  When using --partial, if a single primer has a hit but the contig starts/ends within the length of the expected product size, we call it a hit.

A minimal "filename.fasta    ClermontType" output can be piped to a results file using a shell loop:

```
for i in strain1 strain2 strain3;
	do
	python3 clermontPCR.py ${i} >> results.txt
done
```

This is provided with no guarantee that it will correlate nicely with your experimental data.  If your PCRs work nicely in the lab and disagree with what this program tells you, I would say 'bands don't lie' and move on.


Have fun!


### Testing
The script also contains two really basic tests that can be run by either unittests or nosetests.

### Requirements
Biopython
