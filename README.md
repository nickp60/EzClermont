[![Build Status](https://travis-ci.org/nickp60/clermontpcr.svg?branch=master)](https://travis-ci.org/nickp60/clermontpcr.svg?branch=master)[![Coverage Status](https://coveralls.io/repos/github/nickp60/clermontpcr/badge.svg?branch=master)](https://coveralls.io/github/nickp60/clermontpcr?branch=master)[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/clermontpcr.svg)](https://badge.fury.io/py/clermontpcr)

![Icon](https://github.com/nickp60/clermontpcr/blob/master/icon/clermontPCR-0.png)
# Clermont PCR typing tool

## Description

This is a crude tool for using the Clermont 2013 PCR typing method for *in silico* analysis of *E. coli* whole genomes or assembled contigs. This has no fancy bells, whistles, mismatch handling, fuzzy searches, or any of that nonsense.  This is the work of Sunday night and Monday morning, December 11-12 2016, and whatever modification need to be done to keep it chugging along.

I updated on August 2, 2017 to add reactions that differentiate A/C, D/E/cryptic, and to add more robust tests.

The concept is loosely based on the SMS2 webtool for doing PCR product size prediction.

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


It prints out the presense or absence of the PCR product to stderr, and the resulting phylotype to stdout.  It checks the length, accepting fragments that are within 20bp of the expected size.  When using --partial, if a single primer has a hit but the contig starts/ends within the length of the expected product size, we call it a hit.

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
The tests can be run by either unittests or nosetests.  I dont have a validating C phylotype example, due to some [discrepancies](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0105395) noted in the literature

### Requirements
Biopython
