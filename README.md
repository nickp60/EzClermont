[![Build Status](https://travis-ci.org/nickp60/clermontpcr.svg?branch=master)](https://travis-ci.org/nickp60/EzClermont.svg?branch=master)[![Coverage Status](https://coveralls.io/repos/github/nickp60/EzClermont/badge.svg?branch=master)](https://coveralls.io/github/nickp60/EzClermont?branch=master)[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/clermontpcr.svg)](https://badge.fury.io/py/ezclermont)

![Icon](https://github.com/nickp60/EzClermont/blob/master/icon/clermontPCR-0.png)
# EzClermont: The *E. coli* Clermont PCR phylotyping tool

## Description

This is a tool for using the Clermont 2013 PCR typing method for *in silico* analysis of *E. coli* whole genomes or assembled contigs.

### Changelog
 - bump to version 0.4 in May 2018; improved handling of partial matches
 - made a webapp on April 19th, 2018 after requests from several to make the tool more user friendly.
 - updated on August 2, 2017 to add reactions that differentiate A/C, D/E/cryptic, and to add more robust tests.
 - released Dec. 2016


## Usage
There are (at this time) two command line options:

```
usage: ezclermont [-p] [-c] [-h] contigs

run a 'PCR' to get clermont types

positional arguments:
  contigs           FASTA formatted genome or set of contigs

optional arguments:
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        minimum contig length to consider.default: 500
  -n, --no_partial      If scanning contigs, breaks between contigs could
                        potentially contain your sequence of interest. if
                        --no_partial, these plausible partial matches will NOT
                        be reported; default behaviour is to consider partial
                        hits if the assembly has more than 4 sequnces(ie, no
                        partial matches for complete genomes, allowing for 1
                        chromasome and several plasmids)
  -h, --help            Displays this help message
```


It prints out the presense or absence of the PCR product to stderr, and the resulting phylotype to stdout.  It checks the length, accepting fragments that are within 20bp of the expected size.  When using --partial, if a single primer has a hit but the contig starts/ends within the length of the expected product size, we call it a hit.

A minimal "filename.fasta    ClermontType" output can be piped to a results file using a shell loop:

```
for i in strain1 strain2 strain3;
	do
	  ezclermont ${i} >> results.txt
done
```
or, using GNU parallel, and saving a log file:
```
ls ./folder/with/assemblies/*.fa | parallel "ezclermont {} 1>> results.txt  2>> results.log"
```


Have fun!


### Testing
The tests can be run by either unittests or nosetests.

### Requirements
#### commandline tool
Biopython
#### webapp
flask
biopython


## Acknowledgements
Thanks to [Dave Gamache]( https://github.com/dhg/Skeleton) for Skeleton, the webapp CSS theme.

## Name note
The name of this repo (and pypi package was changed on April 21 from ClermontPCR to EzClermont.
