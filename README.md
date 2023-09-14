# GeRM - Generalised RNA multivalency

**[Mutual homeostasis of charged proteins](https://doi.org/10.1101/2023.08.21.554177)**

Rupert Faraway, Neve Costello Heaven, Holly Digby, Oscar G Wilkins, Anob M Chakrabarti, Ira A Iosub, Lea Knez, Stefan L Ameres, Clemens Plaschka, Jernej Ule

bioRxiv (2023) https://doi.org/10.1101/2023.08.21.554177

## Table of contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Testing](#testing)
4. [Quickstart](#quickstart)
5. [Parameters](#parameters)

## Introduction

GeRM is a command-line tool written in R and Rcpp to calculate **Ge**neralised **R**NA **M**ultivalency Scores for user-supplied sequences. The custom functions are contained within the GeRM R package.

### The algorithm

GeRM is calculated from a string of consecutive overlapping nucleotide sequences of length
_k_ (_k_-mers).

In non-mathematical terms, the GeRM score is calculated by comparing a k-mer to all the other k-mers that surround it in a fixed window. For each of the surrounding k-mers, the sequence similarity to the central k-mer is calculated from the negative exponent of the Hamming distance, such that k-mers with identical sequences have a high score and those with unrelated sequences have a low score. The constant λ determines how quickly this similarity score decays as sequences become more dissimilar to the central k-mer. This sequence similarity score is multiplied by a distance score, which decays linearly from 1 to 0 with distance from the central k-mer. k-mers that overlap with the central k-mer are ignored. For k-mers at the edges of transcripts, where the window exceeds the end of the transcript, all positions that fall outside of the transcript are given a score of 0. The sum of all the distance-weighted sequence similarities is summed to give the GeRM score.

For more details please see the _Methods_ section of the manuscript.

## Installation

Installation has been tested using Mac OS Ventura 13.5.2 and CentOS Linux 7 Core. GeRM requires only a standard desktop computer or laptop, but the RAM requirements will scale with the lenth of the input sequences.

To install GeRM, first clone the repository to your local computer with:
```
git clone https://github.com/ulelab/germ.git
```

Then, there are two options for installing the dependencies.

### 1. Conda option (recommended)

If you have Conda on your system you can create a virtual environment which installs R and all the dependencies using the provided YAML. First move into the directory into which you cloned GeRM and then run:

```
bash create_env.sh
```

You can then activate the environment using:
```
conda activate germs
```

The environment should take approximately 10 minutes to build.

### 2. R option

GeRM requires R to be installed on your system (tested with R 4.1.2 and 4.2.0) and uses some R (`optparse`, `devtools`, `data.table`, `tidyverse`, `scales`, `ggthemes`, `cowplot`, `patchwork`, `logger`) and Bioconductor packages (`Biostrings`). If you have R already installed, you can install the GeRM R package by moving to the directory into which you cloned GeRM and then run:

```
R -e 'devtools:install()'
```

### 3. Docker option

We will soon have a GeRM Docker container available for use.

## Testing

To test the installation has worked you can run the test script. This runs three sets of GeRM test for different sequences and parameters:

```
bash testrun.sh
```

## Quickstart

GeRM can be run from the command line using:

```
Rscript germs.R --help
```

This will output the help for all the parameters that can be supplied to GeRM. The minimum is to provide a FASTA file with sequences for which to calculate GeRM scores (`--fasta`, `-f`)

## Parameters

### Basic

* `--fasta` or `-f` is used to supply the input FASTA file with the sequences for which GeRM scores will be calculated.

* `--k_length` or `-k` is used to supply the _k_-mer length with which to assess multivalency (default: 5).

* `--window_size` or `-w` is used to supply the window size for calculating multivalency (default: 123).

* `--smoothing_size` or `-s` is used to supply the smoothing window size (default: 123).

* `--output` or `-o` is used to supply the output TSV filename. If one is not supplied, then it is generated using the fasta filename, _k_-mer length, window size and smoothing window size.

### Customise GeRM calculation

* `--lambda` is used to supply the lamba value for exponential decay scaling (default: 1).

* `--scaling_function` is used to to supply a custom scaling function.

### Visualisation

* `--transcripts` or `-t` is used to provide either a comma-separated list of sequence names or a text file with one sequence name per line to plot.

* `--plot_folder` or `-p` is used to specify the folder in which to output the plots (default: plots).

### Other

* `--cores` or `-c` is used to specify the number of cores to use for parallel processing (default: 4).

* `--logging` or `-l` is used to specify the level of logging (default: INFO).
