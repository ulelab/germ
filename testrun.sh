#!/bin/bash

# Test run for standard sequences
Rscript germs.R -f test_data/test.fasta -o output.tsv.gz -w 122 -s 122 -t transcripts.txt -p plots

# Test run for sequences non-standard characters
Rscript germs.R -f test_data/test_nonstdchars.fasta -o output_nonstdchars.tsv.gz -w 122 -s 122 -t transcripts_nonstdchars.txt -p plots_nonstdchars

# Test run for standard sequences (no output filename provided)
Rscript germs.R -f test_data/test.fasta -w 122 -s 122 -t transcripts.txt