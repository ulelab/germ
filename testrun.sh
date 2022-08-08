#!/bin/bash

Rscript germs.R -f test_data/test.fasta -o output.tsv.gz -w 122 -s 122 -t transcripts.txt -p plots
