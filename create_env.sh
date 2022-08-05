#!/bin/bash

# Script to create conda environment and install germs
# A. M. Chakrabarti
# 5th August 2022

conda env create -f environment.yml
conda run -n germs bash -c "R -e 'devtools::install()'"