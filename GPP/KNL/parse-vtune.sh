#!/bin/bash 
# 
# Parses a VTune summary report for uncore memory access counts

module load python
python ./parse-vtune.py
