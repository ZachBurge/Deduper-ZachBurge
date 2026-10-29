#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8

/usr/bin/time -v python ../burge_deduper.py -u ../STL96.txt -f ../C1_SE_uniqAlign.sorted.sam -o ../output.sam