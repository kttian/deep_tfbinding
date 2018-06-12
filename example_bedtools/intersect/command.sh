#!/usr/bin/env bash

#Get the core regions that are positives for task 1:
#-wa writes out the entry from file "a" in the output
#-u makes sure that each entry from file "a" is written out at most once,
#even if it overlaps with multiple entries from file "b"
#-f 0.5 means that the overlap fraction must be atleast 50%
bedtools intersect -a core_regions.bed -b task1_peaks.bed -wa -u -f 0.5


#Get the core regions that are positives for task 2:
bedtools intersect -a core_regions.bed -b task2_peaks.bed -wa -u -f 0.5
