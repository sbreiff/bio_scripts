#!/bin/bash

################################################################################
#                                                                              #
#  script for using TrimAl to generate several different trimmed alignments    #
#                                                                              #
#  first argument is input file (aligned), second argument is output prefix    #
#                                                                              #
#  Example Usage:                                                              #  
#  $ bash trim_alignment.sh Aur_for_tree2aligned.fas Aur2                      #
#                                                                              #
################################################################################


~/bin/trimal -in $1 -out $2'_gappyout.fasta' -gappyout

~/bin/trimal -in $1 -out $2'_auto.fasta' -automated1

~/bin/trimal -in $1 -out $2'_strict.fasta' -strict

~/bin/trimal -in $1 -out $2'_strictplus.fasta' -strictplus

~/bin/trimal -in $1 -out $2'_gt75.fasta' -gt 0.75

~/bin/trimal -in $2'_gt75.fasta' -out $2'_trim1.fasta' -st 0.001 -cons 75
