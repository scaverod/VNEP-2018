#!/bin/bash
# Copyright 2014 Leonardo Moura
# example usage: ./% sub 10 1 50 
# --> generates a graf named sub.gtt with 10 nodes, capacities in [1,50] alpha = .5 and beta = .2

outfile=$1
n=$2
minc=$3
maxc=$4

alpha="0.5"
beta="0.2"

./complete_simplegen.sh $outfile $n 0 $minc $maxc $minc $maxc $alpha $beta
