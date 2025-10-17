#!/bin/bash
# Copyright 2014 Leonardo Moura
# example usage: ./% sub 10 1 50 
# --> generates a graf named sub.gtt with 10 nodes, capacities in [1,50] alpha = .5 and beta = .2

outfile=$1
n=$2
mincn=$3  # node range
maxcn=$4
mince=$5  # edge range
maxce=$6 

alpha="0.5"
beta="0.2"

./complete_simplegen.sh $outfile $n 0 $mincn $maxcn $mince $maxce $alpha $beta
