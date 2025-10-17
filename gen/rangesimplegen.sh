#!/bin/bash
# Copyright 2013 Leonardo Moura

folder=$1
reps=$2
prefix=$3
minsize=$4
maxsize=$5
cmin=$6
cmax=$7

for i in $(seq 1 $reps)
  do
    for size in $(seq $minsize $maxsize)
      do
        ./simplegen.sh "${folder}/${prefix}${size}_${i}" $size 0 $cmin $cmax
      done
  done
