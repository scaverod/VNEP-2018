#!/bin/bash
# Copyright 2013 Leonardo Moura

folder=$1
reps=$2
prefix=$3
range=$4
cmin=$5
cmax=$6
topology=$7 # either simple or hier or ts

for i in $(seq 1 $reps)
  do
    for size in $range
      do
        case $topology in
          "simple") ./simplegen.sh "${folder}/${prefix}${size}_${i}" $size $cmin $cmax ;;
          "hier") ./hiergen.sh "${folder}/${prefix}${size}_${i}" $size $cmin $cmax ;;
          "ts") ./tsgen.sh "${folder}/${prefix}${size}_${i}" $size $cmin $cmax ;;
        esac
      done
  done
