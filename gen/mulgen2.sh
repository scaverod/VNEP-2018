#!/bin/bash
# Copyright 2013 Leonardo Moura

folder=$1
reps=$2
prefix=$3
range=$4
cmin=$5 # node range
cmax=$6
emin=$7 # edge range
emax=$8
topology=$9 # either simple or hier or ts

for i in $(seq 1 $reps)
  do
    for size in $range
      do
        case $topology in
          "simple") ./simplegen2.sh "${folder}/${prefix}${size}_${i}" $size $cmin $cmax $emin $emax ;;
          #"hier") ./hiergen.sh "${folder}/${prefix}${size}_${i}" $size $cmin $cmax ;;
          "ts") ./tsgen2.sh "${folder}/${prefix}${size}_${i}" $size $cmin $cmax $emin $emax ;;
        esac
      done
  done
