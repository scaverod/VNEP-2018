#!/bin/bash
# Copyright 2013 Leo Moura

reps=$1

smin=20
smax=40
vmin=10
vmax=25

for i in $(seq 1 $reps)
  do
    for size in $(seq 10 100)
      do
        ./amplgen.sh $size 0 $smin $smax $(( ($size / 2) - 2 )) 0 $vmin $vmax > "../data/ampl2/m${size}_${i}.data"
      done
  done
