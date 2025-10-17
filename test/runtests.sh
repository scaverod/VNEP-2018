#!/bin/bash
# Copyright 2013 Leonardo Moura

folders=$1
methods=$2
timelimit_in_s=$3

for folder in $folders
do
  echo "../data/${folder} - timelimit: ${timelimit_in_s}"
  ./batchtest.sh "../data/${folder}/" $methods $timelimit_in_s &> "t${folder}_${timelimit_in_s}.raw"
done
