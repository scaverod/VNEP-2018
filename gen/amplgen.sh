#!/bin/bash

snmin=$1
srange=$2
smin=$3
smax=$4

vnmin=$5
vrange=$6
vmin=$7
vmax=$8

./simplegen.sh sub $1 $2 $3 $4 1
./simplegen.sh vir $5 $6 $7 $8 1

./gtt2mod.sh sub1.gtt vir1.gtt
rm sub1.gtt vir1.gtt
