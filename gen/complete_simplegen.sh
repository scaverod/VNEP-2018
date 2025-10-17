#!/bin/bash
# Copyright 2014 Leonardo Moura
# example usage: ./% sub 10 10 1 50 1 100 .5 .2
# --> generates a graf named sub.gtt with 10 to 20 nodes, node capacities in [1,50], and edge capacities in [1,100] with alpha = .5 and beta = .2

outfile=$1 # sans .gtt
nmin=$2
range=$(( $3 + 1 ))
mincapn=$4 # node capacities
maxcapn=$5
mincape=$6 # edge capacities
mincape=$7
alpha=$8
beta=$9

n=$((RANDOM%$range+$nmin))
seed=`date +%N`
echo "geo 1 $seed" > topo_spec
echo "$n 100 1 $alpha $beta" >> topo_spec
./itm topo_spec
./sgb2alt topo_spec-0.gb out.gtts
cat out.gtts |
gawk '
/^[0-9]/ { print prefix " " $0; }
/\:/ { prefix = substr($1, 1, 1); } 
' |
gawk -v maxn=$maxcapn -v minn=$mincapn -v maxe=$maxcape -v mine=$mincape '
BEGIN { 
  srand(systime() + PROCINFO["pid"]); 
  nlen=maxn-minn+1
  elen=maxe-mine+1
}
function n_random() { return int(rand()*nlen) + minn; }
function e_random() { return int(rand()*elen) + mine; }
/^G/ { print $0; }
/^V/ { print $1 " " $2 " " n_random();  }
/^E/ { print $1 " " $2 " " $3 " " e_random(); }
' > "${outfile}${i}.gtt"
rm topo_spec topo_spec-0.gb out.gtts
