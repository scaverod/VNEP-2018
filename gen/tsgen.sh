#!/bin/bash
# Copyright 2014 Leonardo Moura

outfile=$1 # without .gtt
n=$2
min=$3
max=$4

alpha="0.4"

seed=`date +%N`
tsize=`echo "100 / ${n}" | bc -l`
echo "ts 1 $seed" > topo_spec
echo "4 0 0" >> topo_spec  # stub domains / transit node
echo "${n} ${tsize} 3 ${alpha}" >> topo_spec # transit domains
echo "${n} ${tsize} 3 ${alpha}" >> topo_spec # nodes/transit domain
echo "6 10 3 ${alpha}" >> topo_spec # nodes/stub domains
./itm topo_spec
./sgb2alt topo_spec-0.gb out.gtts
cat out.gtts |
awk '
/^[0-9]/ { print prefix " " $0; }
!/^[ 0-9]/ { prefix = substr($1, 1, 1); } 
' |
gawk -v max=$max -v min=$min '
BEGIN { srand(systime() + PROCINFO["pid"]); rlen=max-min+1 }
function random() { return int(rand()*rlen) + min; }
/^G/ { print $0; }
/^V/ { print $1 " " $2 " " random();  }
/^E/ { print $1 " " $2 " " $3 " " random(); }
' > "${outfile}${i}.gtt"
rm topo_spec topo_spec-0.gb out.gtts
