#!/bin/bash
# Copyright 2013 Leonardo Moura

outfile=$1
size=$2
min=$3
max=$4

alpha="0.5"
beta="0.2"

seed=`date +%N`
echo "hier 1 $seed" > topo_spec
echo "2 0 0" >> topo_spec
echo "$size 100 1 $alpha $beta" >> topo_spec
echo "5 10 1 $alpha $beta" >> topo_spec

cat topo_spec

./itm topo_spec
./sgb2alt topo_spec-0.gb out.gtts
cat out.gtts |
awk '
/^[0-9]/ { print prefix " " $0; }
/\:/ { prefix = substr($1, 1, 1); } 
' |
gawk -v max=$max -v min=$min '
BEGIN { srand(systime() + PROCINFO["pid"]); rlen=max-min+1 }
function random() { return int(rand()*rlen) + min; }
/^G/ { print $0; }
/^V/ { print $1 " " $2 " " random();  }
/^E/ { print $1 " " $2 " " $3 " " random(); }
' > "${outfile}${i}.gtt"
rm topo_spec topo_spec-0.gb out.gtts
