#!/bin/bash

outfile=$1

nmin=$2
range=$3
alpha="0.5"
beta="0.2"
min=$4
max=$5
access_prob=$6

num_graphs=$7

for i in $(seq 1 $num_graphs)
do
  m=$((RANDOM%$range+$nmin))
  echo $m
seed=`date +%N`
echo "geo 1 $seed" > topo_spec
echo "$m 100 2 $alpha $beta" >> topo_spec
../gtitm/gt-itm/bin/itm topo_spec
../gtitm/gt-itm/bin/sgb2alt topo_spec-0.gb out.gtts
cat out.gtts |
awk '
/^[0-9]/ { print prefix " " $0; }
/\:/ { prefix = substr($1, 1, 1); } 
' |
gawk -v max=$max -v min=$min -v ap=$access_prob '
BEGIN { srand(systime() + PROCINFO["pid"]); rlen=max-min+1 }
function random() { return int(rand()*rlen) + min; }
function access() { return (rand() < ap) ? 1 : 0; }
/^G/ { print $0; }
/^V/ { print $1 " " $2 " " $4 " " $5 " " random() " " random() " " access();  }
/^E/ { print $1 " " $2 " " $3 " " random(); }
' > "${outfile}${i}.gtt"
rm topo_spec topo_spec-0.gb out.gtts
done
