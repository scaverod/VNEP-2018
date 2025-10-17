#!/bin/bash
# Copyright 2013 Leonardo Moura
# runs the method $1 on the substrate graph in file $2 and virtual graph in file $3
# with a time limit of $4

method=$1
subfile=$2
virfile=$3
timelimit_in_s=$4

subnodes=`./gttnvertices.sh $subfile`
virnodes=`./gttnvertices.sh $virfile`

case $method in 
  "ampl")
    ../gen/gtt2mod.sh $subfile $subfile > temp.data
    ../model/solve.sh temp.data
    rm temp.data
  ;;
  "backtrack"|0|2|3|4|5|6|7|8|9|10|11|12|13|14|15)
    if [ "$method" = "backtrack" ]
    then
      method=0
    fi
    ./genparamfile.sh ../vne/params.ini substrate_file $subfile virtual_file $virfile timelimit_in_s $timelimit_in_s
    /usr/bin/time -f "time_back${method}:%e:${subnodes}_${virnodes}" ../vne/bin/vne.out $method tparam.ini 2>&1 |
    awk -v suffix="${subnodes}_${virnodes}" '/^@/ { print $0 ":" suffix; } !/^@/ { print $0; }'
    rm -f tparam.ini /tmp/output
  ;;
  "cplex")
    ./genparamfile.sh ../vne/params.ini substrate_file $subfile virtual_file $virfile timelimit_in_s $timelimit_in_s
    /usr/bin/time -f "time_cplex:%e:${subnodes}_${virnodes}" ../vne/bin/vne.out 1 tparam.ini 2>&1 |
    awk -v suffix="${subnodes}_${virnodes}" ' /^@/ { print $0 ":" suffix; } !/^@/ { print $0; }'
    rm -f tparam.ini /tmp/output
  ;;
  "bb")
    ./genparamfile.sh ../vne/params.ini substrate_file $subfile virtual_file $virfile timelimit_in_s $timelimit_in_s
    /usr/bin/time -f "time_bb:%e:${subnodes}_${virnodes}" ../vne/bin/vne.out 2 tparam.ini 2>&1 
    rm -f tparam.ini
  ;;
  "cg")
    method=10
    ./genparamfile.sh ../vne/params.ini substrate_file $subfile virtual_file $virfile timelimit_in_s $timelimit_in_s
    /usr/bin/time -f "time_cg:%e:${subnodes}_${virnodes}" ../vne/bin/vne.out $method tparam.ini 2>&1 |
    awk -v suffix="${subnodes}_${virnodes}" ' /(Cost)|^@/ { print $0 ":" suffix; } !/^@/ { print $0; }'
    rm -f tparam.ini /tmp/output
  ;;
esac
