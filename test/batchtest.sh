#!/bin/bash
# Copyright 2013 Leonardo Moura
# Run the methods in $2 on all combinations of files sub*.gtt, vir*.gtt in the folder $1 with a time limit of $3
# eg: ./batchtest.sh ../data/large "backtrack cplex" 100

folder=$1
methods=$2
timelimit_in_s=$3

for subfile in $(ls $folder*sub*)
do
  for virfile in $(ls $folder*vir*)
  do
    subnodes=`./gttnvertices.sh $subfile`
    virnodes=`./gttnvertices.sh $virfile`

    if [ $(( $subnodes / 2 )) -ge "$virnodes" ] 
      then
        echo $subfile $virfile
        for method in $methods
          do
            ./test.sh $method $subfile $virfile $timelimit_in_s
          done
      fi
  done
done

