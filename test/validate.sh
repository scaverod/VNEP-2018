#!/bin/bash
# Copyright 2013 Leonardo Moura
# Run methods $2 and $3 on all combinations of files sub*.gtt, vir*.gtt in the folder $1
#     fail if result of $2 is greater than result of $3
# eg: ./validate.sh ../data/large cg cplex

folder=$1
method1=$2
method2=$3
timelimit_in_s=$4

for subfile in $(ls $folder*sub*)
do
  for virfile in $(ls $folder*vir*)
  do
    subnodes=`./gttnvertices.sh $subfile`
    virnodes=`./gttnvertices.sh $virfile`

    if [ $(( $subnodes / 2 )) -ge "$virnodes" ] 
      then
        echo $subfile $virfile
        cost2=`./test.sh $method1 $subfile $virfile $timelimit_in_s | egrep Cost | awk 'BEGIN{ FS=":"} { print $2;}'`
        cost3=`./test.sh $method2 $subfile $virfile $timelimit_in_s | egrep Cost | awk 'BEGIN{ FS=":"} { print $2;}'`
        printf "%i == %i .... " $cost2 $cost3
        if [ "$cost2" -gt "$cost3" ]; then
          echo "erro"
          echo $subfile $virfile
        else
          echo "OK"
        fi
      fi
  done
done

