#!/bin/bash

echo "not tested"
exit 1

count=0
for sub in 1 2 3 4 5
do
  for sim in 2 3 
  do
    count=$(( $count + 1 ))
for tmin in 1 10 100 1000
  do
  for tmax in 100 1000 10000
    do
    if [ "$tmin" -lt "$tmax" ]
    then
      ./genparamfile.sh tmin $tmin tmax $tmax substrate_file "../data/sub${sub}.gtt" virtual_file_prefix "../data/sim${sim}/virtual" > tparam.ini
      #cost="${RANDOM}"
      cost=`./build/vne.out tparam.ini | tail -n 1 | cut -f 2 -d ":"`
      #echo "${tmin}${tmax}:${RANDOM}:${count}"
      echo "${tmin},${tmax}:${cost}:${count}"
      rm -rf tparam.ini
    fi
  done
done
done
done

