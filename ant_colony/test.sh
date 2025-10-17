#!/bin/bash

function printRes() {
  echo "${name}:${cost}:${count}"
}

count=0
for sub in 1 2 3 4 5
do
  for sim in 2 3 
  do
    count=$(( $count + 1 ))
    #cost="${RANDOM}"
    ./genparamfile.sh substrate_file "../data/sub${sub}.gtt" virtual_file_prefix "../data/sim${sim}/virtual"
    cost=`./build/vne.out tparam.ini | tail -n 1 | cut -f 2 -d ":"`
    #echo "${tmin}${tmax}:${RANDOM}:${count}"
    name="complete"
    printRes

    ./genparamfile.sh ants 1 maxiter 1 substrate_file "../data/sub${sub}.gtt" virtual_file_prefix "../data/sim${sim}/virtual"
    cost=`./build/vne.out tparam.ini | tail -n 1 | cut -f 2 -d ":"`
    name="greedy"
    printRes

    ./genparamfile.sh tmin 0 tmax 50202185000 substrate_file "../data/sub${sub}.gtt" virtual_file_prefix "../data/sim${sim}/virtual"
    cost=`./build/vne.out tparam.ini | tail -n 1 | cut -f 2 -d ":"`
    name="nomaxmin"
    printRes

    ./genparamfile.sh onlythebestreinforces false  substrate_file "../data/sub${sub}.gtt" virtual_file_prefix "../data/sim${sim}/virtual"
    cost=`./build/vne.out tparam.ini | tail -n 1 | cut -f 2 -d ":"`
    name="allants"
    printRes

    rm -rf tparam.ini
  done
done

