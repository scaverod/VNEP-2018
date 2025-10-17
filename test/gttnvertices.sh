#!/bin/bash
# Copyright 2013 Leonardo Moura
# outputs how many nodes a graph described in the .gtt file $1 has
# eg: ./gttnvertices.sh ../data/simple/sub1.gtt

gttfile=$1
cat $gttfile | egrep "^G" | cut -f 2 -d " "
