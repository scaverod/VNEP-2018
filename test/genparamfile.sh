#!/bin/bash
# Copyright 2013 Leonardo Moura
# this script generates a file named tparam.ini from a template file
# replacing values
# example usage:
#   ./genparamfile.sh ../vne/params.ini virtual_file vir2.gtt iterations 12

genparamfile=$1
shift 1

trap "rm tparam2.ini" SIGHUP SIGINT SIGTERM

cat $genparamfile > tparam.ini
while [ "$1" ]
  do
    cat tparam.ini |
    awk -v var=$1 -v value=$2 '{ if( match($0, var) > 0) { print var " = " value; } else { print $0; } }'  > tparam2.ini
    cp tparam2.ini tparam.ini
    rm tparam2.ini
    shift 2;
  done
