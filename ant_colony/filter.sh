#!/bin/bash 
while read LINE
do
  for filter in $@
  do  
    echo $LINE | egrep "^\["$filter"\]" #| sed 's/\[.*[^ ]\]//g'
  done
done
