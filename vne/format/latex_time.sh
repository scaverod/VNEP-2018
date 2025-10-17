#!/usr/bin/env bash
time_limit=$1

awk -v t=$time_limit '
  /@time/ { if (int($3) >= t) $3 = t;print; print $1"_tl:" ((int($3) >= t) ? 1 : 0) " : "$5 " " $6 " " $7; }
  !/@time/ { print;} ' |
tformat --latex --mean {@time,@firstint,@cost}{1,13} --sum "@time"{1,13,14,15,16}"_tl" |
./format/sort.sh
