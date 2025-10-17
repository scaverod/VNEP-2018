#!/usr/bin/env bash

sed 's/vir[0-9]//g' |
sort -nk 2 |
sort -ns |
sed "s/\([0-9]\+\)/\1 \&/"
