cat $1 | awk ' /^E/ { sum += $4;} END{print sum;}'
