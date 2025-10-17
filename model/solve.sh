#!/bin/sh

file=$1

echo "model embed.mod;" > run.ampl
echo "data ${file};" >> run.ampl
echo "option solver cplexamp;" >> run.ampl
echo "option send_statuses 0;" >> run.ampl
echo "solve;" >> run.ampl
echo "display solve_result;" >> run.ampl
echo "display Cost;" >> run.ampl

/usr/bin/time -f %U -o time.txt /opt/ampl/bin/ampl run.ampl | sed 's/\:/\.\./g'

echo "time:"`cat time.txt`":"`echo $file | sed 's/\_.*//g'`
#rm -f run.ampl time.txt
