folders=$1

for folder in $folders
do
  ./batchtest.sh "../data/${folder}/" "cg cplex" 900 > "tr${folder}900.raw"
done
