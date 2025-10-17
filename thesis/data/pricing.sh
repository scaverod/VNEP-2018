
for sset in sparse dense hier ts
do
  for vset in vir1 vir4
  do
    echo $sset
    cat "pricingmaster${sset}_3600.raw" | 
    egrep $vset |
    tformat --sum @primal_time13  @addc_time13  @aux_time13 @pricing_time13 | 
    ./../../vne/format/sort.sh | 
    awk '
      BEGIN { a = b = c = d = 0.0; }
      { if(NR == 1) print $0; 
        else if (int($3) > 2 && int($1) > 20) {
          printf("%s \t & %s \t & %.2f \t & %.2f \t & %.2f \t & %.2f \\\\ \n", $1, $3, $4 * 100, $5 * 100, $6 * 100, $7 * 100);
          a += $4;
          b += $5;
          c += $6;
          d += $7;
        }}
        END { printf("a & b \t & %.2f \t & %.2f \t & %.2f \t & %.2f \\\\ \n ", (a * 100) / (NR-1), (b * 100) / (NR-1), (c * 100) / (NR-1), (d * 100) / (NR-1));}
    '  | 
    column -t
  done
done
