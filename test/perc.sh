 awk '{ 
  if (NR == 1) {
    print $0;
  } else {
  for (i =1; i < NF; i++)
    if (i == 5 || i == 11) {
      time = $15
      if (time > 0.0001)
        printf("%.2f ", ($i/1000000.0)/time*100);
      else
        printf("%.2f ", 0.0);
    } else
      printf("%s ", $i);
  print "" }
 }' |
column -t

