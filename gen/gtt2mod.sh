#!/bin/bash

subfile=$1
virfile=$2

printf "set VV := "
cat $virfile | gawk ' /^V/ { printf("%s ", $2); }' 
echo ";"

printf "set VS := "
cat $subfile | gawk ' /^V/ { printf("%s ", $2); }' 
echo ";"

echo "param cv := "
cat $virfile | gawk ' /^V/ { printf("%s %s\n", $2, $3); }' 
echo ";"

echo "param cs := "
cat $subfile | gawk ' /^V/ { printf("%s %s\n", $2, $3); }' 
echo ";"

cat $virfile |
gawk '
/^G/ { num_vertices=int($2);  
  for(i = 0; i < num_vertices; i++) 
    for(j = 0; j < num_vertices; j++)
      e[i, j] = 0;
    }
/^E/ {
  e[$2, $3] = $4;
}
END{ 
  printf("param bv:\t");
  for(i = 0; i < num_vertices; i++)
    printf(" %s", i);
  print " := ";

  for(i = 0; i < num_vertices; i++) {
    printf("\t\t%s ", i);
    for(j = 0; j < num_vertices; j++)
      printf("%s ", e[i, j]);
    print "";
  }
  print ";";
}'

cat $subfile |
gawk '
/^G/ { num_vertices=int($2);  
  for(i = 0; i < num_vertices; i++) 
    for(j = 0; j < num_vertices; j++)
      e[i, j] = 0;
    }
/^E/ {
  e[$2, $3] = $4;
  e[$3, $2] = $4;
}
END{ 
  printf("param bs:\t");
  for(i = 0; i < num_vertices; i++)
    printf(" %s", i);
  print " := ";
  for(i = 0; i < num_vertices; i++) {
    printf("\t\t%s ", i);
    for(j = 0; j < num_vertices; j++)
      printf("%s ", e[i, j]);
    print "";
  }
  print ";";
}'

