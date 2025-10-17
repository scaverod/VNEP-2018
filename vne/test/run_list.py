#!/usr/bin/env python3

import sys
import os.path
import argparse
from collections import defaultdict

from vallib import *

parser = argparse.ArgumentParser(description="test VNEP implementations")
parser.add_argument('files', nargs='+', help='files')
parser.add_argument('-v', nargs='+', help='list of versions to run')
parser.add_argument('-t', '--timeout', type=int, default=900, help='time limit')
args = parser.parse_args()

def test(subfile, virfile, time_s, idf, versions):
  print('%s %s ... ' % (subfile, virfile), end=' ')

  results = { v:runvnep(int(v), subfile, virfile, time_s, v) for v in versions}
  optimal = { v:int(results[v]['@cost' + str(v)]) for v in versions if int(results[v]['@optimal' + v]) == 1}
  infeasible = { v:int(results[v]['@infeasible' + str(v)]) for v in versions}
    
  for v,cost in optimal.items():
    for u,cost2 in optimal.items():
      assert cost == cost2, ("wrong optimal value %s => %d and %s => %d" % (v, cost, u, cost2))
      assert not (infeasible[u] == 1 and cost > 0) and not (infeasible[v] == 1 and cost2 > 0), "infeasible"

  print('OK')

  for v in versions:
    print_output(results[v], idf)

for filename in args.files:
  for line in open(filename, 'r'):
    subfile, virfile = line.rstrip().split(' ')
    print(subfile, virfile)
    snodes = countnodes(subfile)
    vnodes = countnodes(virfile)
    idf = "%d %d %s" % (snodes,vnodes,filename)
    test(subfile, virfile, args.timeout, idf, args.v)
