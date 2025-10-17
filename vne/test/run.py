#!/usr/bin/env python3

import sys
import os.path
import argparse
from collections import defaultdict

from vallib import *

parser = argparse.ArgumentParser(description="test VNEP implementations")
parser.add_argument('folders', nargs='+', help='folders containing the virtual networks')
parser.add_argument('-S', '--sub', help='folder containing the physical substrate networks')
parser.add_argument('-v', nargs='+', help='list of versions to run')
parser.add_argument('-t', '--timeout', type=int, default=900, help='time limit')
parser.add_argument('-s', default='', help='suffix for matching files in the network folders')
parser.add_argument('-m', type=int, default=300, help='maximum number of nodes in the substrate')
parser.add_argument('-i', default='', help='start with this sub file (for interrupted files)')
parser.add_argument('-j', default='', help='start with this vir file (for interrupted files)')
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

for folder in args.folders:
  virfolder = folder
  if not os.path.exists(folder):
    virfolder = '../data/%s' % folder
  
  subfolder = ''
  if not args.sub:
    subfolder = virfolder
  else:
    subfolder = args.sub
    if not os.path.exists(subfolder):
      subfolder = '../data/%s' % subfolder

  for subfile in glob.glob(subfolder + "/sub*" + args.s + ".gtt"):
    # start from file i
    if args.i != '':
      if subfile.find(args.i) >= 0:
        args.i = ''
      else:
        continue

    for virfile in glob.glob(virfolder + "/vir*" + args.s + ".gtt"):
      # start from file j
      if args.j != '':
        if virfile.find(args.j) >= 0:
          args.j = ''
        else:
          continue

      snodes = countnodes(subfile)
      vnodes = countnodes(virfile)
      if snodes < args.m:
        idf = "%d %d %s" % (snodes,vnodes,folder)
        test(subfile, virfile, args.timeout, idf, args.v)
