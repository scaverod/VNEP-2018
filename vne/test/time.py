#!/usr/bin/env python3

import sys
import os.path
import argparse

from vallib import *

parser = argparse.ArgumentParser()
parser.add_argument('folders', nargs='+')
parser.add_argument('-t', '--timeout', type=int, default=900)
parser.add_argument('-s', default='_1')
parser.add_argument('-m', type=int, default=300)
args = parser.parse_args()

def test(subfile, virfile, time_s, idf):
  print('%s %s ... ' % (subfile, virfile), end='')
  cplex = runvnep(1, subfile, virfile, time_s)
  dfs = runvnep(12, subfile, virfile, time_s)
  bfs = runvnep(13, subfile, virfile, time_s)

  print(cplex, dfs, bfs)

  #if (int(bp['@bpTimeLimit']) == 0 and int(bp['@bpInfeasible']) == 1):
  #  assert int(bp['@cpUnknown']) == 1 or  int(bp['@bpInfeasible']) == 1
  if (int(cplex['@cpOptimal']) == 1 and int(dfs['@0bpOptimal']) == 1):
    assert float(cplex['@cpCost']) == float(dfs['@0bpCost'])

  if (int(cplex['@cpOptimal']) == 1 and int(bfs['@0bpOptimal']) == 1):
    assert float(cplex['@cpCost']) == float(bfs['@0bpCost'])

  print('OK')
  print_output(cplex, idf)
  print_output(dfs, idf)
  print_output(bfs, idf)

for folder in args.folders:
  urlfolder = folder
  if not os.path.exists(folder):
    urlfolder = '../data/%s' % folder
  for subfile in glob.glob(urlfolder + "/sub*" + args.s + ".gtt"):
    for virfile in glob.glob(urlfolder + "/vir*" + args.s + ".gtt"):
      snodes = countnodes(subfile)
      vnodes = countnodes(virfile)
      if snodes < args.m:
        idf = "%d,%d,%s" % (snodes,vnodes,folder)
        test(subfile, virfile, args.timeout, idf)
