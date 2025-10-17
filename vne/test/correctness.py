#!/usr/bin/env python3

import sys

from vallib import *

def test(subfile, virfile, time_s, idf):
  print('%s %s ... ' % (subfile, virfile), end='')
  cplex = runvnep(1, subfile, virfile, time_s)
  bp = runvnep(11, subfile, virfile, time_s)

  print(cplex, bp)

  if (int(bp['@bpTimeLimit']) == 0 and int(bp['@bpInfeasible']) == 1):
    assert int(bp['@cpUnknown']) == 1 or  int(bp['@bpInfeasible']) == 1
  if (int(cplex['@cpOptimal']) == 1 and int(bp['@bpOptimal']) == 1):
    assert float(cplex['@cpCost']) == float(bp['@bpCost'])

  print('OK')
  print_output(cplex, idf)
  print_output(bp, idf)

assert len(sys.argv) == 3, "wrong nb of args: % folder time_s"
folder = sys.argv[1]
time_s = int(sys.argv[2])

for subfile in glob.glob(folder + "/sub*.gtt"):
  for virfile in glob.glob(folder + "/vir*.gtt"):
    idf = getnodes(subfile) + "_" + getnodes(virfile)
    test(subfile, virfile, time_s, idf)
