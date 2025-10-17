#!/usr/bin/env python3

from collections import defaultdict
import os
import sys
import glob
import time

def countnodes(gttfilename):
  with open(gttfilename, 'r') as f:
    return int(f.readline().split(' ')[1])

def runvnep(code, subfile, virfile, time_s, fieldsuffix=''):
  inifilename = 'p' + str(time.time()) + '.ini'
  try:
    inifile = open(inifilename, 'w')
    inifile.write('substrate_file = %s \n' % subfile)
    inifile.write('virtual_file = %s \n' % virfile)
    inifile.write('timelimit_in_s = %i \n' % time_s)
    inifile.close()
    output = os.popen('bin/vne.out %d %s' % (code, inifilename), 'r')

    data = defaultdict(lambda: '0')
    for line in output:
      if line[0] == '@':
        s = line.split(':')
        data[s[0] + fieldsuffix] = s[1].rstrip()
    return data
  finally:
    os.system("rm -f " + inifilename)

def print_output(output, idf):
  for key,value in output.items():
    print(key, ":", value, ":", idf)

