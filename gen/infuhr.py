#!/usr/bin/env python3

import os

root = '../data/infuhr2011_original'
output = '../data/infuhr2011'

def gtt(subfile, outputfile):
  for line in subfile:
    if line.startswith('#'):
      state = line.split(':')[0][2:].rstrip()
    else:
      info = line.rstrip().split(' ')
      if state == 'noVertices':
        v = int(line)
      elif state == 'noArcs':
        print('G', v, int(line), 'infuhr2011', file=outputfile)
      elif state == 'Vertices':
        print('V', int(info[0]), int(info[1]), file=outputfile)
      elif state == 'Arcs':
        print('E', int(info[0]), int(info[1]), int(info[3]), file=outputfile)
    

for dirname, dirnames, files in os.walk(root):
  for filename in files:
    filedesc = filename.split('.')[0]
    if filename.endswith('.substrate'):
      infile = open(dirname + '/' + filename, 'r')
      outfile = open(output + '/sub_' + filedesc + ".gtt" , 'w')
      gtt(infile, outfile)
    elif filename.endswith('.slices'):
      infile = open(dirname + '/' + filename, 'r')
      outfile = open(output + '/vir_' + filedesc + ".gtt" , 'w')
      gtt(infile, outfile)
