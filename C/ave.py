#!/usr/bin/env python3.7
# coding: utf-8

import numpy as np
import sys

filename = sys.argv[1]

f  = open(filename)
line = f.readline()
print(line)

data = []
line = f.readline()
while(line):
    data.append(float(line[:-2]))
    line = f.readline()

data = np.array(data)
print('%.6e +- %.6e'%(np.mean(data),np.std(data)/np.sqrt(len(data))))
