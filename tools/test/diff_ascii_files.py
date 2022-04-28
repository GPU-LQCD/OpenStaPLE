#!/usr/bin/env python

from sys import argv, stdout
import numpy as np

n1 = np.genfromtxt(argv[1],dtype=float,filling_values= 0) 
n2 = np.genfromtxt(argv[2],dtype=float,filling_values= 0) 

stdout.write("L2 difference: {}".format(np.sum((n1-n2)**2)))
stdout.write(", Max difference: {}\n".format(np.max(np.abs(n1-n2))))



