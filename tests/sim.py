"""
create test data
"""

import copy
import sys
import os
import random
import numpy as np


def read_wig_sim(wig_file, start=0, end=10**10, locus_size=50):
    loci = []
    signal = []
    with open(wig_file, mode='r', encoding='utf-8') as wig:
        for line in wig:
            line = line.strip()
            try:
                line = line.split('\t')
                line[0] = int(line[0])
                line[1] = float(line[1])
                if line[1] < 0:
                    line[1] = 0
                if line[0] < start:
                    continue
                if line[0] > end:
                    break
                if len(loci) > 0 and line[0] - loci[-1] > locus_size:
                    loci.append(loci[-1] + locus_size)
                    signal.append(0)
                    continue

                loci.append(line[0])
                signal.append(line[1])

            except ValueError:
                continue
            except IndexError:
                continue
    return loci, signal

loci_,orig_ = read_wig_sim(sys.argv[1])

for i in range(int(sys.argv[3])):
    loci = copy.copy(loci_)
    orig = copy.copy(orig_)
    for j in range(0,len(orig)):
        if orig[j] >= 1:
            if random.random() < .10:
                orig[j] = max(int(0),int(orig[j] + np.random.normal(0,float(sys.argv[2]),size=1)))
    fname = sys.argv[1] + f'.c{str(i)}.wig'
    f_ = open(fname,'w',encoding='utf-8')
    for x,y in zip(loci,orig):
        f_.write(f'{int(x)}\t{int(y)}\n')
    f_.close()

