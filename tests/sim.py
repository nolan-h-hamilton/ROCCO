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

locus_size=50
try:
    locus_size = loci_[1]-loci_[0]
except:
    print(f'sim.py: using default locus size {locus_size}')
    pass

for i in range(int(sys.argv[3])):
    loci = copy.copy(loci_)
    orig = copy.copy(orig_)
    for j in range(0,len(orig)):
        if orig[j] > 0:
            if random.random() < .10:
                orig[j] = max(int(0),int(orig[j] + np.random.normal(0,float(sys.argv[2]),size=1)))
    chrom = os.path.basename(sys.argv[1]).split('_')[0]
    odir = os.path.dirname(sys.argv[1])
    fname = f'{odir}/{chrom}_samp{str(i+1)}.wig'
    f_ = open(fname,'w',encoding='utf-8')
    fixedStep = False
    try:
        if sys.argv[4] == 'fs':
            fixedStep = True
    except:
        pass
    if fixedStep:
        f_.write(f'fixedStep chrom={chrom} start={str(loci[0])} step={locus_size}\n')
        for x,y in zip(loci,orig):
            f_.write(f'{int(y)}\n')
        f_.close()
    else:
        for x,y in zip(loci,orig):
            f_.write(f'{int(x)}\t{int(y)}\n')
        f_.close()

