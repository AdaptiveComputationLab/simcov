#!/usr/bin/env python3

import time
import os
import sys
import argparse
import signal
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import pandas
import numpy


fname = sys.argv[1]

style.use('fivethirtyeight')

plt.rc('axes', titlesize=12)
plt.rc('axes', labelsize=10)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('legend', fontsize=8)


df = pandas.read_csv(fname, '\t')
t = df[df.columns[1]]
df_t = df[t == 99]
x = df_t[df.columns[2]]
y = df_t[df.columns[3]]
zlog10 = numpy.log10(df_t[df.columns[5]])
zlog2 = numpy.log2(df_t[df.columns[5]])

print('Area:', len(zlog2))
    
fig = plt.figure()
ax = fig.add_subplot(111)
#plt.xlim(min(30, 10 * int(min(x) / 10)), max(70, 10 * int(max(x + 10) / 10)))
#plt.ylim(min(30, 10 * int(min(y) / 10)), max(70, 10 * int(max(x + 10) / 10)))
plt.xlim(30, 70)
plt.ylim(30, 70)

if False:
    z = zlog10
    levels = numpy.arange(int(min(z) - 10), int(max(z) + 10), 1)
else:
    z = zlog2
    levels = numpy.arange(int(min(z) - 2), int(max(z) + 2), 1)
    
cs = ax.tricontour(x, y, z, linewidths=0.5, colors='k', levels=levels)
ax.clabel(cs, cs.levels, fmt='%0d', inline=True, fontsize=8)
cf = ax.tricontourf(x, y, z, levels=80)
fig.colorbar(cf, ax=ax, ticks=levels)
ax.set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.savefig(fname + '.png')
plt.show()

