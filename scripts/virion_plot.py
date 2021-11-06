#!/usr/bin/env python3

import time
import os
import sys
import argparse
import signal
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import numpy

style.use('fivethirtyeight')
#style.use('qpaper')

plt.rc('axes', titlesize=12)
plt.rc('axes', labelsize=10)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('legend', fontsize=8)

argparser = argparse.ArgumentParser()
argparser.add_argument("-f", "--stats-file", required=True, help="Stats file containing simcov results")
argparser.add_argument("-o", "--output", required=True, help="Output file for images")
argparser.add_argument("-r", "--resolution", type=int, dest='resolution', default=1440, help='Resolution: number of time steps per day') 
argparser.add_argument("--log", dest='log_scale', action="store_true", help='Use log scale for epicells and tcells')
options = argparser.parse_args()

fig = plt.figure(figsize=(12, 6))
ax_epicells = fig.add_subplot(1, 2, 1)
ax_virus = fig.add_subplot(1, 2, 2)

def load_data(fname, columns, log_scale=False, scale=1.0, xcol=0, xscale=1):
    graph_data = open(fname,'r').read()
    lines = graph_data.split('\n')
    xs = []
    ys = []
    fields = lines[0][2:].split('\t')
    zero_max = True
    j = 0
    for line in lines[1:]:
        if len(line) == 0:
            continue
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        if xscale == 1:
            xs.append((float(fields[xcol]) + 1) * xscale / options.resolution)
        else:
            xs.append((float(fields[xcol]) * xscale) / options.resolution)
        if log_scale:
            ys.append(1.0 + scale * float(fields[columns[j]]))
        else:
            ys.append(scale * float(fields[columns[j]]))
        if ys[-1] > 0:
            zero_max = False
    if zero_max:
        return
    return xs, ys

def animate(i):
    try:
        maxX = 0
        maxY1 = 0
        maxY2 = 0
        j = 0
        colors = ['blue', 'red', 'orange', 'black']
        # plot dead epicells
        xs, ys = load_data(options.stats_file,
            [4],
            log_scale=options.log_scale)
        ax_epicells.plot(xs, ys, label='dead', lw=2, alpha=1.0, color=colors[j])
        maxX = max(maxX, int(max(xs)))
        maxY1 = max(maxY1, numpy.max(ys))
        # plot virus
        xs, ys = load_data(options.stats_file,
            [8],
            log_scale=options.log_scale)
        maxX = max(maxX, int(max(xs)))
        maxY2 = max(maxY2, numpy.max(ys))
        ax_virus.plot(xs, ys, label='avg virions per cell', lw=2, alpha=1.0, color=colors[j])
    except (ValueError, IndexError) as e:
        print(e)
        return 0
    return maxX, maxY1, maxY2

    
#ani = animation.FuncAnimation(fig, animate, interval=1440)
maxX, maxY1, maxY2 = animate(0)
# plot settings
ax_epicells.legend(loc='upper left')
ax_epicells.set_xlabel('Time (days)')
ax_epicells.set_title('dead epicells')
xticks = ax_epicells.get_xticks()
ax_epicells.set_xticks(range(0, maxX + 1, 1))
if options.log_scale:
    ax_epicells.set_ylim(0.5, 10 * maxY1)
    ax_epicells.set_yscale('log')
ax_virus.legend(loc='upper left')
ax_virus.set_xlabel('Time (days)')
ax_virus.set_title('avg virions per cell')
xticks = ax_virus.get_xticks()
ax_virus.set_xticks(range(0, maxX + 1, 1))
if options.log_scale:
    ax_virus.set_ylim(0.5, 10 * maxY2)
    ax_virus.set_yscale('log')
plt.tight_layout()
# display plot
plt.show()
#plt.savefig(options.output + '.png')
