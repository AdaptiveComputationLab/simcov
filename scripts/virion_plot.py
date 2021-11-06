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
argparser.add_argument("-o", "--output", required=True, help="Output file for images")
argparser.add_argument("-r", "--resolution", type=int, dest='resolution', default=1440, help='Resolution: number of time steps per day') 
argparser.add_argument("--log", dest='log_scale', action="store_true", help='Use log scale for epicells and tcells')
options = argparser.parse_args()

fig = plt.figure(figsize=(12, 6))
ax_virus = fig.add_subplot(1, 2, 1)
ax_virus_ds = fig.add_subplot(1, 2, 2)

def load_data(fname, columns, log_scale=False, scale=1.0, xcol=0, xscale=1):
    graph_data = open(fname,'r').read()
    lines = graph_data.split('\n')
    xs = []
    ys = []
    fields = lines[0][2:].split('\t')
    zero_max = True
    j = 0
    count = 0
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
        count = count + 1
        if count > 10080:
            break
    if zero_max:
        return
    return xs, ys

def animate(i):
    try:
        maxY = 0
        j = 0
        colors = ['orange', 'purple', 'red']
        labels = ['2D Simulation', '3D Branching Topology', '3D Simuation']
        for fl in ['2d', 'topo', '3d']:
            # plot virus
            fname = 'simcov.' + fl + '.stats'
            xs, ys = load_data(fname, [8], log_scale=options.log_scale)
            ax_virus.plot(xs, ys, label=labels[j], lw=2, alpha=1.0, color=colors[j])
            maxY = max(maxY, numpy.max(ys))
            # plot virus downscaled axes
            ax_virus_ds.plot(xs, ys, label=labels[j], lw=2, alpha=1.0, color=colors[j])
            j = j + 1
    except (ValueError, IndexError) as e:
        print(e)
        return
    return maxY

    
#ani = animation.FuncAnimation(fig, animate, interval=1440)
maxY = animate(0)
maxY1 = 5e8 #TODO 5160000000.0
maxY2 = 5e7
# plot settings
days = 7
# configure virus plot
ax_virus.legend(loc='upper left')
ax_virus.set_xlabel('Time (days)')
ax_virus.set_ylabel('Virion Count')
xticks = ax_virus.get_xticks()
ax_virus.set_xticks(range(0, days + 1, 1))
if options.log_scale:
    ax_virus.set_ylim(0.5, 10 * maxY1)
    ax_virus.set_yscale('log')
else:
    ax_virus.set_ylim(0.5, maxY1)
# configure virus plot downscaled
ax_virus_ds.legend(loc='upper right')
ax_virus_ds.set_xlabel('Time (days)')
ax_virus_ds.set_ylabel('Virion Count')
xticks = ax_virus_ds.get_xticks()
ax_virus_ds.set_xticks(range(0, days + 1, 1))
if options.log_scale:
    ax_virus_ds.set_ylim(0.5, 10 * maxY2)
    ax_virus_ds.set_yscale('log')
else:
    ax_virus_ds.set_ylim(0.5, maxY2)
# display plot
plt.tight_layout()
plt.show()
#plt.savefig(options.output + '.png')
