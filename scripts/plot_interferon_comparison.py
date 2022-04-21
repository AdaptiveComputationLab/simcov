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
ax_prod = fig.add_subplot(1, 2, 1)
ax_infect = fig.add_subplot(1, 2, 2)

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
        maxY1 = 0
        j = 0
        colors = ['orange', 'purple', 'red', 'black', 'blue']
        for fl in ['05', '1', '125', '15', '75']:
            # plot virus
            fname = 'simcov.' + fl + '.stats'
            xs, ys = load_data(fname, [8], log_scale=options.log_scale)
            ax_prod.plot(xs, ys, label=fl, lw=2, alpha=1.0, color=colors[j])
            maxY1 = max(maxY1, numpy.max(ys))
            j = j + 1
        maxY2 = 0
        j = 0
        for fl in ['7', '8', '9', '75i', '85']:
            # plot virus
            fname = 'simcov.' + fl + '.stats'
            xs, ys = load_data(fname, [8], log_scale=options.log_scale)
            ax_infect.plot(xs, ys, label=fl, lw=2, alpha=1.0, color=colors[j])
            maxY2 = max(maxY2, numpy.max(ys))
            j = j + 1
    except (ValueError, IndexError) as e:
        print(e)
        return
    return maxY1, maxY2

    
#ani = animation.FuncAnimation(fig, animate, interval=1440)
maxY1, maxY2 = animate(0)
#maxY1 = 5e8 #TODO 5160000000.0
#maxY2 = 5e7
# plot settings
days = 7
# configure virus production plot
ax_prod.legend(loc='upper left')
ax_prod.set_xlabel('Time (days)')
ax_prod.set_ylabel('Virion Count')
ax_prod.set_title('Virus Production')
xticks = ax_prod.get_xticks()
ax_prod.set_xticks(range(0, days + 1, 1))
if options.log_scale:
    ax_prod.set_ylim(0.5, 10 * maxY1)
    ax_prod.set_yscale('log')
else:
    ax_prod.set_ylim(0.5, maxY1)
# configure virus infectivity plot
ax_infect.legend(loc='upper right')
ax_infect.set_xlabel('Time (days)')
ax_infect.set_ylabel('Virion Count')
ax_infect.set_title('Virus Infectivity')
xticks = ax_infect.get_xticks()
ax_infect.set_xticks(range(0, days + 1, 1))
if options.log_scale:
    ax_infect.set_ylim(0.5, 10 * maxY2)
    ax_infect.set_yscale('log')
else:
    ax_infect.set_ylim(0.5, maxY2)
# display plot
plt.tight_layout()
plt.show()
#plt.savefig(options.output + '.png')
