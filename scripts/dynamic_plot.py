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
#argparser.add_argument("-i", "--indexes", required=True, help="Comma-separated list of column indexes for plotting")
argparser.add_argument("-o", "--output", required=True, help="Output file for images")
argparser.add_argument("-c", "--compare-file", default='', help="File for comparisons")
argparser.add_argument("-e", "--empirical-data", default='', help="File containing empirical data on virus levels")
argparser.add_argument("-r", "--resolution", type=int, dest='resolution', default=1440, help='Resolution: number of time steps per day') 
argparser.add_argument("--virus-scale", type=float, dest='virus_scale', default=1.0, help='Factor to scale comparison virus levels')
argparser.add_argument("--chemo-scale", type=float, dest='chemo_scale', default=1.0, help='Factor to scale comparison chemokine levels')
argparser.add_argument("--log", dest='log_scale', action="store_true", help='Use log scale for epicells and tcells')
options = argparser.parse_args()

#columns = [int(i) for i in options.indexes.split(',')]

fig = plt.figure(figsize=(12, 6))
ax_epicells = fig.add_subplot(2, 2, 1)
ax_tcells = fig.add_subplot(2, 2, 2)
ax_virus = fig.add_subplot(2, 2, 3)
ax_chemo = fig.add_subplot(2, 2, 4)

moddate = None
try:
    moddate = os.stat(options.stats_file)[8]
except FileNotFoundError as err:
    print("File not found: {0}".format(err))
    
unchanged = 0
first = True


def plot_subplot(fname, ax, columns, title, lw=2, alpha=1.0, clear=True, log_scale=False, scale=1.0, xcol=0, xscale=1):
    graph_data = open(fname,'r').read()
    lines = graph_data.split('\n')
    xs = []
    ys = []
    labels = []
    fields = lines[0][2:].split('\t')
    zero_max = True
    for j in range(len(columns)):
        ys.append([])
        labels.append(fields[columns[j]])
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
        for j in range(len(columns)):
            if log_scale:
                ys[j].append(1.0 + scale * float(fields[columns[j]]))
            else:
                ys[j].append(scale * float(fields[columns[j]]))
            if ys[j][-1] > 0:
                zero_max = False
    if clear:
        ax.clear()
    if zero_max:
        return
    colors = ['blue', 'red', 'orange', 'black']
    for j in range(len(columns)):
        #print(title, labels[j], 'max ys', max(ys[j]))
        #if labels[j].startswith('Linear'):
        #    print(xs)
        #    print(ys)
        if labels[j].startswith('Linear'):
            label = 'RNA copies/swab'
        else:
            label = labels[j][0:7]
        if lw > 0:
            ax.plot(xs, ys[j], label=label, lw=lw, alpha=alpha, color=colors[j])
        else:
            ax.plot(xs, ys[j], label=label, lw=0, marker='x', color='red')
    ax.legend(loc='upper left')
    ax.set_xlabel('Time (days)')
    ax.set_title(title)
    xticks = ax.get_xticks()
    if xticks[1] - xticks[0] > 1 and len(xs) > 0:
        ax.set_xticks(range(0, int(max(xs)) + 1, 1))
    if log_scale:
        #if numpy.min(ys) < 0.0001:
        ax.set_ylim(0.5, 10 * numpy.max(ys))
        ax.set_yscale('log')
    plt.tight_layout()
    return xs, ys[j]


def compute_rmsle(emp_xs, emp_ys, sim_xs, sim_ys):
    diffs = []
    for i in range(len(emp_xs)):
        for j in range(len(sim_xs)):
            if int(sim_xs[j]) == int(emp_xs[i]):
                sim = numpy.log(sim_ys[j])
                emp = numpy.log(emp_ys[i])
                diffs.append(numpy.square(sim - emp))
                print(emp_xs[i], '%.2f %.2f %.2f' % (sim, emp, diffs[-1]))
                break

    #print(diffs)
    return numpy.sqrt(numpy.mean(diffs))


def animate(i):
    global moddate
    global unchanged
    global first
    try:
        new_moddate = os.stat(options.stats_file)[8]
    except FileNotFoundError:
        return
    
    try:    
        if new_moddate != moddate or first:
            moddate = new_moddate
            first = False

            plot_subplot(options.stats_file, ax_epicells, [1, 2, 3, 4], 'epicells', log_scale=options.log_scale)
            if options.compare_file != '':
                plot_subplot(options.compare_file, ax_epicells, [1, 2, 3, 4], 'epicells', lw=4, alpha=0.3, clear=False,
                             log_scale=options.log_scale)

            plot_subplot(options.stats_file, ax_tcells, [6, 5], 'tcells', log_scale=options.log_scale)
            if options.compare_file != '':
                plot_subplot(options.compare_file, ax_tcells, [6, 5], 'tcells', lw=4, alpha=0.3, clear=False,
                             log_scale=options.log_scale)

            sim_xs, sim_ys = plot_subplot(options.stats_file, ax_virus, [8], 'avg virions per cell', log_scale=options.log_scale,
                                          scale=options.virus_scale)
            if options.compare_file != '':
                plot_subplot(options.compare_file, ax_virus, [8], 'avg virions per cell', lw=4, alpha=0.3, clear=False,
                             log_scale=options.log_scale, scale=options.virus_scale)
            if options.empirical_data != '':
                emp_xs, emp_ys = plot_subplot(options.empirical_data, ax_virus, [3], 'avg virions per cell', lw=0, alpha=0.3,
                                              clear=False, log_scale=options.log_scale, scale=1.0, xcol=1,
                                              xscale=24*60)
                print('rmsle %.3f' % compute_rmsle(emp_xs, emp_ys, sim_xs, sim_ys))

            plot_subplot(options.stats_file, ax_chemo, [7], 'avg chemokines per cell')
            if options.compare_file != '':
                plot_subplot(options.compare_file, ax_chemo, [7], 'avg chemokines per cell', lw=4, alpha=0.3, clear=False,
                             log_scale=False, scale=options.chemo_scale)
        else:
            unchanged += 1
            if unchanged > 4:
                moddate = new_moddate
                unchanged = 0
                plt.savefig(options.output + '.pdf')
                plt.savefig(options.output + '.png')
    except (ValueError, IndexError) as e:
        print(e)
        return

    
#ani = animation.FuncAnimation(fig, animate, interval=1000)
animate(0)
plt.show()
