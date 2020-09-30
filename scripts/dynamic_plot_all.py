#!/usr/bin/env python3

import time
import os
import sys
import argparse
import signal
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style


style.use('fivethirtyeight')
#style.use('qpaper')

plt.rc('axes', titlesize=12)
plt.rc('axes', labelsize=10)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('legend', fontsize=10)


argparser = argparse.ArgumentParser()
argparser.add_argument("-f", "--file", required=True, help="Stats file containing simcov results")
#argparser.add_argument("-i", "--indexes", required=True, help="Comma-separated list of column indexes for plotting")
argparser.add_argument("-o", "--output", required=True, help="Output file for images")
argparser.add_argument("-c", "--compare-file", help="File for comparisons")
options = argparser.parse_args()

#columns = [int(i) for i in options.indexes.split(',')]

fig = plt.figure(figsize=(12, 6))
ax_epicells = fig.add_subplot(1, 3, 1)
ax_tcells = fig.add_subplot(1, 3, 2)
ax_virus = fig.add_subplot(1, 3, 3)

moddate = os.stat('cycells-test/simcov.stats')[8]
unchanged = 0
first = True

def plot_subplot(fname, ax, columns, title, clear=True, log_scale=False):
    graph_data = open(fname,'r').read()
    lines = graph_data.split('\n')
    xs = []
    ys = []
    labels = []
    fields = lines[0][2:].split('\t')
    for j in range(len(columns)):
        ys.append([])
        labels.append(fields[columns[j]])
    for line in lines[1:]:
        if len(line) == 0:
            continue
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        xs.append(float(fields[0]) / (24 * 60))
        for j in range(len(columns)):
            ys[j].append(float(fields[columns[j]]))
    if clear:
        ax.clear()
    for j in range(len(columns)):
        #print(title, labels[j], 'max ys', max(ys[j]))
        ax.plot(xs, ys[j], label=labels[j], lw=2)
    ax.legend()
    ax.set_xlabel('Time (days)')
    ax.set_title(title)
    if log_scale:
        ax.set_yscale('log')
    plt.tight_layout()

def animate(i):
    global moddate
    global unchanged
    global first
    new_moddate = os.stat('cycells-test/simcov.stats')[8]
    if new_moddate != moddate or first:
        moddate = new_moddate
        first = False
        plot_subplot('cycells-test/simcov.stats', ax_epicells, [1, 2, 3], 'epicells')
        if options.compare_file != '':
            plot_subplot(options.compare_file, ax_epicells, [2, 3, 5], 'epicells', clear=False)
        plot_subplot('cycells-test/simcov.stats', ax_tcells, [5, 6], 'tcells')
        if options.compare_file != '':
            plot_subplot(options.compare_file, ax_tcells, [6, 7, 8], 'tcells', clear=False)
        plot_subplot('cycells-test/simcov.stats', ax_virus, [9], 'virus')#, log_scale=True)
        if options.compare_file != '':
            plot_subplot(options.compare_file, ax_virus, [10], 'virus', clear=False)
    else:
        unchanged += 1
        if unchanged > 4:
            moddate = new_moddate
            unchanged = 0
            plt.savefig(options.output + '.pdf')
            plt.savefig(options.output + '.png')


    
ani = animation.FuncAnimation(fig, animate, interval=1000)
plt.show()

