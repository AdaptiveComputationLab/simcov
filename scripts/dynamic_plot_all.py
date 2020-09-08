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

argparser = argparse.ArgumentParser()
argparser.add_argument("-f", "--file", required=True, help="Stats file containing simcov results")
#argparser.add_argument("-i", "--indexes", required=True, help="Comma-separated list of column indexes for plotting")
argparser.add_argument("-o", "--output", required=True, help="Output file for images")
options = argparser.parse_args()

#columns = [int(i) for i in options.indexes.split(',')]

fig = plt.figure(figsize=(12, 6))
ax_epicells = fig.add_subplot(1, 3, 1)
ax_tcells = fig.add_subplot(1, 3, 2)
ax_virus = fig.add_subplot(1, 3, 3)

moddate = os.stat('cycells-test/simcov.stats')[8]
unchanged = 0
first = True

def plot_subplot(ax, columns, title):
    graph_data = open('cycells-test/simcov.stats','r').read()
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
    ax.clear()
    for j in range(len(columns)):
        ax.plot(xs, ys[j], label=labels[j], lw=2)
    ax.legend()
    ax.set_xlabel('Time (days)')
    ax.set_title(title)
    plt.tight_layout()

def animate(i):
    global moddate
    global unchanged
    global first
    new_moddate = os.stat('cycells-test/simcov.stats')[8]
    if new_moddate != moddate or first:
        moddate = new_moddate
        first = False
        plot_subplot(ax_epicells, [1, 2, 3], 'epithelial cells')
        plot_subplot(ax_tcells, [5, 6], 'T cells')
        plot_subplot(ax_virus, [9], 'viral load')
    else:
        unchanged += 1
        if unchanged > 4:
            moddate = new_moddate
            unchanged = 0
            plt.savefig(options.output + '.pdf')
            plt.savefig(options.output + '.png')


    
ani = animation.FuncAnimation(fig, animate, interval=1000)
plt.show()

