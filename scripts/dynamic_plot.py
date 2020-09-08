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
argparser.add_argument("-i", "--indexes", required=True, help="Comma-separated list of column indexes for plotting")
argparser.add_argument("-o", "--output", required=True, help="Output file for images")
options = argparser.parse_args()

columns = [int(i) for i in options.indexes.split(',')]

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

moddate = os.stat('cycells-test/simcov.stats')[8]
unchanged = 0

def animate(i):
    global moddate
    global unchanged
    new_moddate = os.stat('cycells-test/simcov.stats')[8]
    if new_moddate != moddate:
        moddate = new_moddate
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
        ax1.clear()
        for j in range(len(columns)):
            ax1.plot(xs, ys[j], label=labels[j])
        plt.legend()
        plt.xlabel('Time (days)')
        plt.tight_layout()
    else:
        unchanged += 1
        if unchanged > 4:
            moddate = new_moddate
            unchanged = 0
            plt.savefig(options.output + '.pdf')


    
ani = animation.FuncAnimation(fig, animate, interval=1000)
plt.show()

