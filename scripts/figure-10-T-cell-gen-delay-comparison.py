#!/usr/bin/env python3

# HOW TO USE:
# - I have hardcoded the sames of the output directories, and I assume that they are in the same directory as this script
# current directory: figure-10-T-cell-gen-delay-comparison.py <first_output_dir> <second_output_dir> ...
# inside the output dirs, we have the simcov.log, simcov.stats, etc. files
# you can change the directory names and plotting styles on lines 48-50

import matplotlib.pyplot as plt
from matplotlib import style
import argparse
import numpy
import math

style.use('fivethirtyeight')

argparser = argparse.ArgumentParser()
#argparser.add_argument("-f", "--stats-file", required=True)
argparser.add_argument("-o", "--output", required=True)
argparser.add_argument("--log", dest='log_scale', action="store_true")
options = argparser.parse_args()

fig = plt.figure(figsize=(18,7))
font_size = 20

plt.rc('axes', titlesize=font_size)
plt.rc('axes', labelsize=font_size)
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)
plt.rc('legend', fontsize=font_size)

def get_data(fname, xcol, ycol, xscale, yscale):
    graph_data = open(fname,'r').read()
    lines = graph_data.split('\n')
    xs = []
    ys = []
    for line in lines[1:]:
        if len(line) == 0:
            continue
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        xs.append((float(fields[xcol]) * xscale) / 1440)
        ys.append(float(fields[ycol]) * yscale)
    return [xs, ys]

def plot_data(fig, sim_fname, line_style, line_width, legend_labels, color):
    sim_virions = get_data(sim_fname, 0, 8, 1, 1)
    max_sim_virions = int(max(sim_virions[1]))

    fig.plot(sim_virions[0], sim_virions[1], line_style, color=color, lw=line_width)

    fig.legend(loc='upper left', labels=legend_labels, frameon=True, edgecolor='black')
    fig.set_xlabel('Days Post Infection')
    fig.set_ylabel('Number of Virions')
    
    if options.log_scale:
        fig.set_yscale('log')
        
    fig.set_ylim(bottom=1, top=5e7)
    
c = 0
    
runs = ['default-100-3', 'default-50-3', 'default-100-7', 'default-50-7', 'default-100-0', 'default-50-0']
line_styles = ['-', ':', '-', ':', '-', ':']
line_width = [10, 10, 10, 10, 10, 10]
colors = ['red', 'red', 'green', 'green', 'blue', 'blue']
legend_labels = ['100K at day 3', '50K at day 3', '100K at day 7', '50K at day 7', '100K at day 0', '50K at day 0']
ax = fig.add_subplot(1,1,1)
for run in runs:
    f = open(run + '/simcov.log')
    tcell_rate = 0
    tcell_delay = 0
    num_infections = 0
    for line in f.readlines():
        if 'tcell-generation' in line:
            tcell_rate = float(line.split()[-1]) / 100000
        if 'tcell-initial' in line:
            tcell_delay = int(int(line.split()[-1]) / 1440)
        if 'infection-coords' in line:
            num_infections = line.split(':')[-1]
    f.close()
    f = open(run + '/simcov.stats')
    line = f.readlines()[-1]
    perc_infected = line.split()[10]
    plot_data(ax, run + '/simcov.stats', line_styles[c], line_width[c], legend_labels, colors[c])
    f.close()
    c += 1

plt.tight_layout()
plt.savefig(options.output + '.pdf', dpi=500)
plt.savefig(options.output + '.png', dpi=500)
plt.show()


