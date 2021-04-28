import numpy as np
import vtk
import pyvista as pv
import matplotlib.pyplot as plt
import seaborn as sb
import sys
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib
import argparse

def plot_virus_single(filename, output):
    fig = plt.figure(figsize=(12,10), dpi=300)
    eps = 10e-7
    mesh = pv.read(filename)
    data = np.reshape(mesh["virus"], (mesh.dimensions[0]-1, mesh.dimensions[1]-1))
    ax = sb.heatmap(data, cmap="viridis", cbar=True, vmin=0.0, vmax=255.0)
    # ax.set_xlim([800,2800])
    # ax.set_ylim([400,2400])
    ax.set_xticks([])
    ax.set_yticks([])
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([0, 255*np.log2(100000)/np.log2(125000)])
    colorbar.set_ticklabels(["0",r"$10^{5}$"])
    colorbar.ax.set_yticklabels(["0",r"$10^{5}$"])
    colorbar.set_label("Virion Count",labelpad=-30)
    plt.savefig(output)
    
def plot_chemokine_single(filename_ch, filename_tc, output):
    plt.figure(figsize=(10,10), dpi=300)
    mesh1 = pv.read(filename_ch)
    data1 = np.reshape(mesh1["chemokine"], (mesh1.dimensions[0]-1, mesh1.dimensions[1]-1))
    mesh2 = pv.read(filename_tc)
    data2 = np.reshape(mesh2["t-cell-tissue"], (mesh2.dimensions[0]-1, mesh2.dimensions[1]-1))
    x = []
    y = []
    for i in range(mesh2.dimensions[0]-1):
        for j in range(mesh2.dimensions[1]-1):
            if(data2[i,j] > 0):
                x.append(j)
                y.append(i)
    ax = sb.heatmap(data1, cmap="Greys", vmax=254, vmin=-254, cbar=False)
    # ax.set_xlim([800,2800])
    # ax.set_ylim([400,2400])
    plt.scatter(x,y,color="#009900",s=15,alpha=0.8)
    cmap2 = sb.color_palette([(0xD1/0xFF, 0xEC/0xFF, 0x9C/0xFF, 0), (0xF1/0xFF, 0xEB/0xFF, 0xF4/0xFF, 1.0)])
    ax.set_xticks([])
    ax.set_yticks([])
    plt.savefig(output)
    
def plot_epicells_single(filename, output):
    plt.figure(figsize=(10,10), dpi=300)
    mesh = pv.read(filename)
    data = np.reshape(mesh["epicell"], (mesh.dimensions[0]-1, mesh.dimensions[1]-1))
    n = 5
    colors = ["#b8b8b8", "#ffca4f", "#ff3b05", "#197dff", "#000000"]
    cmap = sb.color_palette(colors)
    ax = sb.heatmap(data, cmap = cmap, cbar=False)
    # colorbar = ax.collections[0].colorbar
    # r = colorbar.vmax - colorbar.vmin 
    # colorbar.set_ticks([colorbar.vmin + r / n * (0.5 + i) for i in range(n)])
    # colorbar.set_ticklabels(["Healthy", "Incubating", "Expressing", "Apoptotic", "Dead"])  
    # ax.set_xlim([800,2800])
    # ax.set_ylim([400,2400])
    ax.set_xticks([])
    ax.set_yticks([])
    plt.savefig(output)

parser = argparse.ArgumentParser()
parser.add_argument("--plot",
    metavar="P",
    type=str,
    help="What plot to generate? [epicell, virus, chemokine]",
    choices=["epicell", "virus", "chemokine"],
    required=True)
parser.add_argument("--input",
    metavar="I",
    type=str,
    help="Input VTK sample File String",
    required=True)
parser.add_argument("--output",
    metavar="I",
    type=str,
    help="Output image file",
    required=True)
parser.add_argument("--input_tc",
    metavar="IT",
    type=str,
    help="Input T-Cell VTK if plotting chemokine plot")

args = parser.parse_args()
plottype = args.plot
inputfile = args.input
outputfile = args.output
input_tc = args.input_tc

if plottype == "epicell":
    plot_epicells_single(inputfile, outputfile)
elif plottype == "virus":
    plot_virus_single(inputfile, outputfile)
elif plottype == "chemokine":
    plot_chemokine_single(inputfile, input_tc, outputfile)