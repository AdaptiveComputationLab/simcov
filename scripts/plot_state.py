import numpy as np
import vtk
import pyvista as pv
import matplotlib.pyplot as plt
import seaborn as sb
import sys

def plot_virus(filename, outputName):
    mesh = pv.read(filename)
    data = np.reshape(mesh["virus"], (mesh.dimensions[0]-1, mesh.dimensions[1]-1))
    plt.figure()
    ax = sb.heatmap(data, cmap="viridis")
    plt.title("Viral Concentration")
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([])
    plt.xlabel("{} Microns".format((mesh.dimensions[0]-1)*5))
    plt.ylabel("{} Microns".format((mesh.dimensions[0]-1)*5))
    plt.xticks([])
    plt.yticks([])
    plt.savefig(outputName + "_virus.png", dpi=300)

def plot_chemokine(filename_ch, filename_tc, outputName):
    mesh1 = pv.read(filename_ch)
    data1 = np.reshape(mesh1["chemokine"], (mesh1.dimensions[0]-1, mesh1.dimensions[1]-1))
    plt.figure()
    ax = sb.heatmap(data1, cmap="Blues")
    cmap2 = sb.color_palette([(0xD1/0xFF, 0xEC/0xFF, 0x9C/0xFF, 0), (0xF1/0xFF, 0xEB/0xFF, 0xF4/0xFF, 1.0)])
    plt.title("Chemokines and t-cells in tissue")
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([])
    plt.xlabel("{} Microns".format((mesh1.dimensions[0]-1)*5))
    plt.ylabel("{} Microns".format((mesh1.dimensions[0]-1)*5))
    plt.xticks([])
    plt.yticks([])
    plt.savefig(outputName + "_chemokine.png", dpi=300)

def plot_epicells(filename, outputName):
    mesh = pv.read(filename)
    data = np.reshape(mesh["epicell"], (mesh.dimensions[0]-1, mesh.dimensions[1]-1))
    plt.figure()
    n = 5
    colors = ["#b8b8b8", "#ffca4f", "#ff3b05", "#197dff", "#000000"]
    cmap = sb.color_palette(colors)
    ax = sb.heatmap(data, cmap = cmap)
    colorbar = ax.collections[0].colorbar
    r = colorbar.vmax - colorbar.vmin 
    colorbar.set_ticks([colorbar.vmin + r / n * (0.5 + i) for i in range(n)])
    colorbar.set_ticklabels(["Healthy", "Incubating", "Expressing", "Apoptotic", "Dead"])  
    plt.xlabel("{} Microns".format((mesh.dimensions[0]-1)*5))
    plt.ylabel("{} Microns".format((mesh.dimensions[0]-1)*5))
    plt.title("Epithelial Cell Status")
    plt.xticks([])
    plt.yticks([])
    plt.savefig(outputName + "_epicells.png", dpi=300)

sample_directory = sys.argv[1]
timestep = sys.argv[2]
outputName = sys.argv[3]

epifile = "{}/sample_epicell_{}.vtk".format(sample_directory, timestep)
virusfile = "{}/sample_virus_{}.vtk".format(sample_directory, timestep)
chemokinefile = "{}/sample_chemokine_{}.vtk".format(sample_directory, timestep)
tcellfile = "{}/sample_tcelltissue_{}.vtk".format(sample_directory, timestep)


plot_epicells(epifile, outputName)
plot_virus(virusfile, outputName)
plot_chemokine(chemokinefile, tcellfile, outputName)