import numpy as np
import vtk
import pyvista as pv
import matplotlib.pyplot as plt
import seaborn as sb
import sys
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib


def plot_virus(filename, axes):
    eps = 10e-7
    mesh = pv.read(filename)
    data = np.reshape(mesh["virus"], (mesh.dimensions[0]-1, mesh.dimensions[1]-1))
    ax = sb.heatmap(data, cmap="viridis", ax=axes, cbar=False)
    axes.set_xticks([])
    axes.set_yticks([])

def plot_chemokine(filename_ch, filename_tc, axes):
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
    ax = sb.heatmap(data1, cmap="Greys", vmax=254, vmin=-254, ax=axes, cbar=False)
    axes.scatter(x,y,color="#009900",s=5,alpha=0.15)
    cmap2 = sb.color_palette([(0xD1/0xFF, 0xEC/0xFF, 0x9C/0xFF, 0), (0xF1/0xFF, 0xEB/0xFF, 0xF4/0xFF, 1.0)])
    axes.set_xticks([])
    axes.set_yticks([])

def plot_epicells(filename, axes):
    mesh = pv.read(filename)
    data = np.reshape(mesh["epicell"], (mesh.dimensions[0]-1, mesh.dimensions[1]-1))
    n = 5
    colors = ["#b8b8b8", "#ffca4f", "#ff3b05", "#197dff", "#000000"]
    cmap = sb.color_palette(colors)
    ax = sb.heatmap(data, cmap = cmap, ax=axes, cbar=False)
    # colorbar = ax.collections[0].colorbar
    # r = colorbar.vmax - colorbar.vmin 
    # colorbar.set_ticks([colorbar.vmin + r / n * (0.5 + i) for i in range(n)])
    # colorbar.set_ticklabels(["Healthy", "Incubating", "Expressing", "Apoptotic", "Dead"])  
    axes.set_xticks([])
    axes.set_yticks([])


def plot_3_panel(epifile, virusfile, chemokinefile, tcellfile, outputName):
	fig, axs = plt.subplots(1,3,figsize=(16,5))
	plot_epicells(epifile, axs[0])
	plot_virus(virusfile, axs[1])
	plot_chemokine(chemokinefile, tcellfile, axs[2])
	plt.tight_layout()
	plt.savefig("{}_3panel.png".format(outputName), dpi=300)


matplotlib.rcParams.update({"font.size":22})

sample_directory = sys.argv[1]
timestep = sys.argv[2]
outputName = sys.argv[3]

epifile = "{}/sample_epicell_{}.vtk".format(sample_directory, timestep)
virusfile = "{}/sample_virus_{}.vtk".format(sample_directory, timestep)
chemokinefile = "{}/sample_chemokine_{}.vtk".format(sample_directory, timestep)
tcellfile = "{}/sample_tcelltissue_{}.vtk".format(sample_directory, timestep)

plot_3_panel(epifile, virusfile, chemokinefile, tcellfile, outputName)


# plot_epicells(epifile, outputName)
# plot_virus(virusfile, outputName)
# plot_chemokine(chemokinefile, tcellfile, outputName)
