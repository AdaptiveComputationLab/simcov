## SIMCoV Hopper Tutorial

This tutorial contains instructions for compiling and running the SIMCoV immunology model on the Hopper cluster.

## Download the source code from GitHub

Change directory to your home: 
```
cd ~
```

Clone the simcov Github repository into your home directory:
```
git clone --recurse-submodules https://github.com/AdaptiveComputationLab/simcov.git
```

## Build SIMCoV

Change to the root directory just cloned: 
```
cd simcov
```

Load Hopper modules, set UPCXX variables, and build SIMCoV:
```
./build_hopper.sh Release
```

## Configure SIMCoV

Simulation parameters with descriptions are listed in Table 2 on page 15 of the paper with parameter derivations on pages 16-18. The configuration files used in research paper are in SIMCoV root directory and end with ".config". You can edit them with a text editor. Reporting on modifying parameters can be seen in Figure 5 page 8 of the paper (config files used are in 'configs/Figure-5' folder).  

## Run SIMCoV 

On Hopper use the Simple Linux Utility for Resource Management (SLURM) to run SIMCoV on compute nodes. The Hopper slurm script 'slurm_hopper.sh' is provided that loads 'covid_default.config' configuration file.

To run SIMCoV on a compute node enter:
```
sbatch slurm_hopper.sh
```

Check the status of your simuation in the batch system:
```
squeue -u <your_user_name>
```

View job console outputs (e.g. cat slurm-35489.out):
```
cat slurm-<your_job_id>.out
```

## Examine Results

By default simulation outputs to 'results' folder in the SIMCoV root directory by default. Change directory to the results folder and examine the simulation log:
```
cd results
cat simcov.log
```

The 'simcov.stats' file contains simulation outputs in column format (time step, incubating, expressing, apoptotic, dead, vasculature tcells, tissue tcells, chemokines, virus, chemokine pts, percentage of infected cells).  

The 'samples' folder contains .vtk files that can be visualized in ParaView.  

## Visualize results using ParaView

Enable .vtk outputs for visualiztions by editing the sampling period parameter in the config file (e.g. 1 day == 1440 mins):
```
; Number of timesteps between samples (set to 0 to disable sampling)
sample-period =               1440
```

Download ParaView software from [here](https://www.paraview.org/download/) by selecting the version then your platform (not MPI versions). In the example below are brief instructions for running ParaView on your local computer. ParaView is also installed on Hopper allowing multiple nodes to be used for visualizing large datasets. You can find the full instructions on [CARC QuickBytes Tutorials](https://github.com/UNM-CARC/QuickBytes/blob/master/paraview.md).

## Examples: Modifying parameters, plotting results with Python, and visualizing results with ParaView

Simulation parameters with descriptions are listed in Table 2 on page 15 of the paper with parameter derivations on pages 16-18. The configuration files used in research paper are in SIMCoV root directory and end with ".config". You can edit them with a text editor. Reporting on modifying parameters can be seen in Figure 5 page 8 of the paper (config files used are in 'configs/Figure-5' folder).  

Create three simulations with the following parameter changes:
1. Default
```
dim =                         1500 1500 1
timesteps =                   14400
sample-period =               1440
output =                      results-tutorial-1
```
2. Default and 10x virus clearance
```
virion-clearance =            0.04
output =                      results-tutorial-2
```
3. Default and /100 tcell production
```
tcell-generation-rate = 	1050
output =                      results-tutorial-3
```

The 'scripts' directory contains Python scripts used to generate plots in the research paper. For more information on using Python on Hopper see full instructions on [CARC QuickBytes Tutorials](https://github.com/UNM-CARC/QuickBytes/blob/master/anaconda_intro.md). Use Python to generate linear plots of results data contained in 'simcov.stats' file and save to 'Figure.png' file:
```
module purge
module load miniconda3-4.10.3-gcc-11.2.0-tl4rbd6
conda create --name SimCovTutorial python=3 matplotlib numpy -y
conda init bash (OPTIONAL IF NOT ALREADY DONE. MUST RESTART SHELL!!!)
source activate SimCovTutorial
python scripts/dynamic_plot.py -f results/simcov.stats -o Figure --log
conda deactivate
```

The 'results/samples' directory contains simulation results for virus, chemokine (inflammatory signal), T cell, and epithelial cells stored in .vtk file format. Use ParaView to generate figures of results data contained in *.vtk files:
1. File->Open
2. Select 'Pipeline Browser' tab to display loaded data files (e.g. sample_virus_*)
3. Select a sample file to display (the eye icon appears next to file name)
4. Select 'Properties' Tab
5. Click 'Apply' button
6. Select the Representation: Outline -> Volume  

On the top menu bar click 'Next Frame' button to display each timestep and 'Play' button to animate frames. Using the top file menu option create a screenshot using 'File->Save ScreenShot...' and a video using 'File->Save Animation...'.

The 'scripts' directory also contains a ParaView Python (pvpyton) script 'generate_paraview_state.py' that generates a ParaView state file from simulation results that can be loaded using the top file menu option 'File->Load State...'. You can find full instructions on [SIMCoV Readme.md](https://github.com/AdaptiveComputationLab/simcov).

## Running SIMCoV with Branching Airways

A description of the Branching Airways is on pages 10-12 of the paper and results visualized in Figures 7 and 8. Do the following to configure simcov to run with Branching Airways:

1. Change to simcov root directory and checkout the latest branch
```
cd simcov
git checkout branching_airways_3min_model
```

2. Edit the covid_ba.config file to point to your simcov directory (__*Note:__ Current BA build is configured for a 3 minute timestep where simcov baseline uses 1 minute timestep)
```
; *Directory containing files for lung model
lung-model =                  /users/<Your-User-Name>/simcov
```

3. Build and submit batch job
```
./build_hopper.sh Release
sbatch slurm_hopper.sh
```

4. Check the status of your simuation in the batch system:
```
squeue -u <your_user_name>
```

5. View job console outputs (e.g. cat slurm-35489.out):
```
cat slurm-<your_job_id>.out
```

6. Launch ParaView to visualize results (see [Examples](#examples-modifying-parameters-plotting-results-with-python-and-visualizing-results-with-paraview) Section)
