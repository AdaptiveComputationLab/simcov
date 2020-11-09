# SimCov #

This is a model for simulating the response of the immune system to infections in the lungs.

## Installing and building

It requires [UPC++](https://bitbucket.org/berkeleylab/upcxx/wiki/Home), C++ and cmake.

This repo contains a submodule, so to install, it's best to run

`git clone --recurse-submodules git@github.com:AdaptiveComputationLab/simcov.git`

to fully initialize the submodules.

Once downloaded, build it with

`./build.sh Release`

or

`./build.sh Debug`

The executable will be installed in

`<simcov-repo-directory>/install/bin`

## Running

To run, execute

`upcxx-run -n <number_processes> -N <number_nodes> -- simcov`

To see the parameters available, run with `-h`.

It will create an output directory, which will contain a detailed log file (`simcov.log`). It will also create a file containing
all the configuration parameters (`simcov.config`).

A run can also be executed with a config file as

`upcxx-run -n <number_processes> -N <number_nodes> -- simcov --config <config_file>`

The config file consists of a list of all the command line options as key-value pairs, with semi-colons denoting comments.
For example:

```
; Dimensions: x y z
  dim =                    100,100,100

; Number of timesteps
  timesteps =                   14000

```

Any options not specified in the config file will be set to the defaults (seen when run with `-h`).

When running with a config file, any parameters passed on the command line will override parameters read from the config file.


If running in debug mode, a subdirectory
called `per_thread` will appear, with one directory per process that contains debug information produced by that proc

Once the run has completed, the outputs can be viewed in paraview by opening the samples subdirectory of the output directory. To help with viewing, a python script can generate a paraview state file:

```
pvpython sample-state.py --data <output_dir>/samples
```

Note that the `pvpython` wrapper is needed (not just plain python).

The above command will generate a file called `samples-state.pvsm` and this can be loaded into paraview, e.g.:

```
paraview sample-state.pvsm
```
