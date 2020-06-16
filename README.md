# SimCov #

This is a model for simulating the response of the immune system to infections in the lungs.

## Installing and building

It requires [UPC++](https://bitbucket.org/berkeleylab/upcxx/wiki/Home), C++ and cmake.

This repo contains a submodule, so to install, it's best to run

`git clone --recurse-submodules git@bitbucket.org:shofmeyr/simcov.git`

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

It will create an output directory, which will contain a detailed log file (`simcov.log`). If running in debug mode, a subdirectory
called `per_thread` will appear, with one directory per process that contains debug information produced by that process.

