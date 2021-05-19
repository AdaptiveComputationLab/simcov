# Lung Model #

## Installing and building

Build it with

`./build.sh`

The executable will be installed in this working directory.

## Running

To run, execute

`./lungmodel --dim 300 300 300 --levels 3 --scale 2000 --output lung_model_data`

It will create an 'lung_model_data' output directory whose full filepath is required in the simcov config file. For example:

```
; Directory containing files for lung model
  lung-model =                  /users/projects/simcov/lungmodel/lung_model_data

```
