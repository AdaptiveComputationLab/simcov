#!/usr/bin/env bash

if [ -n "$LUNGMODEL_BUILD_ENV" ]; then
    source $LUNGMODEL_BUILD_ENV
fi

set -e

SECONDS=0

rootdir=`pwd`

INSTALL_PATH=${LUNGMODEL_INSTALL_PATH:=$rootdir}

if [ "$1" == "clean" ]; then
    rm -rf lungmodel $rootdir/.build $INSTALL_PATH/cmake lung_model_data
else
    mkdir -p $rootdir/.build
    cd $rootdir/.build
    rm -rf *
    cmake $rootdir -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
    make -j install
    echo "Build took $((SECONDS))s"
fi
