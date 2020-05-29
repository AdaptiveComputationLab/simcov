#!/usr/bin/env bash

if [ -n "$SIMCOV_BUILD_ENV" ]; then
    source $SIMCOV_BUILD_ENV
fi

set -e

SECONDS=0

rootdir=`pwd`

INSTALL_PATH=${SIMCOV_INSTALL_PATH:=$rootdir/install}

rm -rf $INSTALL_PATH/bin/simcov

if [ "$1" == "clean" ]; then
    rm -rf .build/*
    # if this isn't removed then the the rebuild will not work
    rm -rf $INSTALL_PATH/cmake
    exit 0
else
    mkdir -p $rootdir/.build
    cd $rootdir/.build
    if [ "$1" == "Debug" ] || [ "$1" == "Release" ]; then
        rm -rf .build/*
        rm -rf $INSTALL_PATH/cmake
        cmake $rootdir -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
    fi
    make -j install
fi

if [ -n "$SIMCOV_BUILD_ENV" ]; then
    env_id=`echo $SIMCOV_BUILD_ENV|cut -d '.' -f2|cut -d '-' -f2-`
    cd $INSTALL_PATH/bin
    rm -f simcov.$env_id
    mv simcov simcov.$env_id
    ln -s simcov.$env_id simcov
fi

echo "Build took $((SECONDS))s"
