#!/usr/bin/env bash

if [ -n "$SIMCOV_BUILD_ENV" ]; then
    source $SIMCOV_BUILD_ENV
fi

upcxx_exec=`which upcxx`

if [ -z "$upcxx_exec" ]; then
    echo "upcxx not found. Please install or set path."
    exit 1
fi

upcxx_exec_canonical=$(readlink -f $upcxx_exec)
if [ "$upcxx_exec_canonical" != "$upcxx_exec" ]; then
    echo "Found symlink for upcxx - using target at $upcxx_exec_canonical"
    export PATH=`dirname $upcxx_exec_canonical`:$PATH
fi

set -e

SECONDS=0

rootdir=`pwd`

INSTALL_PATH=${SIMCOV_INSTALL_PATH:=$rootdir/install}

echo "Installing to $INSTALL_PATH"

rm -rf $INSTALL_PATH/bin/simcov

if [ "$1" == "clean" ]; then
    rm -rf .build
    # if this isn't removed then the the rebuild will not work
    rm -rf $INSTALL_PATH/cmake
    exit 0
else
    mkdir -p $rootdir/.build
    cd $rootdir/.build
    if [ "$1" == "Debug" ] || [ "$1" == "Release" ]; then
        rm -rf *
        rm -rf $INSTALL_PATH/cmake
        #cmake $rootdir -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH -DCMAKE_CXX_COMPILER=/opt/intel/compilers_and_libraries_2020.2.254/linux/mpi/intel64/bin/mpicxx
        cmake $rootdir -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
    fi
    make -j install
fi

echo "Build took $((SECONDS))s"
