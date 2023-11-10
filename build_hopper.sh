#!/usr/bin/env bash

export UPCXX_THREADMODE=seq
export UPCXX_CODEMODE=opt
export UPCXX_NETWORK=ibv
module purge
module load gcc/12.1.0-crtl
module load cmake/3.11.4-qkyj
module load openmpi/4.1.3-j6zb
module load upcxx/2020.10.0-6eh2

if [ -n "$SIMFORAGER_BUILD_ENV" ]; then
    source $SIMFORAGER_BUILD_ENV
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

INSTALL_PATH=${SIMFORAGER_INSTALL_PATH:=$rootdir/install}

echo "Installing to $INSTALL_PATH"

rm -rf $INSTALL_PATH/bin/simforager

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
    	cmake $rootdir -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH -DCMAKE_CXX_COMPILER=/opt/spack/opt/spack/linux-rocky8-cascadelake/gcc-12.1.0/openmpi-4.1.3-j6zbgs4rx7w7mb4imwl6fqk2wxvglehb/bin/mpicxx
    fi
    make -j install
fi

echo "Build took $((SECONDS))s"
