#!/bin/bash

#set CMSSW environment
. /cvmfs/cms.cern.ch/cmsset_default.sh

RUN_KERNEL=$(uname -r | cut -d '-' -f1)

if [ "$RUN_KERNEL" == "3.10.0" ]; then
  export SCRAM_ARCH=slc7_amd64_gcc700
  cd /net/cms11/data/pico/cc7/CMSSW_10_2_13/src
elif [ "$RUN_KERNEL" == "2.6.32" ]; then
  cd /net/cms11/data/pico/cc7/CMSSW_10_2_13/src
fi

. /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -

#set up paths
export SMALL_PHYS_UTILS_DIR=$(dirname $(readlink -e "$BASH_SOURCE"))
export PATH=$SMALL_PHYS_UTILS_DIR/bin::$PATH
export PATH=$SMALL_PHYS_UTILS_DIR/scripts::$PATH
export PYTHONPATH=$SMALL_PHYS_UTILS_DIR/lib:$PYTHONPATH
export LD_LIBRARY_PATH=$SMALL_PHYS_UTILS_DIR/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$SMALL_PHYS_UTILS_DIR/inc:$ROOT_INCLUDE_PATH

#set up newest root
source /data1/jbkim/Linux/el7_v1/root-6.24.02/bin/thisroot.sh
export LD_LIBRARY_PATH=/data1/jbkim/Linux/el7_v1/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/data1/jbkim/Linux/el7_v1/lib64:$LD_LIBRARY_PATH

#workaround for strange Python HDF5 bug
export HDF5_USE_FILE_LOCKING='FALSE'
