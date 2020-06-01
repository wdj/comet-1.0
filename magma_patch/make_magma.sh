#!/bin/bash
#==============================================================================
#
# Script to build a modified version of the MAGMA library.
#
# Usage: type the following from the respective Magma root directory, e.g.,
# from magma_minproduct:
#
# ../make_magma.sh
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function do_make
{
  local IS_CRAY_XK7 # OLCF Titan or Chester
  [[ -n "${CRAYOS_VERSION:-}" ]] && IS_CRAY_XK7="YES" || IS_CRAY_XK7="NO"
  local IS_IBM_AC922 # OLCF Summit or Peak
  [[ -n "${LSF_BINDIR:-}" ]] && IS_IBM_AC922="YES" || IS_IBM_AC922="NO"
  local IS_DGX2
  [[ "$(uname -n)" = "dgx2-b" ]] && IS_DGX2="YES" || IS_DGX2="NO"
  local IS_GPUSYS2
  [[ "$(uname -n)" = "gpusys2" ]] && IS_GPUSYS2="YES" || IS_GPUSYS2="NO"
  local IS_EXPERIMENTAL
  [[ "${COMET_BUILD_EXPERIMENTAL:-}" = YES ]] && IS_EXPERIMENTAL="YES" || \
                                                 IS_EXPERIMENTAL="NO"

  # WARNING: these module loads MUST match those in scripts/cmake.sh
  if [ $IS_EXPERIMENTAL = YES ] ; then
    true # skip for now
  elif [ $IS_CRAY_XK7 = YES ] ; then
    if [ "$PE_ENV" = "PGI" ] ; then
      module unload PrgEnv-pgi
    fi
    module load PrgEnv-gnu
    module load cudatoolkit
    module load acml
    module list 2>&1
    cp ../make.inc.titan make.inc
    export GPU_TARGET=sm35 NV_SM=" -gencode arch=compute_35,code=sm_35" NV_COMP="-gencode arch=compute_35,code=compute_35" MIN_ARCH=350
    MAKE_ARGS=""
  elif [ $IS_IBM_AC922 = YES ] ; then
    module load gcc/6.4.0
    module load cuda
    module list 2>&1
    export CUDA_DIR="${CUDA_DIR:-$OLCF_CUDA_ROOT}"
    cp ../make.inc.summit make.inc
    export GPU_TARGET=sm70 NV_SM=" -gencode arch=compute_70,code=sm_70" NV_COMP="-gencode arch=compute_70,code=compute_70" MIN_ARCH=350
    MAKE_ARGS=""
  elif [ $IS_DGX2 = YES ] ; then
    cp ../make.inc.summit make.inc
    export CUDA_DIR=$HOME/cuda
    export GPU_TARGET=sm70 NV_SM=" -gencode arch=compute_70,code=sm_70" NV_COMP="-gencode arch=compute_70,code=compute_70" MIN_ARCH=350
    MAKE_ARGS="CC=$HOME/.linuxbrew/bin/gcc-6 CXX=$HOME/.linuxbrew/bin/g++-6"
  elif [ $IS_GPUSYS2 = YES ] ; then
    cp ../make.inc.summit make.inc
    export CUDA_DIR=/usr/local/cuda-10.1
    export GPU_TARGET=sm70 NV_SM=" -gencode arch=compute_75,code=sm_75" NV_COMP="-gencode arch=compute_75,code=compute_75" MIN_ARCH=350
    MAKE_ARGS="CC=$(spack location --install-dir gcc)/bin/g++"
  else
    echo "Unknown platform." 1>&2
    exit 1
  fi

  time make lib $MAKE_ARGS -j8
}

#==============================================================================

function main
{
  do_make "$@" 2>&1 | tee out_make.txt
}

#==============================================================================

main "$@"

#==============================================================================
