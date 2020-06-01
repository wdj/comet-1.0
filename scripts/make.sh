#!/bin/bash -l
#==============================================================================
#
# Script to build genomics metrics code.
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function main
{
  local IS_CRAY_XK7 # OLCF Titan or Chester
  [[ -n "${CRAYOS_VERSION:-}" ]] && IS_CRAY_XK7="YES" || IS_CRAY_XK7="NO"
  local IS_IBM_AC922 # OLCF Summit or Peak
  [[ -n "${LSF_BINDIR:-}" ]] && IS_IBM_AC922="YES" || IS_IBM_AC922="NO"
  local IS_DGX2
  [[ "$(uname -n)" = "dgx2-b" ]] && IS_DGX2="YES" || IS_DGX2="NO"
  local IS_GPUSYS2
  [[ "$(uname -n)" = "gpusys2" ]] && IS_GPUSYS2="YES" || IS_GPUSYS2="NO"
  local IS_EDISON
  [[ "${NERSC_HOST:-}" = "edison" ]] && IS_EDISON="YES" || IS_EDISON="NO"
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
    module list
  elif [ $IS_IBM_AC922 = YES ] ; then
    module -q load gcc/6.4.0
    local CUDA_MODULE=cuda
    module -q load $CUDA_MODULE
    module list
  elif [ $IS_DGX2 = YES ] ; then
    true
  elif [ $IS_GPUSYS2 = YES ] ; then
    true
  elif [ $IS_EDISON = YES ] ; then
    module swap PrgEnv-intel PrgEnv-gnu
  else
    echo "Unknown platform." 1>&2
    exit 1
  fi

  time make -j4 VERBOSE=1

  if [ $? = 0 ] ; then
    time make install
    exit $?
  else
    exit $?
  fi
}

#==============================================================================

main "$@"

#==============================================================================
