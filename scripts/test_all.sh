#!/bin/bash
#==============================================================================
#
# Run code tests on ALL versions.
# This should be executed in the directory containing the genomics_gpu
# repo directory.
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

# qsub -Astf006 -lnodes=4 -lwalltime=2:0:0 -I
# bsub -P stf006 -Is -nnodes 2 -alloc_flags gpumps -W 120 $SHELL

#==============================================================================

function main
{
  local IS_CRAY_XK7 # OLCF Titan or Chester
  [[ -n "${CRAYOS_VERSION:-}" ]] && IS_CRAY_XK7="YES" || IS_CRAY_XK7="NO"
  local IS_IBM_AC922 # OLCF Summit or Peak
  [[ -n "${LSF_BINDIR:-}" ]] && IS_IBM_AC922="YES" || IS_IBM_AC922="NO"
  local IS_EXPERIMENTAL
  [[ "${COMET_BUILD_EXPERIMENTAL:-}" = YES ]] && IS_EXPERIMENTAL="YES" || \
                                                 IS_EXPERIMENTAL="NO"

  local host
  host=$(echo $(hostname -f) | sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' -e 's/[.-].*//' -e 's/[0-9]*$//')
  local DIRNAME_STUB
  [[ $IS_EXPERIMENTAL = YES ]] && DIRNAME_STUB=experimental || DIRNAME_STUB=$host

  # WARNING: these module loads may need to match those in scripts/cmake.sh
  if [ $IS_EXPERIMENTAL = YES ] ; then
    true # skip for now
  elif [ $IS_CRAY_XK7 = YES ] ; then
    if [ "$PE_ENV" = "PGI" ] ; then
      module unload PrgEnv-pgi
    fi
    module load PrgEnv-gnu
    module load cudatoolkit
  elif [ $IS_IBM_AC922 = YES ] ; then
    module -q load gcc/6.4.0
    local CUDA_MODULE=cuda
    module -q load $CUDA_MODULE
  else
    echo "Unknown platform." 1>&2
    exit 1
  fi

  local dirs="build_test_$DIRNAME_STUB build_single_test_$DIRNAME_STUB"

  local dir
  for dir in $dirs ; do
    echo "===================="
    echo $dir
    echo "===================="
    pushd $dir
    time make test ARGS=-V 2>&1 | tee out_test.txt
    if [ $? != 0 ] ; then
      exit 1
    fi
    popd
  done

  echo "-------------------------------------------------------------------------------"
  for dir in $dirs ; do
    grep -H fail $dir/out_test.txt
  done
  out_files="$(for dir in $dirs ; do echo $dir/out_test.txt ; done)"
  if [ $(grep ' tests fail' $out_files <(echo ' tests fail') | wc -l) = \
       $(grep ' 0 tests fail' $out_files <(echo ' 0 tests fail') | wc -l) ] ; then
    echo "!!! All tests PASSED !!!"
  fi

  echo "-------------------------------------------------------------------------------"
} # main

#==============================================================================

main "$@"

#==============================================================================
