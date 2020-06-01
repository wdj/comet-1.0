#!/bin/bash
#==============================================================================
#
# Configure ALL versions prior to build.
# This should be executed in the directory containing the genomics_gpu
# repo directory.
#
# Relevant input variables:
#
# OLCF_PROJECT - OLCF project ID
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
  if [ -z "${OLCF_PROJECT:-}" ] ; then
    local OLCF_PROJECT=stf006
  fi
  local host
  host=$(echo $(hostname -f) | sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' -e 's/[.-].*//' -e 's/[0-9]*$//')
  local DIRNAME_STUB
  [[ $IS_EXPERIMENTAL = YES ]] && DIRNAME_STUB=experimental || DIRNAME_STUB=$host
  #
  if [ $IS_EXPERIMENTAL = YES ] ; then
    true # skip for now
  elif [ $IS_CRAY_XK7 = YES ] ; then
    local INSTALLS_DIR=/lustre/atlas/scratch/$(whoami)/$OLCF_PROJECT/comet
  elif [ $IS_IBM_AC922 = YES ] ; then
    local INSTALLS_DIR=/gpfs/alpine/$OLCF_PROJECT/scratch/$(whoami)/comet
  elif [ $IS_EDISON = YES ] ; then
    local INSTALLS_DIR="$SCRATCH/comet"
  else
    local INSTALLS_DIR="$PWD/installs"
    #echo "Unknown platform." 1>&2
    #exit 1
  fi
  mkdir -p "$INSTALLS_DIR"

  #----------------------------------------------------------------------------
  # test / double precision build

  local DO_BUILD_TEST=YES # NO
  [[ $IS_DGX2 = YES ]] && DO_BUILD_TEST=NO
  [[ $IS_GPUSYS2 = YES ]] && DO_BUILD_TEST=NO
  if [ $DO_BUILD_TEST = YES ] ; then
    local BUILD_DIR=build_test_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma_patch # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_test_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR BUILD_TYPE=Debug TESTING=ON \
        ../genomics_gpu/scripts/cmake.sh
    if [ -e magma_patch -a ! -e ../magma_build_$DIRNAME_STUB ] ; then
      mv magma_patch ../magma_build_$DIRNAME_STUB # share common MAGMA build
      ln -s          ../magma_build_$DIRNAME_STUB magma_patch
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi

  #----------------------------------------------------------------------------
  # test / double precision / nompi build

  local DO_BUILD_TEST_NOMPI=YES # NO
  if [ $DO_BUILD_TEST_NOMPI = YES ] ; then
    local BUILD_DIR=build_test_nompi_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma_patch # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_test_nompi_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR BUILD_TYPE=Debug NOMPI=ON \
        ../genomics_gpu/scripts/cmake.sh
    if [ -e magma_patch -a ! -e ../magma_build_$DIRNAME_STUB ] ; then
      mv magma_patch ../magma_build_$DIRNAME_STUB # share common MAGMA build
      ln -s          ../magma_build_$DIRNAME_STUB magma_patch
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi

  #----------------------------------------------------------------------------
  # test / single precision build

  local DO_BUILD_SINGLE_TEST=YES # NO
  [[ $IS_DGX2 = YES ]] && DO_BUILD_SINGLE_TEST=NO
  [[ $IS_GPUSYS2 = YES ]] && DO_BUILD_SINGLE_TEST=NO
  if [ $DO_BUILD_SINGLE_TEST = YES ] ; then
    local BUILD_DIR=build_single_test_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma_patch # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_single_test_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR FP_PRECISION=SINGLE BUILD_TYPE=Debug TESTING=ON \
        ../genomics_gpu/scripts/cmake.sh
    if [ -e magma_patch -a ! -e ../magma_build_$DIRNAME_STUB ] ; then
      mv magma_patch ../magma_build_$DIRNAME_STUB # share common MAGMA build
      ln -s          ../magma_build_$DIRNAME_STUB magma_patch
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi

  #----------------------------------------------------------------------------
  # release / double precision build

  local DO_BUILD_RELEASE=YES # NO
  [[ $IS_DGX2 = YES ]] && DO_BUILD_RELEASE=NO
  [[ $IS_GPUSYS2 = YES ]] && DO_BUILD_RELEASE=NO
  if [ $DO_BUILD_RELEASE = YES ] ; then
    local BUILD_DIR=build_release_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma_patch # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_release_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR BUILD_TYPE=Release \
        ../genomics_gpu/scripts/cmake.sh
    if [ -e magma_patch -a ! -e ../magma_build_$DIRNAME_STUB ] ; then
      mv magma_patch ../magma_build_$DIRNAME_STUB # share common MAGMA build
      ln -s          ../magma_build_$DIRNAME_STUB magma_patch
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi

  #----------------------------------------------------------------------------
  # release / double precision / nompi build

  local DO_BUILD_RELEASE_NOMPI=YES # NO
  if [ $DO_BUILD_RELEASE_NOMPI = YES ] ; then
    local BUILD_DIR=build_release_nompi_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma_patch # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_release_nompi_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR BUILD_TYPE=Release NOMPI=ON \
        ../genomics_gpu/scripts/cmake.sh
    if [ -e magma_patch -a ! -e ../magma_build_$DIRNAME_STUB ] ; then
      mv magma_patch ../magma_build_$DIRNAME_STUB # share common MAGMA build
      ln -s          ../magma_build_$DIRNAME_STUB magma_patch
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi

  #----------------------------------------------------------------------------
  # release / single precision build

  local DO_BUILD_SINGLE_RELEASE=YES # NO
  [[ $IS_DGX2 = YES ]] && DO_BUILD_SINGLE_RELEASE=NO
  [[ $IS_GPUSYS2 = YES ]] && DO_BUILD_SINGLE_RELEASE=NO
  if [ $DO_BUILD_SINGLE_RELEASE = YES ] ; then
    local BUILD_DIR=build_single_release_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma_patch # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_single_release_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR FP_PRECISION=SINGLE BUILD_TYPE=Release \
        ../genomics_gpu/scripts/cmake.sh
    if [ -e magma_patch -a ! -e ../magma_build_$DIRNAME_STUB ] ; then
      mv magma_patch ../magma_build_$DIRNAME_STUB # share common MAGMA build
      ln -s          ../magma_build_$DIRNAME_STUB magma_patch
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi
}

#==============================================================================

main "$@"

#==============================================================================
