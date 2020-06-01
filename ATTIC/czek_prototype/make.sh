#!/bin/bash -l

module swap PrgEnv-pgi PrgEnv-gnu
module load cudatoolkit
module load acml

MAGMA_DIR=../magma/magma_minproduct-1.6.2

if [ ! -e $MAGMA_DIR/lib/libmagma_minproduct.a ] ; then
  #---Build modified magma library.
  pushd $MAGMA_DIR
  #cp /sw/xk6/magma/1.3/cle4.0_gnu4.7.2_cuda5.0_acml5.2.0/source/make.inc .
  #cp /sw/xk6/magma/1.6.2/sles11.3_gnu4.8.2/source/make.inc .
  make -j8
  popd
fi

make distclean
make 2>&1

