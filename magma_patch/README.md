
Modified MAGMA code
===================

CoMet relies on versions of MAGMA that are slightly modified to support
alternative GEMM-like operations.
This directory contains the required patch files and build scripts.

To create the modified MAGMA source libraries, type:

./create_modified_magmas.sh.

To compile, type:

for tag in minproduct tally2 tally3 tally4 ; do
  pushd magma_$tag;
  ../make_magma.sh;
  popd;
done

These operations are done automatically when the configure/make
of CoMet is done.

FOR DEVELOPERS:
The adaptation of MAGMA is a two step process.
First, running ./clone_magma.sh $tag creates a modified MAGMA version,
magma_${tag}.cloned, with functions and files renamed to avoid namespace
collisions.
Second, this is patched by typing ./patch_magma.sh $tag which applies
the patch file magma_${tag}.patch containing the needed changes to
the source code, creating magma_$tag.
If changes are made to the magma_$tagn directory, the patch file can be
updated for committing to the repo by running ./infer_patch.sh $tag.

