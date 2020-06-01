#!/bin/bash
#==============================================================================
#
# Create modified MAGMA versions that are cloned and then patched
# with modified GEMM code.
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function main
{
  local tag
  #for tag in minproduct tally2 tally3 tally4 ; do
  for tag in minproduct mgemm2 mgemm3 mgemm4 mgemm5 ; do
    ./clone_magma.sh $tag
    ./patch_magma.sh $tag
  done
}

#==============================================================================

main "$@"

#==============================================================================
