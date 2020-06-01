#!/bin/bash
#==============================================================================
#
# Infer a patch file that can be applied to a cloned (unpatched) MAGMA dir to
# get the cloned and patched Magma dir.
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function main
{
  local tag="${1:-}"
  if [ "$tag" = "" -o "$tag" = "-h" -o "$tag" = "--help" ] ; then 
    echo "Usage: ${0##*/} <tag>"
    echo "  where <tag> is an alphanumeric string. magma_<tag> is the new library name."
    exit 0
  fi

  local dir_cloned="magma_${tag}.cloned"
  local dir_cloned_patched="magma_${tag}"
  local patch_file=magma_${tag}.patch

  diff -aur $dir_cloned $dir_cloned_patched | grep -v "^Only in " > $patch_file
}

#==============================================================================

main "$@"

#==============================================================================
