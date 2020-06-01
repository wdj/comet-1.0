#!/bin/bash
#==============================================================================
#
# Apply a patch file to to the cloned (unpatched) MAGMA dir to get
# the cloned and patched Magma dir.
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

  cp -rp $dir_cloned $dir_cloned_patched
  patch -p0 <$patch_file
}
#==============================================================================

main "$@"

#==============================================================================
