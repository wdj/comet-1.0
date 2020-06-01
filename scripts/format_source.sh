#!/bin/bash
#==============================================================================
#
# Source code formatter.
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function main
{
  module load clang

  type clang-format

  if [ $? != 0 ] ; then
    exit 1
  fi

  for dir in src driver testing ; do
    for suffix in h cc ; do
      for file in $(find $dir -name '*'.$suffix -print) ; do
        file_bak="${file}.bak"
        if [ ! -e "$file_bak" ] ; then
          cp "$file" "$file_bak"
          clang-format -style chromium -i "$file"
          echo "$file"
        fi
      done
    done
  done
}

#==============================================================================

main "$@"

#==============================================================================
