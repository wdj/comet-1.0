#!/bin/bash
#------------------------------------------------------------------------------

INSTALL_DIR="/lustre/atlas1/bif102/proj-shared/comet"

for i in preprocess preprocess_validate postprocess_file line_indices ccc_validate ; do
  gcc -std=c99 -o $i ${i}.c
done

#------------------------------------------------------------------------------
