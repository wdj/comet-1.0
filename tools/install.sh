#!/bin/bash
#------------------------------------------------------------------------------

INSTALL_DIR="/lustre/atlas1/bif102/proj-shared/comet"

for i in preprocess preprocess_validate postprocess_file line_indices ccc_validate ; do
  cp $i $INSTALL_DIR
done

cp *.sh README.txt $INSTALL_DIR

chgrp -R bif102 $INSTALL_DIR/*
chmod -R g+rX $INSTALL_DIR/*

#------------------------------------------------------------------------------
