#!/bin/bash -l
#------------------------------------------------------------------------------

# qsub -I -Astf006 -lnodes=1 -lwalltime=2:0:0
# ...

NPROC=1

for PROCESS_VECTORS_FN in process_vectors_alt6 process_vectors_alt ; do
#for NCOPIES_V in 1 2 4 8 16 32 ; do
#for NCOPIES_F in 1 2 4 8 16 32 ; do
for NCOPIES_V in 16 ; do
for NCOPIES_F in 16 ; do

  #---skip cases that are likely not to have enough memory.
  if [ $PROCESS_VECTORS_FN = process_vectors_alt -a \
       $(( $NCOPIES_V * $NCOPIES_F * $NCOPIES_F )) -gt 4096 ] ; then
    continue
  fi

  export PROCESS_VECTORS_FN NCOPIES_V NCOPIES_F

  #---Build code with these hardwired settings.
  #---TODO: make these settings command line arguments to avoid recompile.
  ./make.sh 2>&1 >/dev/null

  WD=$PWD
  #cp -p czek $MEMBERWORK/stf006/czek.exe
  #---cd to Lustre to be able to aprun on titan.
  pushd $MEMBERWORK/stf006 > /dev/null
  echo -n "PROCESS_VECTORS_FN $PROCESS_VECTORS_FN "
  echo -n "NCOPIES_V $NCOPIES_V "
  echo -n "NCOPIES_F $NCOPIES_F "
  #---Perform run, show the significant output.
  OUTFILE=$WD/test.output.${PROCESS_VECTORS_FN}.${NCOPIES_V}.${NCOPIES_F}.txt
  aprun -n$NPROC -N1 $WD/czek <$WD/test.input.txt 2>&1 > $OUTFILE 
  grep -v BESC <$OUTFILE | grep ' time'
  popd > /dev/null

done
done
done

#------------------------------------------------------------------------------
