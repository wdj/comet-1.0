#!/bin/bash -l

NPROC=1

for PROCESS_VECTORS_FN in process_vectors_alt ; do
#for NCOPIES_V in 1 2 4 8 16 32 ; do
#for NCOPIES_F in 1 2 4 8 16 32 ; do
for NCOPIES_V in 32 ; do
for NCOPIES_F in 32 ; do

  #if [ $(( $NCOPIES_V * $NCOPIES_F * $NCOPIES_F )) -le 4096 -a \( $NCOPIES_V != 32 -o $NCOPIES_F != 8 \) ] ; then
  #  continue
  #fi

  export PROCESS_VECTORS_FN NCOPIES_V NCOPIES_F

  ./make.sh 2>&1 >/dev/null

  WD=$PWD
  #cp -p czek $MEMBERWORK/stf006/czek.exe
  pushd $MEMBERWORK/stf006 > /dev/null
  echo -n "PROCESS_VECTORS_FN $PROCESS_VECTORS_FN "
  echo -n "NCOPIES_V $NCOPIES_V "
  echo -n "NCOPIES_F $NCOPIES_F "
  OUTFILE=$WD/test.output.${PROCESS_VECTORS_FN}.${NCOPIES_V}.${NCOPIES_F}.txt
  aprun -n$NPROC -N1 $WD/czek <$WD/test.input.txt 2>&1 > $OUTFILE
  grep -v BESC <$OUTFILE | grep ' time'
  popd > /dev/null

done
done
done

