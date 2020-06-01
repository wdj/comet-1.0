#!/bin/bash -l
#==============================================================================
#
# Script to run test cases of code.
#
# Usage:
#
# qsub -I -Astf006 -lnodes=1 -lwalltime=2:0:0
# ...
# cd genomics_gpu
# cd ../build
# ../genomics_gpu/scripts/run.sh
#
#==============================================================================

ACCOUNT=stf006

WD=$PWD
#EXEC=$WD/../install_debug/bin/genomics_metric
EXEC=$WD/genomics_metric

#---cd to Lustre to be able to aprun on titan.
pushd $MEMBERWORK/$ACCOUNT > /dev/null

#for compute_method in CPU GPU ; do
for compute_method in CPU GPU ; do

  #aprun -n3 -N1 $EXEC \
  #    --num_field 1 --num_vector_local 2 \
  #    --compute_method $compute_method --verbosity 2

#  aprun -n1 -N1 $EXEC \
#      --num_field 1 --num_vector_local 2 \
#      --compute_method $compute_method --verbosity 2
#  aprun -n4 -N1 $EXEC \
#      --num_field 1 --num_vector_local 2 --num_proc 1 \
#      --compute_method $compute_method --verbosity 2



  #echo \
  #aprun -n1 $EXEC \
  #    --num_field 5000 --num_vector_local 6000 \
  #    --compute_method $compute_method --verbosity 1
  #aprun -n1 $EXEC \
  #    --num_field 5000 --num_vector_local 6000 \
  #    --compute_method $compute_method --verbosity 1

#  echo \
#  aprun -n1 $EXEC \
#      --num_field 50 --num_vector_local 60 --num_way 3 \
#      --compute_method $compute_method --verbosity 1
#  aprun -n1 $EXEC \
#      --num_field 50 --num_vector_local 60 --num_way 3 \
#      --compute_method $compute_method --verbosity 1

#  for i in {1..2} ; do
#  aprun -n$i $EXEC \
#      --num_field 1 --num_vector_local $(( 6 / $i )) --num_way 2 \
#      --compute_method $compute_method --verbosity 2 --all2all yes
#  done


#  aprun -n1 -N1 $EXEC \
#      --num_field 1 --num_vector_local 4 --all2all yes --num_way 3 \
#      --compute_method $compute_method --verbosity 2
#
#  aprun -n1 -N1 $EXEC \
#      --num_field 1 --num_vector_local 4 --all2all no --num_way 3 \
#      --compute_method $compute_method --verbosity 2


#  aprun -n1 -N1 $EXEC \
#      --num_field 1 --num_vector_local 6 --all2all no --num_way 3 \
#      --compute_method $compute_method --verbosity 2
#
#  aprun -n1 -N1 $EXEC \
#      --num_field 1 --num_vector_local 6 --all2all yes --num_way 3 \
#      --compute_method $compute_method --verbosity 2
#
#  aprun -n2 -N1 $EXEC \
#      --num_field 1 --num_vector_local 3 --all2all yes --num_way 3 \
#      --compute_method $compute_method --verbosity 2



  aprun -n1 -N1 $EXEC \
      --num_field 1 --num_vector_local 18 --all2all no --num_way 3 \
      --compute_method $compute_method --verbosity 1

  aprun -n1 -N1 $EXEC \
      --num_field 1 --num_vector_local 18 --all2all yes --num_way 3 \
      --compute_method $compute_method --verbosity 1
  aprun -n2 -N1 $EXEC \
      --num_field 1 --num_vector_local 9 --all2all yes --num_way 3 \
      --compute_method $compute_method --verbosity 1
  aprun -n3 -N1 $EXEC \
      --num_field 1 --num_vector_local 6 --all2all yes --num_way 3 \
      --compute_method $compute_method --verbosity 1


  #for nproc in {1..4} ; do
  #  num_vector_local=$(( 6000 / $nproc ))
  #  aprun -n$nproc -N1 $EXEC \
  #    --num_field 500 --num_vector_local $num_vector_local --all2all yes \
  #    --compute_method $compute_method --verbosity 1
  #done

#  for log2_num_vector_local in 14 13 12 11 ; do
#    num_vector_local=$(( 1 << $log2_num_vector_local ))
#  for log2_num_field in {11..11} ; do
#    num_field=$(( 1 << $log2_num_field ))
#
#    log2_nproc=$(( 2 * ( 14 - $log2_num_vector_local ) ))
#    nproc=$(( 1 << $log2_nproc ))
#
#    echo \
#    time aprun -n$nproc -N1 $EXEC \
#      --num_field $num_field --num_vector_local $num_vector_local \
#      --all2all yes --compute_method $compute_method --verbosity 1
#    time aprun -n$nproc -N1 $EXEC \
#      --num_field $num_field --num_vector_local $num_vector_local \
#      --all2all yes --compute_method $compute_method --verbosity 1
#  done
#  done

  #for nproc_sqrt in 1 2 4 8  ; do
  #  nproc=$(( $nproc_sqrt * $nproc_sqrt ))
  #  num_field=3072
  #  num_field=30
  #  num_vector_local=$(( 16384 / $nproc_sqrt ))
  #  time aprun -n$nproc -N1 $EXEC \
  #    --num_field $num_field --num_vector_local $num_vector_local \
  #    --all2all yes --compute_method $compute_method --verbosity 1
  #done

  #aprun -n4 -N1 $EXEC \
  #    --num_field 1 --num_vector_local 1 --all2all yes \
  #    --compute_method $compute_method --verbosity 2

done

popd > /dev/null

#==============================================================================
