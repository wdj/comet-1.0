#!/bin/bash
#==============================================================================
#
# Script for executing 2-way weak scaling study on Summit.
#
# Usage:
#
# env num_node_solve=20   bsub -P stf006 -nnodes 40   -W 60 -J CoMet-benchmark scaling_weak_summit_2way.sh
# env num_node_solve=50   bsub -P stf006 -nnodes 70   -W 60 -J CoMet-benchmark scaling_weak_summit_2way.sh
# env num_node_solve=200  bsub -P stf006 -nnodes 220  -W 60 -J CoMet-benchmark scaling_weak_summit_2way.sh
# env num_node_solve=500  bsub -P stf006 -nnodes 520  -W 60 -J CoMet-benchmark scaling_weak_summit_2way.sh
# env num_node_solve=1000 bsub -P stf006 -nnodes 1020 -W 60 -J CoMet-benchmark scaling_weak_summit_2way.sh
# env num_node_solve=2000 bsub -P stf006 -nnodes 2020 -W 60 -J CoMet-benchmark scaling_weak_summit_2way.sh
# env num_node_solve=3000 bsub -P stf006 -nnodes 3020 -W 60 -J CoMet-benchmark scaling_weak_summit_2way.sh
# env num_node_solve=4400 bsub -P stf006 -nnodes 4420 -W 60 -J CoMet-benchmark scaling_weak_summit_2way.sh
# env num_node_solve=4560 bsub -P stf006 -nnodes 4580 -W 60 -J CoMet-benchmark scaling_weak_summit_2way.sh
#
# Options:
#
# export metric_type=ccc sparse=yes # default
# export metric_type=ccc sparse=yes tc=0
# export metric_type=ccc sparse=no
# export metric_type=ccc sparse=no tc=0
# export metric_type=duo sparse=yes
# export metric_type=czekanowski single=1
# export metric_type=czekanowski single=0
#
# original file: /ccs/home/joubert/proj/genomics/results/gbrev_max
#
#==============================================================================

#------------------------------------------------------------------------------
# Node counts

num_node_job=$(( ( $(echo $LSB_MCPU_HOSTS | wc -w) - 2 ) / 2 ))

[[ -z "${num_node_solve:-}" ]] && num_node_solve=$num_node_job
[[ -z "${num_node_launch:-}" ]] && num_node_launch=$num_node_job

ranks_per_node=6
num_proc_vector=$(( $num_node_solve * $ranks_per_node ))

num_proc_field=1
num_proc_repl=1
num_proc=$(( $num_proc_vector * $num_proc_field * $num_proc_repl ))

#------------------------------------------------------------------------------
# Algorithm settings

[[ "${ccc:-}" = 1 ]] && metric_type=ccc   # legacy settings
[[ "${ccc:-}" = 0 ]] && metric_type=czekanowski   # legacy settings
[[ -z "${metric_type:-}" ]] && metric_type=ccc
[[ -z "${single:-}" ]] && single=1
[[ -z "${sparse:-}" ]] && sparse=yes
[[ -z "${tc:-}" && $metric_type != czekanowski ]] && tc=1
[[ -z "${tc:-}" && $metric_type = czekanowski ]] && tc=0
[[ -z "${cksum:-}" ]] && cksum=yes # alt. cksum=no
[[ -z "${problem_type:-}" ]] && problem_type=random # alt. problem_type=analytic
[[ -z "${debug:-}" ]] && debug=0
[[ -z "${num_tc_steps:-}" ]] && num_tc_steps=4

#------------------------------------------------------------------------------
# Problem sizes

if [ "$metric_type" != czekanowski ] ; then
  num_vector_local=9984
  num_field_local=$(( 98304 * $num_tc_steps ))
elif [ "$single" = 1 ] ; then
  num_vector_local=15000
  num_field_local=56000
else
  num_vector_local=14000
  num_field_local=25000
fi

num_vector=$(( $num_vector_local * $num_proc_vector ))

# Compute one phase of results out of possibly many phases
[[ -z "$num_phase_ratio" ]] && num_phase_ratio=$(( 30 * $num_proc_repl ))
[[ -z "$num_phase" ]] && num_phase=$(( ( $num_proc_vector + $num_phase_ratio - 1 ) / $num_phase_ratio ))
[[ -z "$phase_min" ]] && phase_min=$(( $num_phase - 4 ))
[[ $phase_min < 0 ]] && phase_min=0
[[ -z "$phase_max" ]] && phase_max=$(( $num_phase - 1 ))

#------------------------------------------------------------------------------
# Execution settings

module -q load gcc
module -q load cuda

host=$(echo $(hostname -f) | sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' -e 's/[.-].*//' -e 's/[0-9]*$//')
OLCF_PROJECT="$(echo $LSB_PROJECT_NAME | tr A-Z a-z)"
INSTALLS_DIR=/gpfs/alpine/$OLCF_PROJECT/scratch/$(whoami)/comet

if [ $debug = 1 ] ; then
  executable_double=$INSTALLS_DIR/install_test_$host/bin/genomics_metric
  executable_single=$INSTALLS_DIR/install_single_test_$host/bin/genomics_metric
else
  executable_double=$INSTALLS_DIR/install_release_$host/bin/genomics_metric
  executable_single=$INSTALLS_DIR/install_single_release_$host/bin/genomics_metric
fi

if [ "$metric_type" != czekanowski ] ; then
  executable=$executable_double
  [[ "$sparse" == yes ]] && tag=${metric_type}_sparse || tag=${metric_type}_nonsparse
elif [ "$single" = 1 ] ; then
  executable=$executable_single
  tag=czek_single
else
  executable=$executable_double
  tag=czek_double
fi
[[ $tc != 0 ]] && tag=${tag}_tc

uid=${tag}_$(echo $(( 100000 + $num_proc_field )) | sed 's/.//')
uid=${uid}_$(echo $(( 100000 + $num_proc_vector )) | sed 's/.//')
uid=${uid}_$(echo $(( 100000 + $num_proc_repl )) | sed 's/.//')
uid=${uid}_$(echo $(( 100000 + $num_proc )) | sed 's/.//')
uid=${uid}_${num_field_local}_${num_vector}_${phase_min}_${phase_max}_${num_phase}_$$

# Output file stub
out_stub=out_2way_$uid
outfile=${out_stub}_log.txt

ar_opts="PAMI_IBV_ENABLE_DCT=1 PAMI_ENABLE_STRIPING=1 PAMI_IBV_ADAPTER_AFFINITY=0 PAMI_IBV_QP_SERVICE_LEVEL=8 PAMI_IBV_ENABLE_OOO_AR=1"
launch_command="env OMP_NUM_THREADS=7 $ar_opts jsrun --nrs $(( $num_node_launch * $ranks_per_node )) --bind packed:7 --cpu_per_rs 7 --gpu_per_rs 1 --rs_per_host $ranks_per_node --tasks_per_rs 1 -X 1"

#------------------------------------------------------------------------------

[[ $num_node_solve == $num_node_launch ]] &&  fastnodes_arg="" || fastnodes_arg="--fastnodes"

# Command to execute, with options

if [ $metric_type != czekanowski ] ; then
  exec_command="$launch_command $executable \
    --num_field_local $num_field_local \
    --num_vector_local $num_vector_local \
    --metric_type $metric_type \
    --sparse $sparse \
    --all2all yes \
    --compute_method GPU \
    --problem_type $problem_type \
    --checksum $cksum \
    --num_proc_vector $num_proc_vector \
    --num_proc_field $num_proc_field \
    --num_proc_repl $num_proc_repl \
    --num_phase $num_phase --phase_min $phase_min --phase_max $phase_max \
    --threshold .7 \
    --verbosity 1 $fastnodes_arg \
    --tc $tc --num_tc_steps $num_tc_steps "
else
  exec_command="$launch_command $executable \
    --num_field_local $num_field_local \
    --num_vector_local $num_vector_local \
    --metric_type czekanowski \
    --all2all yes \
    --compute_method GPU \
    --problem_type $problem_type \
    --checksum $cksum \
    --num_proc_vector $num_proc_vector \
    --num_proc_field $num_proc_field \
    --num_proc_repl $num_proc_repl \
    --num_phase $num_phase --phase_min $phase_min --phase_max $phase_max \
    --verbosity 1 $fastnodes_arg"
fi

# Perform run

date | tee -a $outfile
echo "$exec_command" | tee -a $outfile
time $exec_command 2>&1 | tee -a $outfile
date | tee -a $outfile

exit

#==============================================================================
