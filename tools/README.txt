
Genomics calculation workflow
=============================

NOTE: see https://code.ornl.gov/wjd/genomics_gpu/tree/master/tools
for the reference copy of the files in this directory

Preprocessing example
---------------------

ssh rhea

qsub -I -Abif102 -lnodes=1 -lwalltime=8:0:0

export PATH="${PATH}:/lustre/atlas1/bif102/proj-shared/comet"

mkdir 28M_permuted

rhea88$ date ; time shuf 28M_original/28M.txt > 28M_permuted/28M.txt ; date
Tue Jan  9 17:49:40 EST 2018
real	17m20.438s
user	0m53.833s
sys	16m21.303s
Tue Jan  9 18:07:01 EST 2018
[1041s] rhea88$ 

rhea88$ date ; time preprocess 28M_permuted/28M.txt 28M_permuted/28M.bin ; date
Tue Jan  9 18:31:46 EST 2018
real	36m14.392s
user	34m31.482s
sys	1m42.424s
Tue Jan  9 19:08:01 EST 2018
[2175s] rhea88$ 

rhea88$ date ; time allele_labels.sh 28M_permuted/28M.txt 28M_permuted/28M_allele_labels.txt ; date
Tue Jan  9 19:09:20 EST 2018
real	38m43.014s
user	53m8.709s
sys	2m24.280s
Tue Jan  9 19:48:03 EST 2018
[2323s] rhea88$ 

rhea115$ date ; time labels.sh 28M_permuted/28M.txt 28M_permuted/28M_labels.txt ; date
Tue Jan  9 19:30:55 EST 2018
real	11m56.312s
user	7m54.487s
sys	4m42.192s
Tue Jan  9 19:42:52 EST 2018
[717s] rhea115$

rhea114$ date ; time line_indices 28M_permuted/28M.txt 28M_permuted/28M_line_indices.txt ; date
Tue Jan  9 19:35:27 EST 2018
real	16m42.788s
user	15m29.286s
sys	1m13.233s
Tue Jan  9 19:52:10 EST 2018
[1003s] rhea114$ 

OPTIONAL:

rhea114$ date ; time preprocess_validate 28M_permuted/28M.bin 28M_permuted/28M.txt 28M_permuted/28M_allele_labels.txt ; date
Tue Jan  9 19:53:05 EST 2018
Number of elements validated: 3523476076l of 3523476076l
real	22m13.176s
user	21m10.405s
sys	1m1.169s
Tue Jan  9 20:15:18 EST 2018
[1333s] rhea114$ 


How to select CoMet settings, 2-way case
----------------------------------------

When using a GPU-capable system like Titan or Summit, the CoMet code is
set up to use 1 GPU per MPI rank (denoted below as "proc" or "process").
On Titan this means 1 rank per node, on Summit 6 ranks per node.
OpenMP is also used to divide up the CPU cores on the node to
speed up the non-GPU work (16 OpenMP threads per MPI rank on Titan,
7 threads per rank on Summit).

CoMet has several tuning parameters for running in parallel:

  --num_proc_vector: the number of process along the "vector" axis.

  --num_proc_field: the number of process along the "field" axis.
    Note num_proc_vector X num_proc_field = num_proc equals the total
    number of MPI processes used to solve the problem.  num_proc
    nodes should be selected by qsub/aprun on Titan.  On Summit which
    has 6 GPUs per node, one should select ceil(num_proc/6) nodes.

  --num_phase, --phase_min, --phase_max: Since solving an entire
    problem meay require more memory than is available, it is
    possible to divide the computation into smaller "phases,"
    each of which computes only a subset of the result metrics.
    Computing a single phase fully utilizes all GPUs in the
    allocation.  Example: one could set num_phase=200,
    phase_min=0, phase_max=0 to compute the first phase only
    out of 200 phases.  One can num_phase=200, phase_min=0,
    phase_max=199 to compute all phases, while only requiring
    the amount of memory required to compute one phase at a time.
    NOTE: For computing all the metrics for a given fixed problem,
    one must fix the num_proc_vector and num_phase values for all phases
    computed to ensure a consistent definition of which metrics
    are contained in each phase.

The rationale for adjusting the settings is:
  1) to make the number of vectors and vector elements on each GPU as large
     as possible to mximize performance;
  2) to make num_phase large enough to ensure metrics will fit onto memory;
  3) to run at least several phases per individual run to amortize setup
     costs.

The settings advisor tool can be used to assist with determining good
settings (EXPERIMENTAL):

titan-ext1$ ./settings_advisor --help
usage: settings_advisor [-h] --num_node NUM_NODE --platform {Titan,Summit}
                        --num_way {2,3} --metric_type {czekanowski,ccc}
                        --num_vector NUM_VECTOR --num_field NUM_FIELD

Suggest settings for running the CoMet code.

optional arguments:
  -h, --help            show this help message and exit
  --num_node NUM_NODE   The number of compute nodes to run on.
  --platform {Titan,Summit}
                        Platform to be run on (Titan or Summit)
  --num_way {2,3}       dimension of metric to compute (2 for 2-way, 3 for
                        3-way)
  --metric_type {czekanowski,ccc}
                        metric type (czekanowski or ccc)
  --num_vector NUM_VECTOR
                        number of vectors to process
  --num_field NUM_FIELD
                        number of fields to process

titan-ext1$ ./settings_advisor --num_node 6000 --platform Titan --num_vector 28342758 --num_field 882 --metric_type ccc --num_way 2

GPU memory: 1.07430846 GB required out of 5.637144576 GB available: OK.
To fit into memory, suggest using num_phase = 86 (or higher). Then:
CPU memory: 32.317998864 GB required out of 32.641751449 GB available: OK.
Suggest executing at least several phases per run to amortize setup costs.
Local problem size: num_vector_local = 4724, target = 4000.

titan-ext1$ ./settings_advisor --num_node 1000 --platform Titan --num_vector 28342758 --num_field 882 --metric_type ccc --num_way 2

Insufficient GPU memory, please increase num_node.
GPU memory: 38.578422561 GB requested out of 5.637144576 GB available.

NOTE: Currently it is recommended to set --num_field to 1.  A later version
of the tool will give more detailed information on this setting.


CoMet execution example, 2-way case
-----------------------------------

Titan:
ssh titan.ccs.ornl.gov
qsub -Abif102 -lnodes=6000 -lwalltime=0:5:0 mybatchscript.sh

Summit:
ssh summit.olcf.ornl.gov
bsub -P bif102 -Is -nnodes 1000 -W 5 mybatchscript.sh

num_vector=28342758
num_field=882
if [ -n "${CRAYOS_VERSION:-}" ] ; then # Titan
  num_proc_vector=6000
else # Summit
  num_proc_vector=1000
fi
num_phase=200
phase_min=0
phase_max=7
#phase_max=199

out_dir=$PWD/outs_full_${num_proc_vector}_2way_200_000-007
mkdir -p $out_dir

if [ -n "${CRAYOS_VERSION:-}" ] ; then # Titan
  executable=/lustre/atlas1/bif102/proj-shared/comet/genomics_metric_titan
  exec_command="env OMP_NUM_THREADS=16 aprun -cc none -n$num_proc_vector -d16 -N1"
else # Summit
  #TBD
  executable=/lustre/atlas1/bif102/proj-shared/comet/genomics_metric_summit
  exec_command="env OMP_NUM_THREADS=7 jsrun --nrs $num_proc_vector --bind packed:7 --cpu_per_rs 7 --gpu_per_rs 1 --rs_per_host 6 --tasks_per_rs 1"
  module load cuda
fi

time $exec_command $executable \
  --num_field $num_field \
  --num_vector $num_vector \
  --metric_type ccc \
  --sparse yes \
  --all2all yes \
  --compute_method GPU \
  --num_proc_vector $num_proc_vector \
  --num_phase $num_phase --phase_min $phase_min --phase_max $phase_max \
  --input_file $PWD/28M_permuted/28M.bin \
  --threshold .7 \
  --checksum no \
  --output_file_stub $out_dir/out \
  --verbosity 1 \
  | tee ${out_dir}_log.txt

Titan:
ctime 47.856096 ops 0.000000e+00 rate 0.000000e+00 rate/proc 0.000000e+00 cmpout 2.738229e+09 cmp 5.645041e+16 rate 1.179587e+15 rate/proc 1.965978e+11 vctime 0.002566 mctime 7.133783 intime 7.442626 outtime 17.090905 cpumem 1.446527e+10 gpumem 1.074351e+09 tottime 79.595748

# KEY:
#
# ctime - compute time for the metrics computation proper
# vctime - time for vectors setup, teardown
# mctime - time for metrics set up, teardown
# intime - time for read of input file
# outtime - time for preparation of metrics and write to output files
# tottime - total time
# cmpout - number of metrics (vector comparisons) written
# cmp - number of metrics computed
# rate/proc - metrics computed per second per GPU
# cpumem - cpu memory high water mark - max 32 GB
# gpumem - gpu memory high water mark - max 6 GB
#

titan-ext1$ ls -sF $out_dir/*.bin |head
 7344 outs_full_6000_2way_200_000-007/out_0000.bin
 5708 outs_full_6000_2way_200_000-007/out_0001.bin
 5824 outs_full_6000_2way_200_000-007/out_0002.bin
 4404 outs_full_6000_2way_200_000-007/out_0003.bin
 3756 outs_full_6000_2way_200_000-007/out_0004.bin
 5428 outs_full_6000_2way_200_000-007/out_0005.bin
 6424 outs_full_6000_2way_200_000-007/out_0006.bin
 4480 outs_full_6000_2way_200_000-007/out_0007.bin
 3952 outs_full_6000_2way_200_000-007/out_0008.bin
 4652 outs_full_6000_2way_200_000-007/out_0009.bin

titan-ext1$ postprocess_all.sh 2 28M_permuted/28M_allele_labels.txt \
  28M_permuted/28M_labels.txt \
  $out_dir/out_0000.bin

titan-ext1$ head $out_dir/out_0000.txt
318 0 137 0 01:4964015_G 476:10658_A 0.705773
361 0 137 0 18:6183012_C 476:10658_A 0.750155
389 0 181 1 18:10905067_A 08:7856910_T 0.737368
389 0 318 1 18:10905067_A 01:4964015_T 0.771506
389 1 334 0 18:10905067_G 03:1708199_A 0.714580
442 1 389 1 12:8286382_T 18:10905067_G 0.772094
499 1 489 0 01:18324031_G 01:18055010_C 0.802916
614 0 23 1 11:7610195_A 17:11150952_T 0.706933
614 0 137 0 11:7610195_A 476:10658_A 0.908482
614 0 321 1 11:7610195_A 188:26480_T 0.717981


Postprocessing example
----------------------

ssh rhea.ccs.ornl.gov

qsub -Abif102 -lnodes=64 -lwalltime=10:0:0 -I

date ; time postprocess_all.sh 2 28M_permuted/28M_allele_labels.txt \
  28M_permuted/28M_labels.txt \
  $(ls $out_dir/out_????.bin) ; date
Wed Jan 10 21:06:43 EST 2018
real    190m42.712s
user    8m51.149s
sys     2236m21.891s
Thu Jan 11 00:17:26 EST 2018

OPTIONAL:

date; time ccc_validate_all.sh 2 28M_permuted/28M.txt \
  28M_permuted/28M_line_indices.txt \
  $(ls $out_dir/out_????.txt) ; date
Thu Jan 11 18:51:07 EST 2018
COMPLETED outs_full_6000_2way_200_000-007/out_0347.txt
. . .
COMPLETED outs_full_6000_2way_200_000-007/out_5547.txt
real    459m32.097s
user    138m42.621s
sys     5832m1.785s
Fri Jan 12 02:30:39 EST 2018

date; time validate_all.sh 2 $(ls $out_dir/out_????.txt); date
Fri Jan 12 10:25:42 EST 2018
real    0m25.229s
user    2m7.665s
sys     0m11.762s
Fri Jan 12 10:26:07 EST 2018

