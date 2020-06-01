//-----------------------------------------------------------------------------
/*!
 * \file   genomics_metric.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code for genomics metric calculation.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#include "cstdio"
#include "cstdlib"
#include "cstddef"
#include "string.h"

#include "execinfo.h"
#include "signal.h"
#include "unistd.h"
#include "sys/wait.h"

#include "mpi.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics.hh"

#include "driver.hh"
#include "input_output.hh"

//=============================================================================
/* Stack tracing code */

#if 0
// http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes */

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}
#endif

// http://stackoverflow.com/questions/4636456/how-to-get-a-stack-trace-for-c-using-gcc-with-line-number-information
// http://stackoverflow.com/questions/3151779/how-its-better-to-invoke-gdb-from-program-to-print-its-stacktrace

void bt_sighandler(int sig, struct sigcontext ctx) {

  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

  char pid_buf[30];
  sprintf(pid_buf, "%d", getpid());
  char name_buf[512];
  name_buf[readlink("/proc/self/exe", name_buf, 511)] = 0;
  int child_pid = fork();
  if (! child_pid) {           
    dup2(2,1); // redirect output to stderr
    fprintf(stdout,"MPI rank %i: stack trace for %s pid=%s\n",
            proc_num, name_buf, pid_buf);

    //FIX
    //char cmd[512];
    //sprintf(cmd, "gdb --batch -n -ex thread -ex bt %s %s 2>&1"
    //             " | grep '^#' | sed -e 's/^/MPI rank %i: /'",
    //        name_buf, pid_buf, proc_num);
    //system(cmd);

    abort(); /* If gdb failed to start */
  } else {
    waitpid(child_pid,NULL,0);
  }

  exit(0);
}

//-----------------------------------------------------------------------------

void install_handler() {

  //signal(SIGSEGV, handler);

  struct sigaction sa;

  sa.sa_handler = (__sighandler_t)bt_sighandler;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART;

  // http://www.comptechdoc.org/os/linux/programming/linux_pgsignals.html

  sigaction(SIGSEGV, &sa, NULL);
  sigaction(SIGFPE, &sa, NULL);
  sigaction(SIGINT, &sa, NULL);
  sigaction(SIGILL, &sa, NULL);
  sigaction(SIGUSR1, &sa, NULL);
  sigaction(SIGINT, &sa, NULL);
}

//=============================================================================
/*---Inform user of usage of command---*/

void usage() {
  /* clang-format off */
  printf(
  "genomics_metric: calculation of comparison metrics from genomics data\n"
  "\n"
  "Usage:\n"
  "\n"
  "    genomics_metric <option> ...\n"
  "\n"
  "Options:\n"
  "\n"
  "    --num_field <value>\n"
  "        (Required) the number of elements in each vector\n"
  "\n"
  "    --num_field_local <value>\n"
  "        (Required) the number of elements in each vector on each processor\n"
  "\n"
  "    --num_vector <value>\n"
  "        (Required) the number of vectors to be processed"
  "\n"
  "    --num_vector_local <value>\n"
  "        (Required) the number of vectors to be processed on each processor\n"
  "\n"
  "    --metric_type <value>\n"
  "        metric type to compute (czekanowski=Czekanowski (default),\n"
  "        ccc=CCC, duo=DUO)\n"
  "\n"
  "    --ccc_multiplier <value>\n"
  "        fixed front multiplier value used to calculate the CCC metric\n"
  "        (default floating point value is 9/2).\n"
  "\n"
  "    --duo_multiplier <value>\n"
  "        fixed front multiplier value used to calculate the DUO metric\n"
  "        (default floating point value is 4).\n"
  "\n"
  "    --ccc_param <value>\n"
  "        fixed coefficient value used to calculate the CCC or DUO metric\n"
  "        (default floating point value is 2/3).\n"
  "\n"
  "    --sparse <value>\n"
  "        for CCC and DUO, interpret vector entries of binary \"10\"\n"
  "        as empty or incomplete data (yes=yes, no=no (default))\n"
  "\n"
  "    --num_way <value>\n"
  "        dimension of metric to compute (2=2-way (default), 3=3-way)\n"
  "\n"
  "    --all2all <value>\n"
  "        whether to perform global all-to-all rather than computing\n"
  "        on each processor separately (yes=yes, no=no (default))\n"
  "\n"
  "    --compute_method <value>\n"
  "        manner of computing the result (CPU=cpu, GPU=gpu (default),\n"
  "        REF=reference method)\n"
  "\n"
  "    --tc <value>\n"
  "        for GPU method on platforms that support it, use the GPU\n"
  "        tensor cores for CCC and DUO (0=no (default), 1=yes)\n"
  "\n"
  "    --num_tc_steps <value>\n"
  "        for tensor core method, tuning parameter to reduce memory\n"
  "        usage on the GPU (default 1)\n"
//"\n"
//"    --fastnodes\n"
//"        experimental method to preselect the fastest nodes to run on\n"
  "\n"
  "    --num_proc_vector <value>\n"
  "        blocking factor to denote number of blocks used to decompose\n"
  "        the total number of vectors across processors \n"
  "        (default is the total number of procs requested)\n"
  "\n"
  "    --num_proc_field <value>\n"
  "        blocking factor to denote number of blocks used to decompose\n"
  "        each vector across processors (default is 1)\n"
  "\n"
  "    --num_proc_repl <value>\n"
  "        processor replication factor.  For each block along the vector\n"
  "        and field axes, this number of processors is applied to\n"
  "        computations for the block (default is 1)\n"
  "\n"
  "    --num_stage <value>\n"
  "        the number of stages the computation is divided into\n"
  "        (default is 1) (available for 3-way case only)\n"
  "\n"
  "    --stage_min <value>\n"
  "        the lowest stage number of the sequence of stages to be computed\n"
  "        for this run (0-based, default is 0)\n"
  "\n"
  "    --stage_max <value>\n"
  "        the highest stage number of the sequence of stages to be computed\n"
  "        for this run (0-based, default is num_stage-1)\n"
  "\n"
  "    --num_phase <value>\n"
  "        the number of phases the computation is divided into\n"
  "        (default is 1)\n"
  "\n"
  "    --phase_min <value>\n"
  "        the lowest phase number of the sequence of phases to be computed\n"
  "        for this run (0-based, default is 0)\n"
  "\n"
  "    --phase_max <value>\n"
  "        the highest phase number of the sequence of phases to be computed\n"
  "        for this run (0-based, default is num_phase-1)\n"
  "\n"
  "    --input_file <value>\n"
  "        string denoting the filename or pathname file\n"
  "        used to store input vectors.  If this option not present,\n"
  "        a synthetic test case is run.\n"
  "\n"
  "    --problem_type <value>\n"
  "        the kind of synthetic test case to run. Allowed choices are\n"
  "        analytic (default) or random\n"
  "\n"
  "    --output_file_stub <value>\n"
  "        string denoting the filename or pathname stub of filenames\n"
  "        used to store result metrics.  Metric values are stored in files\n"
  "        whose names are formed by appending a unique identifier\n"
  "        (e.g., processor number) to the end of this string.  If this\n"
  "        option is absent, no output files are written.\n"
  "\n"
  "    --threshold <value>\n"
  "        output each result value only if its magnitude is greater than\n"
  "        this threshold.  If set negative, no thresholding is done\n"
  "        (default -1)\n"
  "\n"
  "    --checksum <value>\n"
  "        compute checksum of the metrics results (yes=yes (default), no=no)\n"
  "\n"
  "    --verbosity <value>\n"
  "       verbosity level of output (0=none, 1=some (default) 2,3=more)\n"
  "\n"
  );
  /* clang-format on */
}

//=============================================================================

int get_node_id() {

  //int proc_num = 0;
  //MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

  //char name[MPI_MAX_PROCESSOR_NAME];
  //int len = MPI_MAX_PROCESSOR_NAME;
  //MPI_Get_processor_name(name, &len);

  const size_t len = 256;
  char name[len];
  gethostname(name, len);

  // Ignore trailing nonnumeric chars

  for (size_t i=len-1; i>=0; --i) {
    int c = name[i];
    if (c >= '0' && c <= '9') {
      break;
    }
    name[i] = 0;
  }

  int node_id = 1;
  for (size_t i=0; i<len; ++i) {
    if (! name[i]) {
      break;
    }
    int c = name[i];
    if (c >= '0' && c <= '9') {
      node_id = node_id * 10 + (c - '0');
    }
    if (c >= 'a' && c <= 'z') {
      node_id = node_id * 26 + (c - 'a');
    }
  }

  return node_id;
}

//=============================================================================

double bad_node_penalty() {
  const size_t len = 256;
  char name[len];
  gethostname(name, len);

  return false ? 1

       : strcmp(name, "a04n08") == 0 ? 1e3 - .01
       : strcmp(name, "a07n16") == 0 ? 1e3 - .02
       : strcmp(name, "a13n10") == 0 ? 1e3 - .03
       : strcmp(name, "a19n04") == 0 ? 1e3 - .04
       : strcmp(name, "a19n07") == 0 ? 1e3 - .05
       : strcmp(name, "a20n03") == 0 ? 1e3 - .06
       : strcmp(name, "a22n17") == 0 ? 1e3 - .07
       : strcmp(name, "a24n09") == 0 ? 1e3 - .08
       : strcmp(name, "a26n12") == 0 ? 1e3 - .09
       : strcmp(name, "a34n15") == 0 ? 1e3 - .10
       : strcmp(name, "b07n03") == 0 ? 1e3 - .11
       : strcmp(name, "b12n13") == 0 ? 1e3 - .12
       : strcmp(name, "b14n03") == 0 ? 1e3 - .13
       : strcmp(name, "b31n08") == 0 ? 1e3 - .14
       : strcmp(name, "c27n10") == 0 ? 1e3 - .15
       : strcmp(name, "c29n14") == 0 ? 1e3 - .16
       : strcmp(name, "d01n18") == 0 ? 1e3 - .17
       : strcmp(name, "d16n13") == 0 ? 1e3 - .18
       : strcmp(name, "d20n01") == 0 ? 1e3 - .19
       : strcmp(name, "d24n18") == 0 ? 1e3 - .20
       : strcmp(name, "e01n05") == 0 ? 1e3 - .21
       : strcmp(name, "e08n07") == 0 ? 1e3 - .22
       : strcmp(name, "e08n15") == 0 ? 1e3 - .23
       : strcmp(name, "f17n01") == 0 ? 1e3 - .24
       : strcmp(name, "f25n02") == 0 ? 1e3 - .25
       : strcmp(name, "f31n14") == 0 ? 1e3 - .26
       : strcmp(name, "g03n17") == 0 ? 1e3 - .27
       : strcmp(name, "g26n14") == 0 ? 1e3 - .28
       : strcmp(name, "g33n09") == 0 ? 1e3 - .29
       : strcmp(name, "h12n01") == 0 ? 1e3 - .30
       : strcmp(name, "h14n09") == 0 ? 1e3 - .31
       : strcmp(name, "h24n18") == 0 ? 1e3 - .32

       : strcmp(name, "h25n12") == 0 ? 1e6 - 1 // 4X slower
       : strcmp(name, "d06n12") == 0 ? 1e6 - 1 // 2X slower - once
       : strcmp(name, "b08n02") == 0 ? 1e6 - 1 // 5X slower - twice
       : strcmp(name, "d27n11") == 0 ? 1e6 - 1 // seems to be giving problems
       : strcmp(name, "d36n04") == 0 ? 1e6 - 1 // 40X slower - at least once
       : strcmp(name, "d23n01") == 0 ? 1e6 - 1 // several times slower - once
       : strcmp(name, "h41n11") == 0 ? 1e4 - 1 // slow Peak node

       : strcmp(name, "a12n18") == 0 ? 1e4 - 1 // to help workaround: ERROR:  One or more process terminated with signal 9
       : strcmp(name, "a03n12") == 0 ? 1e4 - 1 // ditto

       //: strcmp(name, "a03n10") == 0 ? 1e3 - 1 // one of these 5 nodes causes
       //: strcmp(name, "a03n11") == 0 ? 1e3 - 1 // ERROR:  One or more process terminated with signal 9
       //: strcmp(name, "a03n12") == 0 ? 1e3 - 1 //
       //: strcmp(name, "a03n13") == 0 ? 1e3 - 1 //
       //: strcmp(name, "a03n14") == 0 ? 1e3 - 1 //

//       : strcmp(name, "d16n06") == 0 ? 1e6 - 1 // at least 4X slower multiple times
//       : strcmp(name, "d15n03") == 0 ? 1e6 - 1
//       : strcmp(name, "f11n11") == 0 ? 1e6 - 1
//       : strcmp(name, "f13n10") == 0 ? 1e6 - 1
       : -1.;
// b35n16:	62.965595
}

//=============================================================================

typedef struct {
  double time;
  int node;
} PFElt;

int pfelt_cmp(const void* e1, const void* e2) {
  PFElt pf_elt1 = *(PFElt*)e1;
  PFElt pf_elt2 = *(PFElt*)e2;
  return pf_elt1.time < pf_elt2.time ? -1 :
         pf_elt1.time > pf_elt2.time ? 1 : 0;
}

//=============================================================================

void perform_run_preflight_2(int argc, char** argv, MPI_Comm* fast_comm) {

  // Create an env just to extract run options

  GMEnv env_val = GMEnv_null(), *env = &env_val;;
  GMEnv_create(env, MPI_COMM_WORLD, argc, (char**)argv, NULL);

  const int num_rank_requested = GMEnv_num_proc(env);

  int num_rank_avail = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &num_rank_avail);

  if (num_rank_requested == num_rank_avail) {
    MPI_Comm_dup(MPI_COMM_WORLD, fast_comm);
    return;
  }

  const int metric_type = GMEnv_metric_type(env);

  // Identify nodes for which can't open output file (if needed), mark down.

  bool outfile_can_open = true;

  for (int i=1; i<argc; ++i) {
    if (strcmp(argv[i], "--output_file_stub") == 0) {
      if (i < argc-1) {
        FILE* const outfile = gm_metrics_file_open(argv[i+1], env);
        if (outfile) {
          fclose(outfile);
        } else {
          outfile_can_open= false;
        }
        break;
      }
    }
  }

  GMEnv_destroy(env);

  // Initialize communicators.

  //int num_proc = 0;
  //MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //char name[MPI_MAX_PROCESSOR_NAME];
    //int len = 0;
    //MPI_Get_processor_name(name, &len);
    //const size_t len2 = 256;
    //char name2[len2];
    //gethostname(name2, len2);
    //printf("%s %s %i\n", name, name2, rank);

  // Identify node

  const int node_id = get_node_id();

  // Make node-related communicators

  MPI_Comm node_comm;
  MPI_Comm_split(MPI_COMM_WORLD, node_id, rank, &node_comm);

  int rank_in_node = 0;
  MPI_Comm_rank(node_comm, &rank_in_node);

  int ranks_in_node = 0;
  MPI_Comm_size(node_comm, &ranks_in_node);

  int max_ranks_in_node = 0;
  MPI_Allreduce(&ranks_in_node, &max_ranks_in_node, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);

  MPI_Comm rank_in_node_comm;
  MPI_Comm_split(MPI_COMM_WORLD, rank_in_node, rank, &rank_in_node_comm);

  int num_node = 0;
  MPI_Comm_size(rank_in_node_comm, &num_node);
  MPI_Bcast(&num_node, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int node_num = 0;
  MPI_Comm_rank(rank_in_node_comm, &node_num);

  // Run single node case on every node

  const char* options_template =
    metric_type == GM_METRIC_TYPE_CZEK && GM_FP_PRECISION_DOUBLE ?
    "--num_field 25000 --num_vector_local 13000 "
    "--metric_type czekanowski --all2all yes --compute_method GPU "
    "--num_proc_vector %i --num_proc_field 1 "
    "--num_phase 1 --phase_min 0 --phase_max 0 --checksum no --verbosity 0"
    : metric_type == GM_METRIC_TYPE_CZEK ?
    //"--num_field 1 --num_vector_local 2 "
    //"--num_field 560 --num_vector_local 150 "
    //"--num_field 5600 --num_vector_local 1500 "
    //"--num_field 56000 --num_vector_local 5000 "
    "--num_field 56000 --num_vector_local 15000 "
    "--metric_type czekanowski --all2all yes --compute_method GPU "
    "--num_proc_vector %i --num_proc_field 1 "
    "--num_phase 1 --phase_min 0 --phase_max 0 --checksum no --verbosity 0"
    :
    "--num_field 1280000 --num_vector_local 4000 --metric_type ccc --sparse no "
    "--all2all yes --compute_method GPU --num_proc_vector %i "
    "--num_proc_field 1 --num_phase 1 --phase_min 0 --phase_max 0 "
    "--checksum no --verbosity 0";

  char options[1024];
  sprintf(options, options_template, ranks_in_node);

  int num_trial = 1; //FIX 3;
  double max_time = 0.;

  for (int trial=0; trial<num_trial; ++trial) {
    const size_t len = 256;
    char name[len];
    gethostname(name, len);
    const double penalty = !outfile_can_open ? 991010. : bad_node_penalty();
    if (penalty > 0.) {
      max_time = penalty;
      continue;
    }
    GMEnv env_val = GMEnv_null(), *env = &env_val;;
    GMEnv_create(env, node_comm, options, NULL);
    double t1 = GMEnv_get_synced_time(env);
    perform_run(options, node_comm, env);
    double t2 = GMEnv_get_synced_time(env);
    double time = t2 - t1;
    if (!(trial == 0 && num_trial != 1)) {
      max_time = time > max_time ? time : max_time;
    }
    GMEnv_destroy(env);
  }
  if (rank_in_node == 0) {
    const size_t len = 256;
    char name[len];
    gethostname(name, len);
    printf("Warmup run: node %s time %f\n", name, max_time);
  }

  // Collect all timings to node 0

  double* max_times = (double*)malloc(num_node * sizeof(*max_times));

  MPI_Gather(&max_time, 1, MPI_DOUBLE,
             max_times, 1, MPI_DOUBLE, 0, rank_in_node_comm);

  // Sort

  int* node_ranking = (int*)malloc(num_node * sizeof(*node_ranking));

  if (rank == 0) {
    PFElt* pf_elt = (PFElt*)malloc(num_node * sizeof(*pf_elt));
    for (int i=0; i<num_node; ++i) {
      pf_elt[i].time = max_times[i];
      pf_elt[i].node = i;
    }
    qsort(pf_elt, num_node, sizeof(PFElt), *pfelt_cmp);

    for (int i=0; i<num_node; ++i) {
      node_ranking[i] = pf_elt[i].node;
      //printf("%i %i %f\n", i, pf_elt[i].node,  pf_elt[i].time);
    }

    free(pf_elt);
  }

  // Broadcast ranking

  MPI_Bcast(node_ranking, num_node, MPI_INT, 0, MPI_COMM_WORLD);

  int node_ranking_this = 0;
  for (int i=0; i<num_node; ++i) {
    if (node_ranking[i] == node_num) {
      node_ranking_this = i;
    }
  }

  //if (rank_in_node == 0) {
  //  printf("%i %f\n", node_ranking_this, max_time);
  //}

  const int proc_ranking_this = rank_in_node +
    max_ranks_in_node * node_ranking_this;

  // Create world communicator with this ranking

  MPI_Comm_split(MPI_COMM_WORLD, 0, proc_ranking_this, fast_comm);

  // Cleanup

  free(node_ranking);
  free(max_times);
  MPI_Comm_free(&rank_in_node_comm);
  MPI_Comm_free(&node_comm);
}

//=============================================================================

void perform_run_preflight(int argc, char** argv) {

  GMEnv env_val = GMEnv_null(), *env = &env_val;;
  GMEnv_create(env, MPI_COMM_WORLD, argc, (char**)argv, NULL);

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {

    /*---Perform preliminary run on GPU since sometimes first use is slower---*/

    int num_proc = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    const char* options_template_1 =
        GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK ?
          "--num_field 1 --num_vector_local 2 "
          "--metric_type ccc "
          "--num_proc_vector %i --all2all no --num_way 2 "
          "--compute_method GPU --verbosity 0" :
          "--num_field 1 --num_vector_local 2 "
          "--metric_type czekanowski "
          "--num_proc_vector %i --all2all no --num_way 2 "
          "--compute_method GPU --verbosity 0";

    char options1[1024];
    sprintf(options1, options_template_1, num_proc);

    perform_run(options1);
  }

  GMEnv_destroy(env);
}

//=============================================================================
/*---Main---*/

int main(int argc, char** argv) {
  /*---Initialize---*/

  const double t1 = GMEnv_get_time();

  MPI_Init(&argc, &argv);

  bool use_fast_nodes = false;
  for (int i=1; i<argc; ++i) {
    if (strcmp(argv[i], "--fastnodes") == 0) {
      use_fast_nodes = true;
    }
  }

  setbuf(stdout, NULL);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc == 1) {
    if (rank == 0) {
      usage();
    }
    MPI_Finalize();
    return 0;
  }

  //install_handler();

  //const bool use_fast_nodes = false;
  //const bool use_fast_nodes = true;

  if (use_fast_nodes) { 

    MPI_Barrier(MPI_COMM_WORLD);
    const double t2 = GMEnv_get_time();
    if (rank == 0) {
      printf("MPI_Init called, %i seconds.\n", (int)(.5+t2-t1));
    }

    MPI_Comm fast_comm;

    /*---Perform preflight warmup---*/

    perform_run_preflight(argc, argv);

    perform_run_preflight_2(argc, argv, &fast_comm);

    /*---Perform actual run---*/

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
      printf("Commencing run.\n");
    }

    perform_run(argc, (char**)argv, NULL, fast_comm);

    MPI_Comm_free(&fast_comm);

  } else {

    /*---Perform preflight warmup---*/

    perform_run_preflight(argc, argv);

    /*---Perform actual run---*/

    perform_run(argc, (char**)argv, NULL, MPI_COMM_WORLD);

  }

  MPI_Finalize();
  return 0;
}

//-----------------------------------------------------------------------------
