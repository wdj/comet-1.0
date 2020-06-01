//-----------------------------------------------------------------------------
/*!
 * \file   driver.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions.
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
#include "float.h"
#include "errno.h"

#include "unistd.h"

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "checksum.hh"
#include "compute_metrics.hh"

#include "test_problems.hh"
#include "input_output.hh"
#include "driver.hh"

//=============================================================================
/*---Parse remaining unprocessed arguments---*/

void finish_parsing(int argc, char** argv, DriverOptions* do_, GMEnv* env) {
  errno = 0;
  int i = 0;
  for (i = 1; i < argc; ++i) {
    /*----------*/
    if (strcmp(argv[i], "--num_field") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_field.");
      const long num_field = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && num_field >= 0
                    && "Invalid setting for num_field.");
      do_->num_field_active = num_field;
      do_->num_field_active_initialized = true;
      do_->num_field_local_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--num_field_local") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_field_local.");
      const long num_field_local = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno && num_field_local >= 0 &&
                    (long)(int)num_field_local == num_field_local &&
                    "Invalid setting for num_field_local.");
      do_->num_field_local = num_field_local;
      do_->num_field_local_initialized = true;
      do_->num_field_active_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--num_vector") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_vector.");
      const long num_vector = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && num_vector >= 0
                    && "Invalid setting for num_vector.");
      do_->num_vector_active = num_vector;
      do_->num_vector_active_initialized = true;
      do_->num_vector_local_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--num_vector_local") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_vector_local.");
      const long num_vector_local = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno && num_vector_local >= 0 &&
                    (long)(int)num_vector_local == num_vector_local &&
                    "Invalid setting for num_vector_local.");
      do_->num_vector_local = num_vector_local;
      do_->num_vector_local_initialized = true;
      do_->num_vector_active_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--verbosity") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for verbosity.");
      const long verbosity = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno && verbosity >= 0 &&
                    "Invalid setting for verbosity.");
      do_->verbosity = verbosity;
      /*--------------------*/
    } else if (strcmp(argv[i], "--checksum") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for checksum.");
      if (strcmp(argv[i], "yes") == 0) {
        do_->checksum = true;
      } else if (strcmp(argv[i], "no") == 0) {
        do_->checksum = false;
      } else {
        GMInsistInterface(env, false && "Invalid setting for checksum.");
      }
    /*----------*/
    } else if (strcmp(argv[i], "--num_stage") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_stage.");
      const long num_stage = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && num_stage >= 1
                    && (long)(int)num_stage == num_stage
                    && "Invalid setting for num_stage.");
      env->num_stage = num_stage;
      do_->stage_min_0based = 0;
      do_->stage_max_0based = env->num_stage - 1;
    /*----------*/
    } else if (strcmp(argv[i], "--stage_min") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for stage_min.");
      const long stage_min_0based = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && stage_min_0based >= 0
                    && (long)(int)stage_min_0based == stage_min_0based
                    && "Invalid setting for stage_min.");
      do_->stage_min_0based = stage_min_0based;
    /*----------*/
    } else if (strcmp(argv[i], "--stage_max") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for stage_max.");
      const long stage_max_0based = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && stage_max_0based < env->num_stage
                    && (long)(int)stage_max_0based == stage_max_0based
                    && "Invalid setting for stage_max.");
      do_->stage_max_0based = stage_max_0based;
    /*----------*/
    } else if (strcmp(argv[i], "--num_phase") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_phase.");
      const long num_phase = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && num_phase >= 1
                    && (long)(int)num_phase == num_phase
                    && "Invalid setting for num_phase.");
      env->num_phase = num_phase;
      do_->phase_min_0based = 0;
      do_->phase_max_0based = env->num_phase - 1;
    /*----------*/
    } else if (strcmp(argv[i], "--phase_min") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for phase_min.");
      const long phase_min_0based = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && phase_min_0based >= 0
                    && (long)(int)phase_min_0based == phase_min_0based
                    && "Invalid setting for phase_min.");
      do_->phase_min_0based = phase_min_0based;
    /*----------*/
    } else if (strcmp(argv[i], "--phase_max") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for phase_max.");
      const long phase_max_0based = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && phase_max_0based < env->num_phase
                    && (long)(int)phase_max_0based == phase_max_0based
                    && "Invalid setting for phase_max.");
      do_->phase_max_0based = phase_max_0based;
    /*----------*/
    } else if (strcmp(argv[i], "--input_file") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for input_file.");
      do_->input_file_path = argv[i];
    /*----------*/
    } else if (strcmp(argv[i], "--output_file_stub") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for output_file_stub.");
      do_->metrics_file_path_stub = argv[i];
      /*--------------------*/
    } else if (strcmp(argv[i], "--problem_type") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for problem_type.");
      if (strcmp(argv[i], "random") == 0) {
        do_->problem_type = GM_PROBLEM_TYPE_RANDOM;
        //GMEnv_set_compute_method(env, GM_PROBLEM_TYPE_RANDOM);
      } else if (strcmp(argv[i], "analytic") == 0) {
        do_->problem_type = GM_PROBLEM_TYPE_ANALYTIC;
        //GMEnv_set_compute_method(env, GM_PROBLEM_TYPE_ANALYTIC);
      } else {
        GMInsistInterface(env, false && "Invalid setting for problem_type.");
      }
    /*----------*/
    } else if (strcmp(argv[i], "--threshold") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for threshold.");
      errno = 0;
      const double threshold = strtod(argv[i], NULL);
      GMInsistInterface(env, 0 == errno && "Invalid setting for threshold.");
      do_->threshold = threshold;
     /*----------*/
    } else if (strcmp(argv[i], "--metric_type") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_way") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--all2all") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--ccc_multiplier") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--duo_multiplier") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--ccc_param") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--sparse") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--fastnodes") == 0) {
      /*---optionally processed by caller---*/
    } else if (strcmp(argv[i], "--tc") == 0) {
      ++i; /*---optionally processed by caller---*/
    } else if (strcmp(argv[i], "--num_tc_steps") == 0) {
      ++i; /*---optionally processed by caller---*/
    } else {
    /*----------*/
      if (GMEnv_proc_num(env) == 0) {
        fprintf(stderr, "Invalid argument \"%s\". ", argv[i]);
      }
      GMInsistInterface(env, false && "Error: argument not recognized.");
    /*----------*/
    } /*---if/else---*/

  } /*---for i---*/

  GMInsistInterface(env, (do_->num_field_local_initialized ||
                do_->num_field_active_initialized)
                && "Error: must set num_field_local or num_field.");
  GMInsistInterface(env, (do_->num_vector_local_initialized ||
                do_->num_vector_active_initialized)
                && "Error: must set num_vector_local or num_vector.");
}

//-----------------------------------------------------------------------------

void set_vectors(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMInsist(vectors && do_ && env);

  if (do_->input_file_path != NULL) {
    set_vectors_from_file(vectors, do_, env);
  } else {
    set_vectors_synthetic(vectors, do_->problem_type, do_->verbosity, env);
  }
}
//=============================================================================
/*---Perform a single metrics computation run---*/

void perform_run(const char* const options, MPI_Comm base_comm, GMEnv* env) {
  GMInsist(options);

  CoMet::Checksum cksum;

  perform_run(cksum, options, base_comm, env);
}

//-----------------------------------------------------------------------------

void perform_run(CoMet::Checksum& cksum, const char* const options,
                 MPI_Comm base_comm, GMEnv* env) {
  GMInsist(options);

  /*---Convert options string to args---*/

  size_t len = strlen(options);
  char argstring[len+1];
  char* argv[len+1];
  int argc = 0;
  strcpy(argstring, options);
  gm_create_args(argstring, &argc, argv);

  return perform_run(cksum, argc, argv, options, base_comm, env);
}

//-----------------------------------------------------------------------------

void perform_run(int argc, char** argv, const char* const description,
                            MPI_Comm base_comm, GMEnv* env) {

  CoMet::Checksum cksum;

  perform_run(cksum, argc, argv, description, base_comm, env);
}

//-----------------------------------------------------------------------------

void perform_run(CoMet::Checksum& cksum_result, int argc, char** argv,
                 const char* const description,
                 MPI_Comm base_comm, GMEnv* env) {

  /*---Initialize environment---*/

  bool create_env = ! env;

  GMEnv env_local = GMEnv_null();

  if (create_env) {
    env = &env_local;
    GMEnv_create(env, base_comm, argc, argv, description);
  }

  if (! GMEnv_is_proc_active(env)) {
    if (create_env) {
      GMEnv_destroy(env);
    }
    return;
  }

  double total_time_beg = GMEnv_get_synced_time(env);

  /*---Parse remaining unprocessed arguments---*/

  DriverOptions do_ = {0};
  do_.num_field_local_initialized = false;
  do_.num_field_active_initialized = false;
  do_.num_vector_local_initialized = false;
  do_.num_vector_active_initialized = false;
  do_.verbosity = 1;
  do_.stage_min_0based = 0;
  do_.stage_max_0based = env->num_stage - 1;
  do_.phase_min_0based = 0;
  do_.phase_max_0based = env->num_phase - 1;
  do_.input_file_path = NULL;
  do_.metrics_file_path_stub = NULL;
  //do_.problem_type = GM_PROBLEM_TYPE_RANDOM;
  //do_.problem_type = GM_PROBLEM_TYPE_ANALYTIC;
  do_.problem_type = problem_type_default();
  do_.threshold = -1.;
  do_.checksum = true;
  do_.num_incorrect = 0;
  do_.max_incorrect_diff = 0.;

  finish_parsing(argc, argv, &do_, env);

  /*---Set up parallel deomp for vectors, metrics---*/

  double vctime = 0;
  double time_beg = GMEnv_get_synced_time(env);
  GMDecompMgr dm_value = GMDecompMgr_null(), *dm = &dm_value;
  GMDecompMgr_create(dm,
    do_.num_field_local_initialized,
    do_.num_vector_local_initialized,
    do_.num_field_local_initialized ? do_.num_field_local
                                    : do_.num_field_active,
    do_.num_vector_local_initialized ? do_.num_vector_local
                                     : do_.num_vector_active,
    GMEnv_data_type_vectors(env), env);
  double time_end = GMEnv_get_synced_time(env);
  vctime += time_end - time_beg;

//TODO: possibly replace this with stuff from dm
  if (do_.num_vector_local_initialized) {
    do_.num_vector = do_.num_vector_local *
      (size_t)GMEnv_num_proc_vector_i(env);
    do_.num_vector_active = do_.num_vector;
  } else {
    /*---Pad up so that every proc has same number of vectors---*/
    do_.num_vector_local = gm_num_vector_local_required(
      gm_ceil_i8(do_.num_vector_active, GMEnv_num_proc_vector_i(env)), env);
    do_.num_vector = do_.num_vector_local *
      (size_t)GMEnv_num_proc_vector_i(env);
  }

  if (do_.num_field_local_initialized) {
    do_.num_field = do_.num_field_local * (size_t) GMEnv_num_proc_field(env);
    do_.num_field_active = do_.num_field;
  } else {
    /*---Pad up so that every proc has same number of fields---*/
    do_.num_field_local = gm_ceil_i8(
        do_.num_field_active, GMEnv_num_proc_field(env));
    do_.num_field = do_.num_field_local * (size_t) GMEnv_num_proc_field(env);
  }

//printf("%i %i %i %i\n", env->proc_num_base_, env->proc_num_, env->proc_num_repl_, env->proc_num_vector_i_);

  const bool do_print = GMEnv_is_proc_active(env) &&
     GMEnv_proc_num(env) == 0 && do_.verbosity > 0;

  /*---Allocate vectors---*/

  time_beg = GMEnv_get_synced_time(env);
  GMVectors vectors_value = GMVectors_null(), *vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), dm, env);
  time_end = GMEnv_get_synced_time(env);
  vctime += time_end - time_beg;

  /*---Set vectors---*/

  double intime = 0;
  time_beg = GMEnv_get_synced_time(env);

//double t1 = GMEnv_get_time(env);

  set_vectors(vectors, &do_, env);

//double t2 = GMEnv_get_time(env);
//printf("TIME %f %i\n", t2-t1, GMEnv_is_proc_active(env));

//  int rank = 0;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  printf("RANK %i\n", rank);

  time_end = GMEnv_get_synced_time(env);
  intime += time_end - time_beg;

  /*---More initializations---*/

  CoMet::Checksum cksum(do_.checksum);
  CoMet::Checksum cksum_local(do_.checksum);

  double outtime = 0;
  double mctime = 0;
  double cktime = 0;

  size_t num_elts_local_computed = 0;
  size_t num_local_written = 0;

  /*---Open output files---*/

  {
  time_beg = GMEnv_get_synced_time(env);
  MetricsFile metric_file(&do_, env);
  time_end = GMEnv_get_synced_time(env);
  outtime += time_end - time_beg;

  {
  GMMetricsMem metrics_mem(env);

  GMComputeMetrics compute_metrics_value = {0},
                  *compute_metrics = &compute_metrics_value;

  GMComputeMetrics_create(compute_metrics, dm, env);

  /*---Loops over phases, stages---*/

  for (env->phase_num=do_.phase_min_0based;
       env->phase_num<=do_.phase_max_0based; ++env->phase_num) {

    for (env->stage_num=do_.stage_min_0based;
         env->stage_num<=do_.stage_max_0based; ++env->stage_num) {

      /*---Set up metrics container for results---*/

      time_beg = GMEnv_get_synced_time(env);
      GMMetrics metrics_value = GMMetrics_null(), *metrics = &metrics_value;
      GMMetrics_create(metrics, GMEnv_data_type_metrics(env), dm,
                       &metrics_mem, env);
      time_end = GMEnv_get_synced_time(env);
      mctime += time_end - time_beg;

      /*---Calculate metrics---*/

      gm_compute_metrics(compute_metrics, metrics, vectors, env);

      num_elts_local_computed += metrics->num_elts_local_computed;

      /*---Output results---*/

      time_beg = GMEnv_get_synced_time(env);
      metric_file.write(metrics, env);
      time_end = GMEnv_get_synced_time(env);
      outtime += time_end - time_beg;

      /*---Check correctness---*/

      if (do_.checksum) {
        time_beg = GMEnv_get_synced_time(env);
        check_metrics(metrics, &do_, env);
        time_end = GMEnv_get_synced_time(env);
        cktime += time_end - time_beg;
      }

      /*---Compute checksum---*/

      if (do_.checksum) {
        time_beg = GMEnv_get_synced_time(env);
        CoMet::Checksum::compute(cksum, cksum_local, *metrics, *env);
        time_end = GMEnv_get_synced_time(env);
        cktime += time_end - time_beg;
      }
      time_beg = GMEnv_get_synced_time(env);
      GMMetrics_destroy(metrics, env);
      time_end = GMEnv_get_synced_time(env);
      mctime += time_end - time_beg;

      if (do_print) {
        if (env->num_phase > 1 && env->num_stage > 1) {
          printf("Completed phase %i stage %i\n",
                 env->phase_num, env->stage_num);
        } else if (env->num_phase > 1) {
          printf("Completed phase %i\n",
                 env->phase_num);
        } else if (env->num_stage > 1) {
          printf("Completed stage %i\n",
                 env->stage_num);
        }
      }

    }

  } /*---End loops over phases, stages---*/

  GMComputeMetrics_destroy(compute_metrics, env);

  /*---Finalize metrics mem---*/

  time_beg = GMEnv_get_synced_time(env);
  }
  time_end = GMEnv_get_synced_time(env);
  mctime += time_end - time_beg;
  /*---Close output files---*/

  num_local_written += metric_file.get_num_written();
  time_beg = GMEnv_get_synced_time(env);
  }
  time_end = GMEnv_get_synced_time(env);
  outtime += time_end - time_beg;

  /*---Deallocate vectors---*/

  time_beg = GMEnv_get_synced_time(env);
  GMVectors_destroy(vectors, env);
  GMDecompMgr_destroy(dm, env);
  time_end = GMEnv_get_synced_time(env);
  vctime += time_end - time_beg;

  /*---Perform some checks---*/

  GMInsist(env->cpu_mem == 0);
  GMInsist(env->gpu_mem == 0);

  size_t num_written = 0;
  if (GMEnv_is_proc_active(env)) {
    int mpi_code = 0;
    size_t num_elts_computed = 0;
    mpi_code = MPI_Allreduce(&num_elts_local_computed, &num_elts_computed, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                             GMEnv_mpi_comm_repl_vector(env));
    GMInsist(mpi_code == MPI_SUCCESS);

    mpi_code = MPI_Allreduce(&num_local_written, &num_written, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                             GMEnv_mpi_comm_repl_vector(env));
    GMInsist(mpi_code == MPI_SUCCESS);

    if (GMEnv_num_way(env) == GM_NUM_WAY_2 && GMEnv_all2all(env) &&
        do_.phase_min_0based==0 && do_.phase_max_0based==env->num_phase - 1) {
      GMInsist(num_elts_computed == (do_.num_vector) * (size_t)
                                          (do_.num_vector - 1) / 2);
    }

    if (GMEnv_num_way(env) == GM_NUM_WAY_3 && GMEnv_all2all(env) &&
        do_.phase_min_0based==0 && do_.phase_max_0based==env->num_phase - 1 &&
        do_.stage_min_0based==0 && do_.stage_max_0based==env->num_stage - 1) {
      GMInsist(num_elts_computed == (do_.num_vector) * (size_t)
                                          (do_.num_vector - 1) * (size_t)
                                          (do_.num_vector - 2) / 6);
    }
  }

  double total_time_end = GMEnv_get_synced_time(env);

  /*---Output run information---*/

  if (do_print) {
    //-----
    if (do_.checksum) {
      printf("metrics checksum ");
      //GMChecksum_print(cksum, env);
      cksum.print(*env);
      printf(" ");
    }
    //-----
    printf("ctime %.6f", env->time);
    //-----
    printf(" ops %e", env->ops);
    if (env->time > 0) {
      printf(" ops_rate %e", env->ops / env->time);
      printf(" ops_rate/proc %e", env->ops / (env->time*GMEnv_num_proc(env)) );
    }
    //-----
    printf(" vcmp %e", env->veccompares);
    if (NULL != do_.metrics_file_path_stub) {
      printf(" vcmpout %e", (double)num_written);
    }
    //-----
    printf(" cmp %e", env->compares);
    printf(" ecmp %e", env->eltcompares);
    if (env->time > 0) {
      printf(" ecmp_rate %e", env->eltcompares / env->time);
      printf(" ecmp_rate/proc %e", env->eltcompares / (env->time*GMEnv_num_proc(env)) );
    }
    //-----
    printf(" vctime %.6f", vctime);
    printf(" mctime %.6f", mctime);
    if (do_.checksum) {
      printf(" cktime %.6f", cktime);
    }
    //if (NULL != do_.input_file_path) {
    printf(" intime %.6f", intime);
    //}
    //if (NULL != do_.metrics_file_path_stub) {
    printf(" outtime %.6f", outtime);
    //}
    //-----
    printf(" cpumem %e", (double)env->cpu_mem_max);
    printf(" gpumem %e", (double)env->gpu_mem_max);
    //-----
    printf(" tottime %.6f", total_time_end - total_time_beg);
    //-----
    printf("\n");
  }

  // Output a local checksum, for testing purposes.

  if (false) {
    // One more sync before checking num_correct, to allow flush of output.
    GMEnv_get_synced_time(env);
    if (do_.checksum && GMEnv_is_proc_active(env) && do_.verbosity > 0) {
      printf("local checksum: ");
      cksum_local.print(*env);
      printf("\n");
    }
  }
  GMEnv_get_synced_time(env);

  // Validation: check for any wrong answers.

  if (do_.num_incorrect) {
    const size_t hnlen = 256;
    char hn[hnlen];
    gethostname(hn, hnlen);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Error: incorrect results found.  num_incorrect  %zu  "
           "max_incorrect_diff  %e  hostname  %s  rank  %i\n",
           do_.num_incorrect, do_.max_incorrect_diff, hn, rank);
  }

  GMInsist(do_.num_incorrect == 0);

  /*---Finalize---*/

  if (create_env) {
    GMEnv_destroy(env);
  }

  cksum_result.copy(cksum);
}

//=============================================================================
