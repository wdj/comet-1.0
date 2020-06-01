//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Top-level function to calculate metrics.
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

#include "stdio.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics_2way.hh"
#include "compute_metrics_3way.hh"
#include "compute_metrics.hh"

//=============================================================================

void GMComputeMetrics_create(
    GMComputeMetrics* this_,
    GMDecompMgr* dm,
    GMEnv* env) {
  GMInsist(this_ && dm && env);

  double time_begin = GMEnv_get_synced_time(env);

  GMComputeMetrics2Way_create(&(this_->compute_metrics_2way), dm, env);
  GMComputeMetrics3Way_create(&(this_->compute_metrics_3way), dm, env);

  double time_end = GMEnv_get_synced_time(env);
  env->time += time_end - time_begin;
}

//-----------------------------------------------------------------------------

void GMComputeMetrics_destroy(
    GMComputeMetrics* this_,
    GMEnv* env) {
  GMInsist(this_ && env);

  double time_begin = GMEnv_get_synced_time(env);

  GMComputeMetrics2Way_destroy(&(this_->compute_metrics_2way), env);
  GMComputeMetrics3Way_destroy(&(this_->compute_metrics_3way), env);

  double time_end = GMEnv_get_synced_time(env);
  env->time += time_end - time_begin;
}

//=============================================================================

void gm_compute_metrics(GMMetrics* metrics, GMVectors* vectors, GMEnv* env) {
  GMInsist(metrics && vectors && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  GMComputeMetrics compute_metrics_value = {0},
                  *compute_metrics = &compute_metrics_value;

  GMComputeMetrics_create(compute_metrics, vectors->dm, env);

  gm_compute_metrics(compute_metrics, metrics, vectors, env);

  GMComputeMetrics_destroy(compute_metrics, env);
}

//=============================================================================

void gm_compute_metrics(GMComputeMetrics* compute_metrics, GMMetrics* metrics,
                        GMVectors* vectors, GMEnv* env) {
  GMInsist(metrics && vectors && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  // Start timer.

  double time_begin = GMEnv_get_synced_time(env);

  // Perform metrics computation.

  if (GMEnv_num_way(env) == 2 && ! GMEnv_all2all(env)) {

    gm_compute_metrics_2way_notall2all(
      &(compute_metrics->compute_metrics_2way), metrics, vectors, env);

  } else if (GMEnv_num_way(env) == 2 && GMEnv_all2all(env)) {

    gm_compute_metrics_2way_all2all(
      &(compute_metrics->compute_metrics_2way), metrics, vectors, env);

  } else if (GMEnv_num_way(env) == 3 && ! GMEnv_all2all(env)) {

    gm_compute_metrics_3way_notall2all(
      &(compute_metrics->compute_metrics_3way), metrics, vectors, env);

  } else if (GMEnv_num_way(env) == 3 && GMEnv_all2all(env)) {

    gm_compute_metrics_3way_all2all(
      &(compute_metrics->compute_metrics_3way), metrics, vectors, env);

  } else {

    GMInsistInterface(env, false &&
                      "This num_way/all2all combination unimplemented.");

  }

  // Stop timer.

  double time_end = GMEnv_get_synced_time(env);
  env->time += time_end - time_begin;

  // Check computed element count.

  GMInsist(metrics->num_elts_local == metrics->num_elts_local_computed &&
           "Failure to compute all requested metrics.");

  // Compute global counts of compares and operations.

  double num_elts_local = metrics->num_elts_local;
  double num_elts = 0;

  int mpi_code = MPI_Allreduce(&num_elts_local, &num_elts, 1,
                           MPI_DOUBLE, MPI_SUM, GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Allreduce.");

  env->compares += metrics->num_field_active * num_elts *
                   metrics->data_type_num_values;

  env->eltcompares += metrics->num_field_active * num_elts;

  env->veccompares += num_elts;

  mpi_code = MPI_Allreduce(&env->ops_local, &env->ops, 1, MPI_DOUBLE, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Allreduce.");

  // Compute global CPU, GPU memory high water marks.

  const size_t cpu_mem_max_local = env->cpu_mem_max;
  mpi_code = MPI_Allreduce(&cpu_mem_max_local, &env->cpu_mem_max, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Allreduce.");

  const size_t gpu_mem_max_local = env->gpu_mem_max;
  mpi_code = MPI_Allreduce(&gpu_mem_max_local, &env->gpu_mem_max, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Allreduce.");
}

//-----------------------------------------------------------------------------
