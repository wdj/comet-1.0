//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block.
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

#include "string.h"

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "comm_xfer_utils.hh"
#include "compute_metrics_3way_block_gpu.hh"
#include "compute_metrics_3way_block_nongpu.hh"
#include "compute_metrics_3way_block.hh"

//=============================================================================

void GMComputeNumerators3Way_create(GMComputeNumerators3Way* this_,
                                    int nvl, int npvfl, GMEnv* env) {
  GMInsist(this_ && env);
  GMInsist(nvl >= 0 && npvfl >= 0);

  this_->matM_ij_buf = GMMirroredBuf_null();
  this_->matM_jk_buf = GMMirroredBuf_null();
  this_->matM_kik_buf = GMMirroredBuf_null();

  for (int i=0; i<2; ++i) {
    this_->tmp_buf[i] = GMMirroredBuf_null();
    this_->matX_buf[i] = GMMirroredBuf_null();
    this_->matB_buf[i] = GMMirroredBuf_null();
  }

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    for (int i=0; i<2; ++i) {
      if (env->do_reduce) {
        GMMirroredBuf_create(&(this_->tmp_buf[i]), nvl, nvl, env);
      }
      GMMirroredBuf_create(&(this_->matX_buf[i]), npvfl, nvl, env);
      GMMirroredBuf_create(&(this_->matB_buf[i]), nvl, nvl, env);
    }
    if (env->need_2way) {
      GMMirroredBuf_create(&(this_->matM_ij_buf), nvl, nvl, env);
      GMMirroredBuf_create(&(this_->matM_jk_buf), nvl, nvl, env);
      GMMirroredBuf_create(&(this_->matM_kik_buf), nvl, nvl, env);
    }
  }
}

//=============================================================================

void GMComputeNumerators3Way_destroy(GMComputeNumerators3Way* this_,
                                     GMEnv* env) {
  GMInsist(this_ && env);

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    for (int i=0; i<2; ++i) {
      if (env->do_reduce) {
        GMMirroredBuf_destroy(&this_->tmp_buf[i], env);
      }
      GMMirroredBuf_destroy(&this_->matX_buf[i], env);
      GMMirroredBuf_destroy(&this_->matB_buf[i], env);
    }
    if (env->need_2way) {
      GMMirroredBuf_destroy(&this_->matM_ij_buf, env);
      GMMirroredBuf_destroy(&this_->matM_jk_buf, env);
      GMMirroredBuf_destroy(&this_->matM_kik_buf, env);
    }
  }
}

//=============================================================================
/*---Start calculation of numerators, 3-way generic---*/

void GMComputeNumerators3Way_start(
    GMComputeNumerators3Way* this_,
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredBuf* vectors_i_buf,
    GMMirroredBuf* vectors_j_buf,
    GMMirroredBuf* vectors_k_buf,
    int j_block,
    int k_block,
    const GMVectorSums* vector_sums_i,
    const GMVectorSums* vector_sums_j,
    const GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env) {
  GMInsist(this_ && metrics && env);
  GMInsist(vectors_i && vectors_j && vectors_k);
  GMInsist(vectors_i_buf && vectors_j_buf && vectors_k_buf);
  GMInsist(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMInsist(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMInsist(! (GMEnv_proc_num_vector_i(env) == j_block &&
              GMEnv_proc_num_vector_i(env) != k_block));
  GMInsist(! (GMEnv_proc_num_vector_i(env) == k_block &&
              GMEnv_proc_num_vector_i(env) != j_block));
  GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMInsist(vector_sums_i && vector_sums_j && vector_sums_k);

  /*----------------------------------------*/
  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/
    gm_compute_3way_nums_gpu_start_(this_,
                                    vectors_i, vectors_j, vectors_k,
                                    metrics, vectors_i_buf, vectors_j_buf,
                                    vectors_k_buf, j_block, k_block,
                                    vector_sums_i, vector_sums_j,
                                    vector_sums_k,
                                    section_step, env);
    /*----------------------------------------*/
  } else /*---(GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU)---*/ {
    /*----------------------------------------*/
    switch (GMEnv_metric_type(env)) {
      case GM_METRIC_TYPE_CZEK: {
        gm_compute_3way_nums_nongpu_czek_start_(this_,
            vectors_i, vectors_j, vectors_k, metrics, vectors_i_buf,
            vectors_j_buf, vectors_k_buf, j_block, k_block,
            vector_sums_i, vector_sums_j, vector_sums_k,
            section_step, env);
      } break;
      case GM_METRIC_TYPE_CCC: {
        gm_compute_3way_nums_nongpu_ccc_start_(this_,
            vectors_i, vectors_j, vectors_k, metrics, vectors_i_buf,
            vectors_j_buf, vectors_k_buf, j_block, k_block,
            vector_sums_i, vector_sums_j, vector_sums_k,
            section_step, env);
      } break;
      default:
        GMInsistInterface(env, false && "Selected metric_type unimplemented.");
    } /*---case---*/
    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================

//-----------------------------------------------------------------------------
