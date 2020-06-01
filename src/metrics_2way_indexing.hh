//-----------------------------------------------------------------------------
/*!
 * \file   metrics_2way_indexing.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 2-way, indexing.
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

#ifndef _gm_metrics_2way_indexing_hh_
#define _gm_metrics_2way_indexing_hh_

//=============================================================================
/*---Helper functions for 2-way case---*/

static int gm_bdiag_computed_max_allphase(GMEnv* env) {
  GMAssert(env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));

  /*---Max number of blocks of any block row computed on all phases---*/
  const int num_block = GMEnv_num_block_vector(env);
  return 1 + num_block / 2;
}

//-----------------------------------------------------------------------------

static int gm_bdiag_computed_min(GMEnv* env) {
  GMAssert(env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));

  // First block diag computed for this phase (min across vec procs)
  const int max_rectangle_width = gm_bdiag_computed_max_allphase(env);
  return (max_rectangle_width*env->phase_num) / env->num_phase;
}

//-----------------------------------------------------------------------------

static int gm_bdiag_computed_max(GMEnv* env) {
  GMAssert(env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));

  // 1 + last block diag computed for this phase (max across vec procs)
  const int max_rectangle_width = gm_bdiag_computed_max_allphase(env);
  return (max_rectangle_width*(env->phase_num+1)) / env->num_phase;
}

//-----------------------------------------------------------------------------

static int gm_block_computed_this_row_min(GMEnv* env) {
  GMAssert(env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));

  return gm_bdiag_computed_min(env);
}

//-----------------------------------------------------------------------------

static int gm_block_computed_this_row_max(GMEnv* env) {
  GMAssert(env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));

  const int num_block = GMEnv_num_block_vector(env);
  const int i_block = GMEnv_proc_num_vector_i(env);

  const bool is_row_short_by_1 = num_block % 2 == 0 && 2*i_block >= num_block;
  const bool is_last_phase = env->phase_num == env->num_phase - 1;

  const int diag_max = gm_bdiag_computed_max(env);

  // 1 + last block diag computed for this phase, all repl procs (this vec proc)
  const int n = is_last_phase && is_row_short_by_1 ? diag_max - 1 : diag_max;
  GMAssert(n >= 0);
  GMAssert(n <= num_block);
  return n;
}

//-----------------------------------------------------------------------------

static int gm_blocks_computed_this_row(GMEnv* env) {
  GMAssert(env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));

  // num block diags computed for this phase, all repl procs (this vec proc)
  const int n = gm_block_computed_this_row_max(env) -
                gm_block_computed_this_row_min(env);
  GMAssert(n >= 0);
  GMAssert(n <= GMEnv_num_block_vector(env));
  return n;
}

//=============================================================================
/*---Accessors: indexing: (contig) index from coord, 2-way---*/

static size_t gm_triang_(int i) {
  return (i * (size_t)(i-1)) >> 1;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_index_from_coord_2(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_proc_num_repl(env) == 0);

  size_t index = gm_triang_(j) + i;
  GMAssert(i + metrics->num_vector_local *
               (size_t)GMEnv_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] % metrics->num_vector);
  GMAssert(j + metrics->num_vector_local *
               (size_t)GMEnv_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] / metrics->num_vector);
  return index;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper2way_maindiag_block_(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int j_block,
                                                   GMEnv* env) {
  GMAssert(j_block == GMEnv_proc_num_vector_i(env));

  return gm_triang_(j) + i;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper2way_offdiag_block_(GMMetrics* metrics,
                                                  int i,
                                                  int j,
                                                  int j_block,
                                                  GMEnv* env) {
  GMAssert(j_block != GMEnv_proc_num_vector_i(env));

  const int num_block = GMEnv_num_block_vector(env);

  const int num_proc_r = GMEnv_num_proc_repl(env);

  const int block_min = metrics->block_min;

  /* clang-format off */
  return metrics->index_offset_0_ +
      i + metrics->num_vector_local * (size_t)(
      j + metrics->num_vector_local * (
      ((j_block - block_min + num_block) % num_block) / num_proc_r ));
  /* clang-format on */
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_index_from_coord_all2all_2(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int j_block,
                                                   GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
//  GMAssert(GMEnv_proc_num_repl(env) == 0 ||
//           j_block != GMEnv_proc_num_vector(env)); // DEFUNCT
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  const int i_block = GMEnv_proc_num_vector_i(env);

  size_t index = j_block == i_block
           ? GMMetrics_helper2way_maindiag_block_(metrics, i, j, j_block, env)
           : GMMetrics_helper2way_offdiag_block_(metrics, i, j, j_block, env);

  GMAssert(index >= 0 && index < metrics->num_elts_local);
  return index;
}

//=============================================================================
//=============================================================================
/*---Accessors: indexing: global coord from (contig) index: 2-way---*/

static int GMMetrics_coord0_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);


  const size_t i64 = metrics->coords_global_from_index[index] %
                     metrics->num_vector;
  const int i = (int)i64;
  GMAssert((size_t)i == i64);

  return i;
}

//-----------------------------------------------------------------------------

static int GMMetrics_coord1_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);

  const size_t j64 = metrics->coords_global_from_index[index] /
                     metrics->num_vector;
  const int j = (int)j64;
  GMAssert((size_t)j == j64);

  return j;
}

//=============================================================================

#endif // _gm_metrics_2way_indexing_hh_

//-----------------------------------------------------------------------------
