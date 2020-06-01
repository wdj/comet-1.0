//-----------------------------------------------------------------------------
/*!
 * \file   metrics_3way_indexing.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 3-way, indexing.
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

#ifndef _gm_metrics_3way_indexing_hh_
#define _gm_metrics_3way_indexing_hh_

#include "cstdint"

//=============================================================================
/*---Helper functions for 3-way case---*/

//-----------------------------------------------------------------------------
/*---NOTE: the following does not specialize based on part1/2/3---*/

static int gm_num_section_steps(const GMEnv* const env, int part_num) {
  GMAssert(env && part_num >= 1 && part_num <= 3);
  // Number of section steps to be executed for a given block.

  const bool collapse = ! GMEnv_all2all(env) || GMEnv_num_proc_repl(env) == 1;
  return collapse || part_num == 3 ? 1 : 6;
}

//-----------------------------------------------------------------------------

static int gm_num_sections(const GMEnv* const env, int part_num) {
  GMAssert(env && part_num >= 1 && part_num <= 3);
  // Number of sections the block is divided into.

  return part_num == 3 ? 6 : gm_num_section_steps(env, part_num);
}

//-----------------------------------------------------------------------------

static int gm_num_section_blocks(const GMEnv* const env) {
  GMAssert(env);
  // Total section steps across all blocks, phases.

  const int npv = GMEnv_num_proc_vector_i(env);

  const bool collapse = ! GMEnv_all2all(env) || GMEnv_num_proc_repl(env) == 1;
  return collapse ? npv*npv - 2*npv + 2 : (npv+1) * (npv+2);
}

//-----------------------------------------------------------------------------

static int gm_section_block_phase_min(const GMEnv* const env) {
  GMAssert(env);

  return (gm_num_section_blocks(env)*env->phase_num) / env->num_phase;
}

//-----------------------------------------------------------------------------

static int gm_section_block_phase_max(const GMEnv* const env) {
  GMAssert(env);

  return (gm_num_section_blocks(env)*(env->phase_num+1)) / env->num_phase;
}


//-----------------------------------------------------------------------------

static bool gm_is_section_block_in_phase(const GMEnv* const env,
                                        int section_block) {
  GMAssert(env);
  GMAssert(section_block >= 0 && section_block < gm_num_section_blocks(env));

  return section_block >= gm_section_block_phase_min(env) &&
         section_block < gm_section_block_phase_max(env);
}

//-----------------------------------------------------------------------------

static bool gm_is_part1(int i_block, int j_block, int k_block) {
  return i_block == j_block && j_block == k_block;
}

//-----------------------------------------------------------------------------

static bool gm_is_part3(int i_block, int j_block, int k_block) {
  return i_block != j_block && j_block != k_block && i_block != k_block;
}

//-----------------------------------------------------------------------------

static int gm_section_axis_part3(int i_block, int j_block, int k_block) {
  /*---NOTE: this could possibly be implemented somewhat more efficiently---*/
  /* clang-format off */
  return i_block < j_block && i_block < k_block ? 0 : /*---i axis---*/
         j_block < i_block && j_block < k_block ? 1 : /*---j axis---*/
                                                  2;  /*---k axis---*/
  /* clang-format on */
}

//-----------------------------------------------------------------------------

static int gm_section_num_part3(int i_block, int j_block, int k_block) {
  /*---NOTE: this could possibly be implemented somewhat more efficiently---*/
  /* clang-format off */
  return i_block < k_block && k_block < j_block ?    0 :
         i_block < j_block && j_block < k_block ?    1 :
         j_block < i_block && i_block < k_block ?    2 :
         j_block < k_block && k_block < i_block ?    3 :
         k_block < j_block && j_block < i_block ?    4 :
       /*k_block < i_block && i_block < j_block ? */ 5;
  /* clang-format on */
}

//-----------------------------------------------------------------------------

static int gm_J_lo(int section_num, int nvl, int part_num, GMEnv* env) {
  GMAssert(env);
  GMAssert(section_num >= 0 && section_num < 6);
  GMAssert(nvl >= 0);
  GMAssert(part_num >= 1 && part_num <= 3);
  const int num_sections = gm_num_sections(env, part_num);
  GMAssert(section_num >= 0 && section_num <= num_sections);
  GMAssert(env->num_stage > 0);
  GMAssert(env->stage_num >= 0 && env->stage_num < env->num_stage);

  const int result = ((env->stage_num + env->num_stage * section_num)*nvl) /
                     (num_sections * env->num_stage);

  return result;
}

//-----------------------------------------------------------------------------

static int gm_J_hi(int section_num, int nvl, int part_num, GMEnv* env) {
  GMAssert(env);
  GMAssert(section_num >= 0 && section_num < 6);
  GMAssert(nvl >= 0);
  GMAssert(part_num >= 1 && part_num <= 3);
  const int num_sections = gm_num_sections(env, part_num);
  GMAssert(section_num >= 0 && section_num <= num_sections);
  GMAssert(env->num_stage > 0);
  GMAssert(env->stage_num >= 0 && env->stage_num < env->num_stage);

  const int result = ((env->stage_num + 1 + env->num_stage * section_num)*nvl) /
                     (num_sections * env->num_stage);

  return result;
}

//=============================================================================
/*---GMSectionInfo---*/

typedef struct {
  bool is_part1;
  bool is_part2;
  bool is_part3;
  bool sax0;
  bool sax1;
  bool sax2;
  int part_num;
  int section_axis;
  int section_num;
  int i_lb;
  int j_lb;
  int k_lb;
  int i_ub;
  int j_ub;
  int k_ub;
  int J_lb;
  int J_ub;
  int num_vector_local;
} GMSectionInfo;

//-----------------------------------------------------------------------------

static void GMSectionInfo_create(
  GMSectionInfo* si,
  int i_block,
  int j_block,
  int k_block,
  int section_step,
  int num_vector_local,
  GMEnv* env) {
  GMInsist(si && env);
  GMInsist(i_block >= 0 && i_block < GMEnv_num_block_vector(env));
  GMInsist(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMInsist(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMInsist(num_vector_local >= 0);

  si->num_vector_local = num_vector_local;

  si->is_part1 = gm_is_part1(i_block, j_block, k_block);
  si->is_part3 = gm_is_part3(i_block, j_block, k_block);
  si->is_part2 = ! si->is_part1 && ! si->is_part3;

  si->part_num = si->is_part1 ? 1 :
                 si->is_part2 ? 2 : 3;

  const int num_section_steps = gm_num_section_steps(env, si->part_num);
  GMInsist(section_step>=0);
  GMInsist(section_step<num_section_steps);

  si->section_axis =
    ! si->is_part3 ? 1 /*---j axis---*/ :
    gm_section_axis_part3(i_block, j_block, k_block);

  si->section_num = ! si->is_part3 ? section_step :
                    gm_section_num_part3(i_block, j_block, k_block);

  /*---Define bounding box containing region to be computed---*/

  si->J_lb = gm_J_lo(si->section_num, num_vector_local, si->part_num, env);
  si->J_ub = gm_J_hi(si->section_num, num_vector_local, si->part_num, env);

  si->i_lb = si->section_axis == 0 ? si->J_lb : 0;
  si->j_lb = si->section_axis == 1 ? si->J_lb : 0;
  si->k_lb = si->section_axis == 2 ? si->J_lb : 0;

  si->i_ub = si->section_axis == 0 ? si->J_ub : num_vector_local;
  si->j_ub = si->section_axis == 1 ? si->J_ub : num_vector_local;
  si->k_ub = si->section_axis == 2 ? si->J_ub : num_vector_local;

  si->sax0 = si->section_axis == 0;
  si->sax1 = si->section_axis == 1;
  si->sax2 = si->section_axis == 2;
}

/*
- I_max = is_part1 ? J : nvl; // XXX can work same way if permuted or not
- K_min = is_part3 ? 0 : J + 1; // XXX can work same way if permuted or not
- put in functions for permuted (I, K) and nonpermuted (i, k)
- store I_ub, etc.
- should we always permute axes (I think not - perhaps only if parallel all2all)
- should we always slice into 6 sections (?)
- do lb/ub values apply for part1/2 - is there a permutation issue - ? OK
- should this be cognizant of all2all value

- * deploy section_num usage for part1/2
*/

//-----------------------------------------------------------------------------

static void GMSectionInfo_destroy(
  GMSectionInfo* si,
  GMEnv* env) {
}

//-----------------------------------------------------------------------------

static int GMSectionInfo_k_min(
  GMSectionInfo* si,
  int j,
  GMEnv* env) {
  GMInsist(si && env);
  GMAssert(j >= 0 && j < si->num_vector_local);

  return si->is_part3 ? si->k_lb : j + 1;
}

//-----------------------------------------------------------------------------

static int GMSectionInfo_i_max(
  GMSectionInfo* si,
  int j,
  GMEnv* env) {
  GMInsist(si && env);
  GMAssert(j >= 0 && j < si->num_vector_local);

  return si->is_part1 ? j : si->i_ub;
}

//=============================================================================
/*---Accessors: indexing: (contig) index from coord, 3-way---*/

//-----------------------------------------------------------------------------
/*---elements in a part of a trapezoid, cut orthog to j axis---*/

static size_t gm_trap_size(int j, int nvl) {
  return ( j *(size_t) (j-1) *(size_t) (3*nvl-2*j-2) ) / 6;
}

//-----------------------------------------------------------------------------
/*---elements in a part of a triang, cut orthog to j axis---*/

static size_t gm_triang_size(int j, int nvl) {
  return gm_triang_(nvl) - gm_triang_(nvl-j);
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_index_from_coord_3(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           int k,
                                           GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(i < j && j < k);

  const int nvl = metrics->num_vector_local;

  /* clang-format off */
  const int64_t index = i +
                        (k-j-1)*(size_t)j +
                        gm_trap_size(j, nvl);
  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < (int64_t)metrics->num_elts_local);

  GMAssert(i + metrics->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] % metrics->num_vector);
  GMAssert(j + metrics->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env) ==
           (metrics->coords_global_from_index[index] / metrics->num_vector) %
               metrics->num_vector);
  GMAssert(k + metrics->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] /
               (metrics->num_vector * (size_t)metrics->num_vector));

  return index;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper3way_part1_(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int k,
                                          int i_block,
                                          int j_block,
                                          int k_block,
                                          GMEnv* env) {
  const int nvl = metrics->num_vector_local;

  const int num_section_steps = gm_num_section_steps(env, 1);
  const int section_num = (j * num_section_steps) / nvl;
  GMAssert(metrics->section_num_valid_part1_[section_num]);

  const int64_t elts_offset = metrics->index_offset_section_part1_[section_num];

  /* clang-format off */
  const int64_t index = elts_offset +
                        i +
                        (k-j-1)*(size_t)j +
                        gm_trap_size(j, nvl);
  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < (int64_t)metrics->num_elts_local);

  return index;
  //return GMMetrics_index_from_coord_3(metrics, i, j, k, env);
}

//-----------------------------------------------------------------------------
/*---Faster version of true mod, needed for special situation---*/

static int gm_mod1_(int i, int n) {
  return (i + n) % n;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper3way_part2_(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int k,
                                          int i_block,
                                          int j_block,
                                          int k_block,
                                          GMEnv* env) {
  const int nvl = metrics->num_vector_local;

  const int num_section_steps = gm_num_section_steps(env, 2);
  const int section_num = (j * num_section_steps) / nvl;
  GMAssert(metrics->section_num_valid_part2_[section_num]);

  const int64_t elts_offset = metrics->index_offset_section_part2_[section_num];

  const int num_block = GMEnv_num_block_vector(env);
  const int j_i_offset = gm_mod1_(j_block - i_block, num_block);
  const int block_num_part2 = j_i_offset - 1
      - metrics->phase_block_start_2_[section_num];

  // Packing offset for multiple section blocks for this proc_r and section_num
  const int num_proc_r = GMEnv_num_proc_repl(env);
  const int blocks_offset = block_num_part2 / num_proc_r;

  const size_t section_size = metrics->section_size_part2[section_num];

  // Ordering: outer loop is section num, inner loop is block num.

  /* clang-format off */
  const int64_t index = elts_offset +
                        i + nvl*(
                        (k-j-1) +
                        gm_triang_size(j, nvl) + section_size*(
                        blocks_offset
                        ));
 /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < (int64_t)metrics->num_elts_local);

  return index;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper3way_part3_(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int k,
                                          int i_block,
                                          int j_block,
                                          int k_block,
                                          GMEnv* env) {
  const int nvl = metrics->num_vector_local;

  //const int num_section_steps = 1;
  const int section_num = gm_section_num_part3(i_block, j_block, k_block);

  const int64_t elts_offset = metrics->index_offset_01_;

  const int num_block = GMEnv_num_block_vector(env);
  const int j_i_offset = gm_mod1_(j_block - i_block, num_block);
  const int k_i_offset = gm_mod1_(k_block - i_block, num_block);
  const int block_num_part3 =
    ((num_block-2) * (k_i_offset - 1)) +
    (j_i_offset - 1 - (j_i_offset > k_i_offset))
    - metrics->phase_block_start_3_;

  // Packing offset for multiple blocks for this proc_r
  const int num_proc_r = GMEnv_num_proc_repl(env);
  const int blocks_offset = block_num_part3 / num_proc_r;

  const int section_axis = gm_section_axis_part3(i_block, j_block, k_block);
  const int J_lo = metrics->J_lo_part3_[section_num];
  const int J_wi = metrics->J_wi_part3_[section_num];

  /* clang-format off */
  const int64_t index = elts_offset +
                        i - ( section_axis == 0 ? J_lo : 0 ) +
                            ( section_axis == 0 ? J_wi : nvl ) * (
                        k - ( section_axis == 2 ? J_lo : 0 ) +
                            ( section_axis == 2 ? J_wi : nvl ) * (
                        j - ( section_axis == 1 ? J_lo : 0 ) +
                            ( section_axis == 1 ? J_wi : nvl ) * (
                        (int64_t)blocks_offset
        )));
  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < (int64_t)metrics->num_elts_local);

  return index;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_index_from_coord_all2all_3(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int k,
                                                   int j_block,
                                                   int k_block,
                                                   GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0 && j >= 0 && k >= 0);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(! (GMEnv_proc_num_vector_i(env) == j_block &&
              GMEnv_proc_num_vector_i(env) != k_block));
  GMAssert(! (GMEnv_proc_num_vector_i(env) == k_block &&
              GMEnv_proc_num_vector_i(env) != j_block));
  /*---WARNING: these conditions are not exhaustive---*/

  const int i_block = GMEnv_proc_num_vector_i(env);

  const int64_t index = j_block == i_block && k_block == i_block ?
    GMMetrics_helper3way_part1_(metrics, i, j, k,
                                i_block, j_block, k_block, env) :
                 j_block == k_block ?
    GMMetrics_helper3way_part2_(metrics, i, j, k,
                                i_block, j_block, k_block, env) :
    GMMetrics_helper3way_part3_(metrics, i, j, k,
                                i_block, j_block, k_block, env);

  GMAssert(index >= 0);
  GMAssert(index < (int64_t)metrics->num_elts_local);

  GMAssert(metrics->coords_global_from_index[index] %
             (metrics->num_vector_local * (size_t)GMEnv_num_block_vector(env)) ==
           i + i_block * (size_t)metrics->num_vector_local);

  GMAssert((metrics->coords_global_from_index[index] /
            (metrics->num_vector_local * (size_t)GMEnv_num_block_vector(env))) %
               (metrics->num_vector_local * GMEnv_num_block_vector(env)) ==
           j + j_block * (size_t)metrics->num_vector_local);

  GMAssert((metrics->coords_global_from_index[index] /
            (metrics->num_vector_local * (size_t)GMEnv_num_block_vector(env))) /
               (metrics->num_vector_local * GMEnv_num_block_vector(env)) ==
           k + k_block * (size_t)metrics->num_vector_local);

  return index;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper3way_part1_permuted_(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int i_block,
    int j_block,
    int k_block,
    GMEnv* env) {

  return GMMetrics_helper3way_part1_(metrics, I, J, K,
                                     i_block, j_block, k_block, env);
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper3way_part2_permuted_(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int i_block,
    int j_block,
    int k_block,
    GMEnv* env) {

  return GMMetrics_helper3way_part2_(metrics, I, J, K,
                                     i_block, j_block, k_block, env);
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper3way_part3_permuted_(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int i_block,
    int j_block,
    int k_block,
    GMEnv* env) {
  const int nvl = metrics->num_vector_local;

  //const int num_section_steps = 1;
  const int section_num = gm_section_num_part3(i_block, j_block, k_block);

  const int64_t elts_offset = metrics->index_offset_01_;

  const int num_block = GMEnv_num_block_vector(env);
  const int j_i_offset = gm_mod1_(j_block - i_block, num_block);
  const int k_i_offset = gm_mod1_(k_block - i_block, num_block);
  const int block_num_part3 =
    ((num_block-2) * (k_i_offset - 1)) +
    (j_i_offset - 1 - (j_i_offset > k_i_offset))
    - metrics->phase_block_start_3_;

  const int num_proc_r = GMEnv_num_proc_repl(env);
  const int blocks_offset = block_num_part3 / num_proc_r;

  const int J_lo = metrics->J_lo_part3_[section_num];
  const int J_wi = metrics->J_wi_part3_[section_num];

  /* clang-format off */

  const int64_t index = elts_offset +
                        I + nvl * (
                        K + nvl * (
                        J - J_lo + J_wi * (
                        (int64_t)blocks_offset
                        )));

  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < (int64_t)metrics->num_elts_local);

  return index;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_index_from_coord_all2all_3_permuted(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && J >= 0 && K >= 0);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(! (GMEnv_proc_num_vector_i(env) == j_block &&
              GMEnv_proc_num_vector_i(env) != k_block));
  GMAssert(! (GMEnv_proc_num_vector_i(env) == k_block &&
              GMEnv_proc_num_vector_i(env) != j_block));
  /*---WARNING: these conditions are not exhaustive---*/

  const int i_block = GMEnv_proc_num_vector_i(env);

  size_t index = j_block == i_block && k_block == i_block ?
    GMMetrics_helper3way_part1_permuted_(metrics, I, J, K,
                                i_block, j_block, k_block, env) :
                 j_block == k_block ?
    GMMetrics_helper3way_part2_permuted_(metrics, I, J, K,
                                i_block, j_block, k_block, env) :
    GMMetrics_helper3way_part3_permuted_(metrics, I, J, K,
                                         i_block, j_block, k_block, env);

  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);

#ifdef GM_ASSERTIONS_ON

  int section_step = 0;

  if (gm_is_part3(i_block, j_block, k_block)) {
    section_step = 0;
  } else if (gm_is_part1(i_block, j_block, k_block)) {
    const int num_section_steps = gm_num_section_steps(env, 1);
    const int section_num = (J * num_section_steps) / metrics->num_vector_local;
    section_step = section_num;
  } else {
    const int num_section_steps = gm_num_section_steps(env, 2);
    const int section_num = (J * num_section_steps) / metrics->num_vector_local;
    section_step = section_num;
  }

  GMSectionInfo si_value_;
  GMSectionInfo* si = &si_value_;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  const bool is_part3 = si->is_part3;
  const bool no_perm = ! is_part3;
  const bool sax0 = si->section_axis == 0;
  const bool sax1 = si->section_axis == 1;
  //const bool sax2 = si->section_axis == 2;

  /* clang-format off */
  const int i = no_perm ? I :
                sax0 ?    J :
                sax1 ?    I :
             /* sax2 ?*/  K;
  const int j = no_perm ? J :
                sax0 ?    K :
                sax1 ?    J :
             /* sax2 ?*/  I;
  const int k = no_perm ? K :
                sax0 ?    I :
                sax1 ?    K :
             /* sax2 ?*/  J;
  /* clang-format on */

  GMSectionInfo_destroy(si, env);

  GMAssert(metrics->coords_global_from_index[index] % metrics->num_vector ==
           i + i_block * (size_t)metrics->num_vector_local);

  GMAssert((metrics->coords_global_from_index[index] / metrics->num_vector) %
               metrics->num_vector ==
           j + j_block * (size_t)metrics->num_vector_local);

  GMAssert((metrics->coords_global_from_index[index] / metrics->num_vector) /
               metrics->num_vector ==
           k + k_block * (size_t)metrics->num_vector_local);
#endif

  return index;
}

//-----------------------------------------------------------------------------

typedef struct {
  bool is_initialized;
  int I;
  int K;
  size_t index;
} GMIndexCache;

//-----------------------------------------------------------------------------

static size_t GMMetrics_index_from_coord_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMIndexCache* index_cache,
    GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && J >= 0 && K >= 0);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(! (GMEnv_proc_num_vector_i(env) == j_block &&
              GMEnv_proc_num_vector_i(env) != k_block));
  GMAssert(! (GMEnv_proc_num_vector_i(env) == k_block &&
              GMEnv_proc_num_vector_i(env) != j_block));
  /*---WARNING: these conditions are not exhaustive---*/

  if (index_cache->is_initialized && K == index_cache->K) {
      const size_t index = index_cache->index + (I-index_cache->I);
      index_cache->index = index;
      index_cache->I = I;
      return index;
  }

  const size_t index = GMMetrics_index_from_coord_all2all_3_permuted(
    metrics, I, J, K, j_block, k_block, env);

  index_cache->I = I;
  index_cache->K = K;
  index_cache->index = index;
  index_cache->is_initialized = true;

  return index;
}
//=============================================================================
/*---Accessors: indexing: global coord from (contig) index: 3-way---*/

static int GMMetrics_coord0_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index >= 0 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);

  const size_t i64 = metrics->coords_global_from_index[index] %
                     metrics->num_vector;
  const int i = (int)i64;
  GMAssert((size_t)i == i64);

  return i;
}

//-----------------------------------------------------------------------------

static int GMMetrics_coord1_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index >= 0 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);

  const size_t j64 =
      (metrics->coords_global_from_index[index] / metrics->num_vector) %
      metrics->num_vector;
  const int j = (int)j64;
  GMAssert((size_t)j == j64);

  return j;
}

//-----------------------------------------------------------------------------

static int GMMetrics_coord2_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index >= 0 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);

  const size_t k64 = metrics->coords_global_from_index[index] /
                     (metrics->num_vector * (size_t)metrics->num_vector);
  const int k = (int)k64;
  GMAssert((size_t)k == k64);

  return k;
}

//=============================================================================

#endif // _gm_metrics_3way_indexing_hh_

//-----------------------------------------------------------------------------
