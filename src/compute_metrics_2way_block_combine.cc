//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way_block_combine.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Combine numerators and denominators, 2-way, for a single block.
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

#include "cstdint"
#include "inttypes.h"

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "compute_metrics_2way_block_combine.hh"

//=============================================================================
/*---Combine nums and denoms on CPU to get final result, 2-way Czek---*/

void gm_compute_2way_proc_combine_czek_(
  GMMetrics* metrics,
  GMMirroredBuf* metrics_buf,
  const GMVectorSums* vector_sums_left,
  const GMVectorSums* vector_sums_right,
  int j_block,
   bool do_compute_triang_only,
   GMEnv* env) {

  GMInsist(metrics && metrics_buf);
  GMInsist(vector_sums_left && vector_sums_right && env);
  GMInsist(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_2);

  // NOTE: here and elsewhere, vector_sums_left and vector_sums_right
  // may sometimes be aliases.  Since they are read-only, there should
  // be no danger of the compiler generating incorrect code.
  // (cf. https://en.wikipedia.org/wiki/Pointer_aliasing)
  // However by accounting for this one might be able to in principle
  // remove a load instruction to improve performance.

  const int nvl = metrics->num_vector_local;
  const GMVectorSums* vs_l = vector_sums_left;
  const GMVectorSums* vs_r = vector_sums_right;

  /*---For CPU case, copy numerator out of metrics struct which is temporarily
       holding numerators.
       For GPU case, directly access the metrics_buf holding the numerators.
  ---*/

  /*----------------------------------------*/
  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU && GMEnv_all2all(env)) {
    /*----------------------------------------*/

    for (int j = 0; j < nvl; ++j) {
      const GMFloat vs_j = GMVectorSums_sum(vs_r, j, env);
      const int i_max = do_compute_triang_only ? j : nvl;
      for (int i = 0; i < i_max; ++i) {
        const GMFloat numer =
            GMMetrics_float_get_all2all_2(metrics, i, j, j_block, env);
        const GMFloat vs_i = GMVectorSums_sum(vs_l, i, env);
        const GMFloat denom = vs_i < vs_j ? vs_i + vs_j : vs_j + vs_i;
        const GMFloat multiplier = (GMFloat)2;
        const GMFloat value = (multiplier * numer) / denom;
        GMMetrics_float_set_all2all_2(metrics, i, j, j_block, value, env);
//const size_t index = GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
//const size_t vi_ = GMMetrics_coord_global_from_index(metrics, index, 0, env);
//const size_t vj_ = GMMetrics_coord_global_from_index(metrics, index, 1, env);
//if (vi_==364001 && vj_==8714000) printf("===========================+++ %.15e %.15e %.20e\n", numer, denom, value);
//if (denom != (GMFloat)(size_t)denom) printf("Bad d %zu %zu %.20e\n", vi_, vj_, denom);
//if (numer != (GMFloat)(size_t)numer) printf("Bad n %zu %zu %.20e\n", vi_, vj_, numer);
      } /*---for i---*/
      metrics->num_elts_local_computed += i_max;
    }   /*---for j---*/
        /*---TODO: here and elsewhere check for unlikely case denom is/nearly
         * zero---*/

    /*----------------------------------------*/
  } else if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    for (int j = 0; j < nvl; ++j) {
      const GMFloat vs_j = GMVectorSums_sum(vs_r, j, env);
      const int i_max = j;
      for (int i = 0; i < i_max; ++i) {
        const GMFloat numer = GMMetrics_float_get_2(metrics, i, j, env);
        const GMFloat vs_i = GMVectorSums_sum(vs_l, i, env);
        const GMFloat denom = vs_i < vs_j ?  vs_i + vs_j : vs_j + vs_i;
        const GMFloat multiplier = (GMFloat)2;
        const GMFloat value = (multiplier * numer) / denom;
        GMMetrics_float_set_2(metrics, i, j, value, env);
      } /*---for i---*/
      metrics->num_elts_local_computed += i_max;
    }   /*---for j---*/

    /*----------------------------------------*/
  } else if (GMEnv_all2all(env)) {
    /*----------------------------------------*/

    if (do_compute_triang_only) {
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const GMFloat vs_j = GMVectorSums_sum(vs_r, j, env);
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const GMFloat numer =
            GMMirroredBuf_elt<GMFloat>(metrics_buf, i, j);
          const GMFloat vs_i = GMVectorSums_sum(vs_l, i, env);
          const GMFloat denom = vs_i < vs_j ? vs_i + vs_j : vs_j + vs_i;
          const GMFloat multiplier = (GMFloat)2;
          const GMFloat value = (multiplier * numer) / denom;
          GMMetrics_float_set_all2all_2(metrics, i, j, j_block, value, env);
        } /*---for i---*/
      }   /*---for j---*/
      for (int j = 0; j < nvl; ++j) {
        const int i_max = j;
        metrics->num_elts_local_computed += i_max;
      }   /*---for j---*/
    } else {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const GMFloat vs_j = GMVectorSums_sum(vs_r, j, env);
          const GMFloat numer =
            GMMirroredBuf_elt<GMFloat>(metrics_buf, i, j);
          const GMFloat vs_i = GMVectorSums_sum(vs_l, i, env);
          const GMFloat denom = vs_i < vs_j ? vs_i + vs_j : vs_j + vs_i;
          const GMFloat multiplier = (GMFloat)2;
          const GMFloat value = (multiplier * numer) / denom;
          GMMetrics_float_set_all2all_2(metrics, i, j, j_block, value, env);
        } /*---for i---*/
      }   /*---for j---*/
      metrics->num_elts_local_computed += nvl * (size_t)nvl;
    }

    /*----------------------------------------*/
  } else {
    /*----------------------------------------*/

    #pragma omp parallel for schedule(dynamic,1000)
    for (int j = 0; j < nvl; ++j) {
      const GMFloat vs_j = GMVectorSums_sum(vs_r, j, env);
      const int i_max = j;
      for (int i = 0; i < i_max; ++i) {
        const GMFloat numer =
          GMMirroredBuf_elt<GMFloat>(metrics_buf, i, j);
        const GMFloat vs_i = GMVectorSums_sum(vs_l, i, env);
        const GMFloat denom = vs_i < vs_j ? vs_i + vs_j : vs_j + vs_i;
        const GMFloat multiplier = (GMFloat)2;
        const GMFloat value = (multiplier * numer) / denom;
        GMMetrics_float_set_2(metrics, i, j, value, env);
      } /*---for i---*/
    }   /*---for j---*/
    for (int j = 0; j < nvl; ++j) {
      const int i_max = j;
      metrics->num_elts_local_computed += i_max;
    }   /*---for j---*/

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================
/*---Combine nums and denoms on CPU to get final result, 2-way CCC---*/

void gm_compute_2way_proc_combine_ccc_(
  GMMetrics* metrics,
  GMMirroredBuf* metrics_buf,
  const GMVectorSums* vector_sums_left,
  const GMVectorSums* vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  GMEnv* env) {

  GMInsist(metrics && metrics_buf);
  GMInsist(vector_sums_left && vector_sums_right && env);
  GMInsist(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_2);

  const int nvl = metrics->num_vector_local;
  const GMVectorSums* vs_l = vector_sums_left;
  const GMVectorSums* vs_r = vector_sums_right;

  /*---Copy from metrics_buffer for GPU case; perform checks---*/

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    /*--------------------*/
    if (GMEnv_all2all(env)) {
      /*--------------------*/

      if (do_compute_triang_only) {
        #pragma omp parallel for schedule(dynamic,1000)
        for (int j = 0; j < nvl; ++j) {
          const int i_max = j;
          for (int i = 0; i < i_max; ++i) {
            const GMTally2x2 value =
              GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
            GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, value, env);
#if 1
            // ISSUE: this check may increase runtime nontrivially
            if (! env->sparse) {
              // 4-sum check.
              const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
              const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
              const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
              const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
              const bool error1 = (uint64_t)r00 + (uint64_t)r01 +
                                  (uint64_t)r10 + (uint64_t)r11 !=
                       (uint64_t)(4 * metrics->num_field_active);
              if (error1) {
                const size_t index = GMMetrics_index_from_coord_all2all_2(metrics,
                  i, j, j_block, env);
                const size_t coords = metrics->coords_global_from_index[index];
                printf("Error: r00 %" PRIu64 " r01 %" PRIu64 " r10 %" PRIu64 " r11 %" PRIu64 " m %" PRIu64 " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)(metrics->num_field_active),
                       coords, GMEnv_proc_num(env));
                GMInsist((! error1) && "Violation of algorithm computational invariant.");
              }
              // 2-sum check.
              const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
              const bool error2 = (uint64_t)r10 + (uint64_t)r11 !=
                                  (uint64_t)(2 * si1);
              if (error2) {
                const size_t index = GMMetrics_index_from_coord_all2all_2(metrics,
                  i, j, j_block, env);
                const size_t coords = metrics->coords_global_from_index[index];
                printf("Error: r00 %" PRIu64 " r01 %" PRIu64 " r10 %" PRIu64 " r11 %" PRIu64 " si1 %" PRIu64 " actual %" PRIu64 " expected %" PRIu64 " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)si1,
                       (uint64_t)r10 + (uint64_t)r11, (uint64_t)(2 * si1),
                       coords, GMEnv_proc_num(env));
                GMInsist((! error2) && "Violation of algorithm computational invariant.");
              }
              // 2-sum check.
              const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
              const bool error3 = (uint64_t)r01 + (uint64_t)r11 !=
                                  (uint64_t)(2 * sj1);
              if (error3) {
                const size_t index = GMMetrics_index_from_coord_all2all_2(metrics,
                  i, j, j_block, env);
                const size_t coords = metrics->coords_global_from_index[index];
                printf("Error: r00 %" PRIu64 " r01 %" PRIu64 " r10 %" PRIu64 " r11 %" PRIu64 " sj1 %" PRIu64 " actual %" PRIu64 " expected %" PRIu64 " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)sj1,
                       (uint64_t)r01 + (uint64_t)r11, (uint64_t)(2 * sj1),
                       coords, GMEnv_proc_num(env));
                GMInsist((! error3) && "Violation of algorithm computational invariant.");
              }
            }
#endif
#ifdef GM_ASSERTIONS_ON
            if (! env->sparse) {
              // 4-sum check.
              const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
              const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
              const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
              const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
              GMAssert((uint64_t)r00 + (uint64_t)r01 + (uint64_t)r10 +
                           (uint64_t)r11 ==
                       (uint64_t)(4 * metrics->num_field_active));
              // 2-sum checks.
              const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
              const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
              GMAssert((uint64_t)r10 + (uint64_t)r11 == (uint64_t)(2 * si1));
              GMAssert((uint64_t)r01 + (uint64_t)r11 == (uint64_t)(2 * sj1));
            }
#endif
          } /*---for i---*/
        }   /*---for j---*/
      } else {
        // don't use collapse because of overflow for large sizes
        //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
        #pragma omp parallel for schedule(dynamic,1000)
        for (int j = 0; j < nvl; ++j) {
          for (int i = 0; i < nvl; ++i) {
            const GMTally2x2 value =
              GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
            GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, value, env);
#if 1
            // ISSUE: this check may increase runtime nontrivially
            if (! env->sparse) {
              // 4-sum check.
              const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
              const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
              const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
              const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
              const bool error1 = (uint64_t)r00 + (uint64_t)r01 +
                                  (uint64_t)r10 + (uint64_t)r11 !=
                       (uint64_t)(4 * metrics->num_field_active);
              if (error1) {
                const size_t index = GMMetrics_index_from_coord_all2all_2(metrics,
                  i, j, j_block, env);
                const size_t coords = metrics->coords_global_from_index[index];
                printf("Error: r00 %" PRIu64 " r01 %" PRIu64 " r10 %" PRIu64 " r11 %" PRIu64 " m %" PRIu64 " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)(metrics->num_field_active),
                       coords, GMEnv_proc_num(env));
                GMInsist((! error1) && "Violation of algorithm computational invariant.");
              }
              // 2-sum check.
              const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
              const bool error2 = (uint64_t)r10 + (uint64_t)r11 !=
                                  (uint64_t)(2 * si1);
              if (error2) {
                const size_t index = GMMetrics_index_from_coord_all2all_2(metrics,
                  i, j, j_block, env);
                const size_t coords = metrics->coords_global_from_index[index];
                printf("Error: r00 %" PRIu64 " r01 %" PRIu64 " r10 %" PRIu64 " r11 %" PRIu64 " si1 %" PRIu64 " actual %" PRIu64 " expected %" PRIu64 " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)si1,
                       (uint64_t)r10 + (uint64_t)r11, (uint64_t)(2 * si1),
                       coords, GMEnv_proc_num(env));
                GMInsist((! error2) && "Violation of algorithm computational invariant.");
              }
              // 2-sum check.
              const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
              const bool error3 = (uint64_t)r01 + (uint64_t)r11 !=
                                  (uint64_t)(2 * sj1);
              if (error3) {
                const size_t index = GMMetrics_index_from_coord_all2all_2(metrics,
                  i, j, j_block, env);
                const size_t coords = metrics->coords_global_from_index[index];
                printf("Error: r00 %" PRIu64 " r01 %" PRIu64 " r10 %" PRIu64 " r11 %" PRIu64 " sj1 %" PRIu64 " actual %" PRIu64 " expected %" PRIu64 " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)sj1,
                       (uint64_t)r01 + (uint64_t)r11, (uint64_t)(2 * sj1),
                       coords, GMEnv_proc_num(env));
                GMInsist((! error3) && "Violation of algorithm computational invariant.");
              }
            }
#endif
#ifdef GM_ASSERTIONS_ON
            if (! env->sparse) {
              // 4-sum check.
              const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
              const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
              const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
              const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
              GMAssert((uint64_t)r00 + (uint64_t)r01 + (uint64_t)r10 +
                           (uint64_t)r11 ==
                       (uint64_t)(4 * metrics->num_field_active));
              // 2-sum checks.
              const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
              const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
              GMAssert((uint64_t)r10 + (uint64_t)r11 == (uint64_t)(2 * si1));
              GMAssert((uint64_t)r01 + (uint64_t)r11 == (uint64_t)(2 * sj1));
            }
#endif
          } /*---for i---*/
        }   /*---for j---*/
     }

      /*--------------------*/
    } else /*---(! GMEnv_all2all(env))---*/ {
      /*--------------------*/
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const int i_max = do_compute_triang_only ? j : nvl;
        for (int i = 0; i < i_max; ++i) {
          const GMTally2x2 value =
              GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
          GMMetrics_tally2x2_set_2(metrics, i, j, value, env);
#ifdef GM_ASSERTIONS_ON
          if (! env->sparse) {
            // 4-sum check.
            const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
            const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
            const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
            const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
            GMAssert((uint64_t)r00 + (uint64_t)r01 + (uint64_t)r10 +
                         (uint64_t)r11 ==
                     (uint64_t)(4 * metrics->num_field_active));
            // 2-sum checks.
            const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
            const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
            GMAssert((uint64_t)r10 + (uint64_t)r11 == (uint64_t)(2 * si1));
            GMAssert((uint64_t)r01 + (uint64_t)r11 == (uint64_t)(2 * sj1));
          }
#endif
        } /*---for i---*/
      }   /*---for j---*/
      /*--------------------*/
    } /*---if---*/
    /*--------------------*/
  }

  /*---Compute multipliers---*/

  /*--------------------*/
  if (GMEnv_all2all(env)) {
    /*--------------------*/

    if (do_compute_triang_only) {
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          GMMetrics_float2_S_set_all2all_2(metrics, i, j, j_block, si1_sj1, env);
          if (env->sparse) {
            const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_r, j, env);
            const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_l, i, env);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            GMMetrics_float2_C_set_all2all_2(metrics, i, j, j_block, ci_cj, env);
          } /*---if sparse---*/
        }   /*---for i---*/
      }   /*---for j---*/
      for (int j = 0; j < nvl; ++j) {
        const int i_max = j;
        metrics->num_elts_local_computed += i_max;
      }   /*---for j---*/
    } else {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
          const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          GMMetrics_float2_S_set_all2all_2(metrics, i, j, j_block, si1_sj1, env);
          if (env->sparse) {
            const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_r, j, env);
            const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_l, i, env);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            GMMetrics_float2_C_set_all2all_2(metrics, i, j, j_block, ci_cj, env);
          } /*---if sparse---*/
        }   /*---for i---*/
      }   /*---for j---*/
      metrics->num_elts_local_computed += nvl * (size_t)nvl;
   }

    /*--------------------*/
  } else /*---(! GMEnv_all2all(env))---*/ {
    /*--------------------*/
    #pragma omp parallel for schedule(dynamic,1000)
    for (int j = 0; j < nvl; ++j) {
      const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
      const int i_max = do_compute_triang_only ? j : nvl;
      for (int i = 0; i < i_max; ++i) {
        const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
        const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
        GMMetrics_float2_S_set_2(metrics, i, j, si1_sj1, env);
        if (env->sparse) {
          const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_r, j, env);
          const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_l, i, env);
          const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
          GMMetrics_float2_C_set_2(metrics, i, j, ci_cj, env);
        } /*---if sparse---*/
      } /*---for i---*/
    }   /*---for j---*/
    for (int j = 0; j < nvl; ++j) {
      const int i_max = do_compute_triang_only ? j : nvl;
      metrics->num_elts_local_computed += i_max;
    }   /*---for j---*/
    /*--------------------*/
  } /*---if---*/
  /*--------------------*/
}

//=============================================================================
/*---Combine nums and denoms on CPU to get final result, 2-way DUO---*/

void gm_compute_2way_proc_combine_duo_(
  GMMetrics* metrics,
  GMMirroredBuf* metrics_buf,
  const GMVectorSums* vector_sums_left,
  const GMVectorSums* vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  GMEnv* env) {

  GMInsist(metrics && metrics_buf);
  GMInsist(vector_sums_left && vector_sums_right && env);
  GMInsist(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_2);

  const int nvl = metrics->num_vector_local;
  const GMVectorSums* vs_l = vector_sums_left;
  const GMVectorSums* vs_r = vector_sums_right;

  /*---Copy from metrics_buffer for GPU case---*/

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    /*--------------------*/
    if (GMEnv_all2all(env)) {
      /*--------------------*/

      if (do_compute_triang_only) {
        #pragma omp parallel for schedule(dynamic,1000)
        for (int j = 0; j < nvl; ++j) {
          const int i_max = j;
          for (int i = 0; i < i_max; ++i) {
            const GMTally2x2 value =
              GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
            GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, value, env);
          } /*---for i---*/
        }   /*---for j---*/
      } else {
        // don't use collapse because of overflow for large sizes
        //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
        #pragma omp parallel for schedule(dynamic,1000)
        for (int j = 0; j < nvl; ++j) {
          for (int i = 0; i < nvl; ++i) {
            const GMTally2x2 value =
              GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
            GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, value, env);
          } /*---for i---*/
        }   /*---for j---*/
     }

      /*--------------------*/
    } else /*---(! GMEnv_all2all(env))---*/ {
      /*--------------------*/
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const int i_max = do_compute_triang_only ? j : nvl;
        for (int i = 0; i < i_max; ++i) {
          const GMTally2x2 value =
              GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
          GMMetrics_tally2x2_set_2(metrics, i, j, value, env);
        } /*---for i---*/
      }   /*---for j---*/
      /*--------------------*/
    } /*---if---*/
    /*--------------------*/
  }

  /*---Compute multipliers---*/

  /*--------------------*/
  if (GMEnv_all2all(env)) {
    /*--------------------*/

    if (do_compute_triang_only) {
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          GMMetrics_float2_S_set_all2all_2(metrics, i, j, j_block, si1_sj1, env);
          if (env->sparse) {
            const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_r, j, env);
            const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_l, i, env);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            GMMetrics_float2_C_set_all2all_2(metrics, i, j, j_block, ci_cj, env);
          } /*---if sparse---*/
        }   /*---for i---*/
      }   /*---for j---*/
      for (int j = 0; j < nvl; ++j) {
        const int i_max = j;
        metrics->num_elts_local_computed += i_max;
      }   /*---for j---*/
    } else {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
          const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          GMMetrics_float2_S_set_all2all_2(metrics, i, j, j_block, si1_sj1, env);
          if (env->sparse) {
            const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_r, j, env);
            const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_l, i, env);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            GMMetrics_float2_C_set_all2all_2(metrics, i, j, j_block, ci_cj, env);
          } /*---if sparse---*/
        }   /*---for i---*/
      }   /*---for j---*/
      metrics->num_elts_local_computed += nvl * (size_t)nvl;
   }

    /*--------------------*/
  } else /*---(! GMEnv_all2all(env))---*/ {
    /*--------------------*/
    #pragma omp parallel for schedule(dynamic,1000)
    for (int j = 0; j < nvl; ++j) {
      const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
      const int i_max = do_compute_triang_only ? j : nvl;
      for (int i = 0; i < i_max; ++i) {
        const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
        const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
        GMMetrics_float2_S_set_2(metrics, i, j, si1_sj1, env);
        if (env->sparse) {
          const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_r, j, env);
          const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_l, i, env);
          const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
          GMMetrics_float2_C_set_2(metrics, i, j, ci_cj, env);
        } /*---if sparse---*/
      } /*---for i---*/
    }   /*---for j---*/
    for (int j = 0; j < nvl; ++j) {
      const int i_max = do_compute_triang_only ? j : nvl;
      metrics->num_elts_local_computed += i_max;
    }   /*---for j---*/
    /*--------------------*/
  } /*---if---*/
  /*--------------------*/
}

//=============================================================================
/*---Combine nums and denoms on CPU to get final result, 2-way generic---*/

void gm_compute_2way_proc_combine(
  GMMetrics* metrics,
  GMMirroredBuf* metrics_buf,
  const GMVectorSums* vector_sums_left,
  const GMVectorSums* vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  GMEnv* env) {

  GMInsist(metrics && metrics_buf);
  GMInsist(vector_sums_left && vector_sums_right && env);
  GMInsist(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_2);

  switch (GMEnv_metric_type(env)) {
    case GM_METRIC_TYPE_CZEK: {
      gm_compute_2way_proc_combine_czek_(metrics, metrics_buf,
                                         vector_sums_left, vector_sums_right,
                                         j_block, do_compute_triang_only, env);
    } break;
    case GM_METRIC_TYPE_CCC: {
      gm_compute_2way_proc_combine_ccc_(metrics, metrics_buf,
                                        vector_sums_left, vector_sums_right,
                                        j_block, do_compute_triang_only, env);
    } break;
    case GM_METRIC_TYPE_DUO: {
      gm_compute_2way_proc_combine_duo_(metrics, metrics_buf,
                                        vector_sums_left, vector_sums_right,
                                        j_block, do_compute_triang_only, env);
    } break;
    default:
      GMInsistInterface(env, false && "Selected metric_type unimplemented.");
  } /*---case---*/
}

//-----------------------------------------------------------------------------
