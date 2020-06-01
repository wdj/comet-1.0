//-----------------------------------------------------------------------------
/*!
 * \file   test_problems.cc
 * \author Wayne Joubert
 * \date   Mon Aug  7 17:02:51 EDT 2017
 * \brief  Generator for synthetic test problems.
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
//#include "stdlib.h"
//#include "stddef.h"
//#include "string.h"
//#include "float.h"
//#include "errno.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"

#include "driver.hh"
#include "test_problems.hh"

//=============================================================================
/*---Set the entries of the vectors---*/

void set_vectors_random_(GMVectors* vectors, int verbosity, GMEnv* env) {
  GMInsist(vectors && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  const size_t nva = vectors->dm->num_vector_active;
  const size_t nfa = vectors->dm->num_field_active;

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        /*---Fill pad vectors with copies of the last vector---*/
        // By construction, active vectors are packed for lower procs.
        const size_t vector_capped = gm_min_i8(vector, nva);
        int fl = 0;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          if (field >= vectors->num_field_active) {
            continue; // These entries will be padded to zero elsewhere.
          }
          /*---Compute element unique id---*/
          const size_t uid = field + nfa * vector_capped;
          /*---Generate large random number---*/
          size_t rand1 = uid;
          rand1 = gm_randomize(rand1);
          rand1 = gm_randomize(rand1);
          size_t rand2 = uid;
          rand2 = gm_randomize(rand2);
          rand2 = gm_randomize(rand2);
          rand2 = gm_randomize(rand2);
          const size_t rand_max = gm_randomize_max();
          size_t rand_value = rand1 + rand_max * rand2;
          /*---Reduce so that after summing num_field times the integer
               still exactly representable by floating point type---*/
          const size_t rand_max2 = rand_max * rand_max;
          const int log2_num_summands_3way_numer = 2;
          const int shift_amount1 = gm_max_i8(0,
             gm_log2(log2_num_summands_3way_numer * rand_max2 * nfa)
             - gm_mant_dig<GMFloat>() + 1);
          // Account for cast to float in magma Volta version.
          const int shift_amount2 = gm_max_i8(0,
                             gm_log2(rand_max2) - gm_mant_dig<float>() + 1);
          const int shift_amount = gm_max_i8(shift_amount1, shift_amount2);
          //const int shift_amount = gm_log2(log2_num_summands_3way_numer*
          //                                 rand_max2*nfa)
          //                         - mant_dig;
          rand_value >>= shift_amount > 0 ? shift_amount : 0;
          /*---Store---*/
          GMFloat float_value = (GMFloat)rand_value;
          GMInsist((size_t)float_value == rand_value);
          GMInsist(float_value * vectors->num_field_active <
                         ((size_t)1)<<gm_mant_dig<GMFloat>());
          GMVectors_float_set(vectors, fl, vl, float_value, env);
        } /*---field_local---*/
      }   /*---vector_local---*/
      /*---Print---*/
//TODO: move this
      if (verbosity > 2) {
        GMVectors_print(vectors, env);
      }
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {

        size_t vector = vl +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        /*---Fill pad vectors with copies of the last vector---*/
        const size_t vector_capped = gm_min_i8(vector, nva);

        // XXX
        /*---Fill pad vectors with copies of the last vector---*/
        // const int nval = vectors->dm->num_vector_active_local;
        // const size_t vector_capped = gm_min_i8(vl, nval) +
        //   nval * (size_t)GMEnv_proc_num_vector_i(env);

        for (int fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          if (field >= vectors->num_field_active) {
            continue; // These entries will be padded to zero elsewhere.
          }
          /*---Compute element unique id---*/
          const size_t uid = field + vectors->num_field_active * vector_capped;
          size_t index = uid;
          /*---Randomize---*/
          index = gm_randomize(index);
          index = gm_randomize(index);
          /*---Calculate random number between 0 and 3---*/
          const float float_rand_value = index / (float)gm_randomize_max();
          /*---Create 2-bit value - make extra sure less than 4---*/
          GMBits2 value = (int)((4. - 1e-5) * float_rand_value);
          /*---Store---*/
          GMVectors_bits2_set(vectors, fl, vl, value, env);
        } /*---fl---*/
      }   /*---vl---*/
      /*---Print---*/
//TODO: move this
      if (verbosity > 2) {
        GMVectors_print(vectors, env);
      }
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
}

//-----------------------------------------------------------------------------

static size_t perm_shuffle(size_t key, size_t i, size_t n) {
  GMAssert((key & (~(size_t)1)) == 0);
  GMAssert(i>=0 && i<n);
  GMAssert(n>=0);

  // For an integer between 0 and n-1, permute it to another such integer.
  // The permutation choice is specified by 1 bit.
  // For an ascending sequence of integers, first output the even values,
  // then the odd values, for key=0. If key=1, same with even/odd reversed.

  const size_t nhalf = (n+1-key)/2;
  const size_t result = i < nhalf ? 2*i + key : 2*(i-nhalf) + 1 - key;
  GMAssert(result>=0 && result<n);
  return result;
}

//-----------------------------------------------------------------------------

enum {NUM_SHUFFLE = 3};

// This is not a totally (pseudo)random permutation.  However it does
// have the advantage that it can be computed formulaically and quickly.

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")

static size_t perm(size_t key, size_t i, size_t n) {
  GMAssert((key & (~(size_t)((1<<NUM_SHUFFLE)-1))) == 0);
  GMAssert(i>=0 && i<n);
  GMAssert(n>=0);

  // For an integer between 0 and n-1, permute it to another such integer.
  // The permutation choice is specified by NUM_SHUFFLE bits.

  size_t result = i;
  size_t key_resid = key;
  for (int shuffle_num = 0; shuffle_num < NUM_SHUFFLE; ++shuffle_num) {
    result = perm_shuffle(key_resid&1, result, n);
    key_resid >>= 1;
  }
  GMAssert(result>=0 && result<n);
  return result;
}

#pragma GCC pop_options

//-----------------------------------------------------------------------------

void set_vectors_analytic_(GMVectors* vectors, int verbosity, GMEnv* env) {
  GMInsist(vectors && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  const size_t nfa = vectors->num_field_active;
  const size_t nva = vectors->dm->num_vector_active;

  // Upper bound on integer representable exactly by floating point type.
  // Account for cast to float in magma Volta version.
  const size_t max_float = ((size_t)1) <<
    (GMEnv_data_type_vectors(env) == GM_DATA_TYPE_FLOAT ?
     gm_mant_dig<float>() : gm_mant_dig<GMFloat>());
  // Czek account for number of terms summed in denom or num
  const size_t overflow_limit =
    GMEnv_data_type_vectors(env) != GM_DATA_TYPE_FLOAT ? 1 :
    GMEnv_num_way(env) == GM_NUM_WAY_2 ? 2 : 4;
  // Sum nfa times down the vector, is it still exact.
  const size_t value_limit = (max_float - 1) / (overflow_limit * nfa);

  const size_t value_min = 1;
  //const size_t value_max = (nva+value_min) < value_limit ?
  //                         (nva+value_min) : value_limit;
  const size_t value_max = gm_min_i8(value_min+nva, value_limit);

  // The elements of a single permuted vector are partitioned into
  // "groups", with all elements in a group contiguous and having
  // the same value.
  // By keeping the number of groups (here = 8) much smaller than
  // the vector length, the calculation of the exact comparisons
  // is much cheaper -- the comparison of 2 or 3 vectors by element
  // is the same across all elements of the group.

  const size_t num_group = 1 << NUM_SHUFFLE;
  //const size_t group_size_max = (nfa+num_group-1) / num_group;
  const size_t group_size_max = gm_ceil_i8(nfa, num_group);

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        /*---Fill pad vectors with copies of the last vector---*/
        const size_t vector_capped = gm_min_i8(vector, nva-1);
        for (int fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          if (field >= nfa) {
            continue; // These entries will be padded to zero elsewhere.
          }
          const size_t f = field; // field number
          const size_t v = vector_capped; // vector number

          const size_t pf = perm(0, f, nfa); // permuted field number
          const size_t g = pf / group_size_max; // group number
          GMAssert(g>=0 && g<num_group);

          const size_t pv = perm(g, v, nva); // permuted vector number

          // Linearly map pv to small interval.
          const size_t value = value_min + (pv * value_max) / (value_min+nva);

          const GMFloat float_value = value;

          /*---Store---*/
          GMInsist(float_value * nfa >= 1);
          GMInsist(float_value * nfa < max_float);
          GMVectors_float_set(vectors, fl, vl, float_value, env);

        } /*---field_local---*/
      }   /*---vector_local---*/
      /*---Print---*/
//TODO: move this
      if (verbosity > 2) {
        GMVectors_print(vectors, env);
      }
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/

#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        /*---Fill pad vectors with copies of the last vector---*/
        const size_t vector_capped = gm_min_i8(vector, nva-1);
        for (int fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          if (field >= nfa) {
            continue; // These entries will be padded to zero elsewhere.
          }
          /*---Create 2-bit value - make extra sure less than 4---*/

          const size_t f = field;
          const size_t v = vector_capped;

          const size_t pf = perm(0, f, nfa);
          const size_t g = pf / group_size_max;
          GMAssert(g>=0 && g<num_group);

          const size_t pv = perm(g, v, nva);

          const size_t value = value_min + ( pv * value_max ) / (nva+value_min);

          const GMBits2 bval = ((size_t)3) & (value - value_min);
//mycount[bval]++;

          /*---Store---*/
          GMVectors_bits2_set(vectors, fl, vl, bval, env);

        } /*---field_local---*/
      }   /*---vector_local---*/
//TODO: move this
      if (verbosity > 2) {
        GMVectors_print(vectors, env);
      }
//printf("%i %i %i %i\n", mycount[0], mycount[1], mycount[2], mycount[3]);
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
}

//=============================================================================

void set_vectors_synthetic(GMVectors* vectors, int problem_type, int verbosity,
                           GMEnv* env) {
  GMInsist(vectors && env);

  if (problem_type == GM_PROBLEM_TYPE_RANDOM) {
    set_vectors_random_(vectors, verbosity, env);
  } else if (problem_type == GM_PROBLEM_TYPE_ANALYTIC) {
    set_vectors_analytic_(vectors, verbosity, env);
  } else {
    GMInsist(false && "Invalid problem_type");
  }
}

//=============================================================================
/*---Check correctness of metrics, if possible---*/

void check_metrics_analytic_(GMMetrics* metrics, DriverOptions* do_,
                             GMEnv* env) {
  GMInsist(metrics && do_ && env);
  GMInsist(GM_PROBLEM_TYPE_ANALYTIC == do_->problem_type);
  GMInsist(NULL == do_->input_file_path);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  const size_t nfa = metrics->num_field_active;
  const size_t nva = metrics->num_vector_active;

  // Upper bound on integer representable exactly by floating point type.
  // Account for cast to float in magma Volta version.
  const size_t max_float = ((size_t)1) <<
    (GMEnv_data_type_vectors(env) == GM_DATA_TYPE_FLOAT ?
     gm_mant_dig<float>() : gm_mant_dig<GMFloat>());
  // Czek account for number of terms summed in denom or num
  const size_t overflow_limit =
    GMEnv_data_type_vectors(env) != GM_DATA_TYPE_FLOAT ? 1 :
    GMEnv_num_way(env) == GM_NUM_WAY_2 ? 2 : 4;
  // Sum nfa times down the vector, is it still exact.
  const size_t value_limit = (max_float - 1) / (overflow_limit * nfa);

  const size_t value_min = 1;
  //const size_t value_max = (nva+value_min) < value_limit ?
  //                         (nva+value_min) : value_limit;
  const size_t value_max = gm_min_i8(value_min+nva, value_limit);

  const size_t num_group = 1 << NUM_SHUFFLE;
  //const size_t group_size_max = (nfa+num_group-1) / num_group;
  const size_t group_size_max = gm_ceil_i8(nfa, num_group);

  size_t num_incorrect = 0;
  const size_t max_to_print = 10;
  double max_incorrect_diff = 0.;

  switch (GMEnv_data_type_metrics(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
//      if (gm_gpu_compute_capability() == 700 &&
//          GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU &&
//          GM_FP_PRECISION_DOUBLE) {
//        // For this case modified MAGMA code casts down to single.
//        break;
//      }
      if (GMEnv_num_way(env) == GM_NUM_WAY_2) {
#pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t vi =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t vj =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          if (vi >= nva || vj >= nva) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czek_get_from_index(metrics, index, env);

          GMFloat float_n = 0;
          GMFloat float_d = 0;

          size_t n = 0;
          size_t d = 0;

          // For each comparison of vectors, the compared/summed
          // elements are treated as num_group groups.  All element
          // comparisons in the group have the same value, so we just
          // compute once and multiply that by the group size.

          for (size_t g=0; g<num_group; ++g) {

            const size_t pf_min = g * group_size_max;
            const size_t pf_max = gm_min_i8((g+1) * group_size_max, nfa);
            const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

            const size_t pvi = perm(g, vi, nva);
            const size_t pvj = perm(g, vj, nva);

            const size_t value_i = value_min + ( pvi * value_max ) /
                                               (value_min+nva);
            const size_t value_j = value_min + ( pvj * value_max ) /
                                               (value_min+nva);
            float_n += gm_min_i8(value_i, value_j) * gs_this;
            float_d += (value_i + value_j) * gs_this;
            n += gm_min_i8(value_i, value_j) * gs_this;
            d += (value_i + value_j) * gs_this;

          } //---g

          GMInsist(n == (size_t)float_n);
          GMInsist(d == (size_t)float_d);

          const GMFloat multiplier = (GMFloat)2;

          const GMFloat value_expected = (multiplier * float_n) / float_d;

          const bool is_incorrect = value_expected != value;
          if (is_incorrect) {
            const double diff = fabs(value - value_expected);
            max_incorrect_diff = diff > max_incorrect_diff ? diff : max_incorrect_diff;
            if (num_incorrect < max_to_print) {
              printf("Error: incorrect result detected.  coords %zu %zu  "
                     "expected %.20e  actual %.20e  diff %.20e\n", vi, vj,
                     (double)value_expected, (double)value,
                     (double)value-(double)value_expected);
            }
          }

          num_incorrect += is_incorrect;
        } //---for index
      } //---if
      if (GMEnv_num_way(env) == GM_NUM_WAY_3) {
#pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t vi =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t vj =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          const size_t vk =
            GMMetrics_coord_global_from_index(metrics, index, 2, env);
          if (vi >= nva || vj >= nva || vk >= nva) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czek_get_from_index(metrics, index, env);

          GMFloat float_n = 0;
          GMFloat float_d = 0;

          for (size_t g=0; g<num_group; ++g) {

            const size_t pf_min = g * group_size_max;
            const size_t pf_max = gm_min_i8((g+1) * group_size_max, nfa);
            const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

            const size_t pvi = perm(g, vi, nva);
            const size_t pvj = perm(g, vj, nva);
            const size_t pvk = perm(g, vk, nva);

            const size_t value_i = value_min + ( pvi * value_max ) /
                                               (nva+value_min);
            const size_t value_j = value_min + ( pvj * value_max ) /
                                               (nva+value_min);
            const size_t value_k = value_min + ( pvk * value_max ) /
                                               (nva+value_min);

            float_n += gm_min_i8(value_i, value_j) * gs_this;
            float_n += gm_min_i8(value_i, value_k) * gs_this;
            float_n += gm_min_i8(value_j, value_k) * gs_this;

            float_n -= gm_min_i8(value_i, gm_min_i8(value_j, value_k)) * gs_this;

            float_d += (value_i + value_j + value_k) * gs_this;

          } //---g

          const GMFloat multiplier = (GMFloat)1.5;

          const GMFloat value_expected = (multiplier * float_n) / float_d;

          const bool is_incorrect = value_expected != value;
          if (is_incorrect) {
            const double diff = fabs(value - value_expected);
            max_incorrect_diff = diff > max_incorrect_diff ? diff : max_incorrect_diff;
            if (num_incorrect < max_to_print) {
              printf("Error: incorrect result detected.  coords %zu %zu %zu  "
                     "expected %.20e  actual %.20e  diff %.20e\n", vi, vj, vk,
                     (double)value_expected, (double)value,
                     (double)value-(double)value_expected);
            }
          }

          num_incorrect += is_incorrect;
        } //---for index
      } //---if
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
    /*--------------------*/

    const int cbpe = GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC ? 2 : 1;

#pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
      for (size_t index = 0; index < metrics->num_elts_local; ++index) {
        const size_t vi =
          GMMetrics_coord_global_from_index(metrics, index, 0, env);
        const size_t vj =
          GMMetrics_coord_global_from_index(metrics, index, 1, env);
        if (vi >= nva || vj >= nva) {
          continue;
        }
        for (int i0 = 0; i0 < 2; ++i0) {
          for (int i1 = 0; i1 < 2; ++i1) {
            const GMFloat value = cbpe == 2 ?
                GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env) :
                GMMetrics_duo_get_from_index_2(metrics, index, i0, i1, env);

            GMTally1 rij = 0;
            GMTally1 si = 0;
            GMTally1 sj = 0;
            GMTally1 ci = 0;
            GMTally1 cj = 0;
            GMTally1 cij = 0;

            for (size_t g=0; g<num_group; ++g) {

              const size_t pf_min = g * group_size_max;
              const size_t pf_max = gm_min_i8((g+1) * group_size_max, nfa);
              const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

              const size_t pvi = perm(g, vi, nva);
              const size_t pvj = perm(g, vj, nva);

              const size_t value_i = value_min + ( pvi * value_max ) /
                                                 (nva+value_min);
              const size_t value_j = value_min + ( pvj * value_max ) /
                                                 (nva+value_min);

              const GMBits2 bval_i = ((size_t)3) & (value_i - value_min);
              const GMBits2 bval_j = ((size_t)3) & (value_j - value_min);

              const int bval_i_0 = !!(bval_i&1);
              const int bval_i_1 = !!(bval_i&2);
              const int bval_j_0 = !!(bval_j&1);
              const int bval_j_1 = !!(bval_j&2);

              const bool unknown_i = env->sparse && bval_i == GM_2BIT_UNKNOWN;
              const bool unknown_j = env->sparse && bval_j == GM_2BIT_UNKNOWN;
              const bool unknown_ij = unknown_i || unknown_j;

              if (! unknown_i) {
                ci += gs_this;
                si += cbpe == 2 ?
                  ((bval_i_0 == i0) + (bval_i_1 == i0)) * gs_this :
                  (bval_i_0 == i0) * gs_this;
              }

              if (! unknown_j) {
                cj += gs_this;
                sj += cbpe == 2 ?
                  ((bval_j_0 == i1) + (bval_j_1 == i1)) * gs_this :
                  (bval_j_0 == i1) * gs_this;
              }

              if (! unknown_ij) {
                cij += cbpe * cbpe * gs_this;
                rij += cbpe == 2 ?
                       (((bval_i_0 == i0) && (bval_j_0 == i1)) +
                        ((bval_i_0 == i0) && (bval_j_1 == i1)) +
                        ((bval_i_1 == i0) && (bval_j_0 == i1)) +
                        ((bval_i_1 == i0) && (bval_j_1 == i1))) *
                       gs_this :
                       ((bval_i_0 == i0) && (bval_j_0 == i1)) *
                       gs_this;
              }
            } //---g

//printf("%i %i %i\n", (int) i0, (int)i1, (int)rij);

            GMFloat value_expected_floatcalc = 0;
            if (!(ci == 0 || cj == 0 || cij == 0)) {
              const GMFloat f_one = 1;

              const GMFloat f_ci = (GMFloat) ci;
              const GMFloat f_cj = (GMFloat) cj;

              const GMFloat f_cicj_min = f_ci < f_cj ? f_ci : f_cj;
              const GMFloat f_cicj_max = f_ci > f_cj ? f_ci : f_cj;

              const GMFloat f_cij = (GMFloat) cij;
              const GMFloat recip_cicjcij = f_one /
                                            (f_cicj_min * f_cicj_max * f_cij);

              const GMFloat recip_ci = env->sparse ?
                f_cj * f_cij * recip_cicjcij : metrics->recip_m;
              const GMFloat recip_cj = env->sparse ?
                f_ci * f_cij * recip_cicjcij : metrics->recip_m;

              const GMFloat recip_sumcij = env->sparse ?
                f_cicj_min * f_cicj_max * recip_cicjcij :
                (f_one / (cbpe * cbpe)) * metrics->recip_m;

              //const GMFloat recip_ci = env->sparse ? f_one/ci : metrics->recip_m;
              //const GMFloat recip_cj = env->sparse ? f_one/cj : metrics->recip_m;

              //const GMFloat recip_sumcij = env->sparse ? f_one/cij :
              //                               (f_one / 4) * metrics->recip_m;

              value_expected_floatcalc = cbpe == 2 ?
                GMMetrics_ccc_duo_value_2<2>(metrics, rij, si, sj,
                                    recip_ci, recip_cj, recip_sumcij, env) :
                GMMetrics_ccc_duo_value_2<1>(metrics, rij, si, sj,
                                    recip_ci, recip_cj, recip_sumcij, env);
            }

            GMFloat value_expected = value_expected_floatcalc;

#if 0
//#ifdef HAVE_INT128
            if (env->are_ccc_params_default) {
            if (!(0 == ci || 0 == cj || 0 == cij)) {
              value_expected = GMMetrics_ccc_value_nofp_2(metrics,
                rij, si, sj, ci, cj, cij, env); 
            }
            }
#endif

            const bool is_incorrect = value_expected != value;
            if (is_incorrect) {
              const double diff = fabs(value - value_expected);
              max_incorrect_diff = diff > max_incorrect_diff ?
                                   diff : max_incorrect_diff;
              if (num_incorrect < max_to_print) {
                printf("Error: incorrect result detected.  coords %zu %zu  "
                       "expected %.20e  actual %.20e  diff %.20e\n", vi, vj,
                       (double)value_expected, (double)value,
                       (double)value-(double)value_expected);
              }
            }

            num_incorrect += is_incorrect;
          } //---j
        } //---i
      } //---for index
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY4X2: {
    /*--------------------*/
    GMInsist(GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC);
#pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
      for (size_t index = 0; index < metrics->num_elts_local; ++index) {
        const size_t vi =
          GMMetrics_coord_global_from_index(metrics, index, 0, env);
        const size_t vj =
          GMMetrics_coord_global_from_index(metrics, index, 1, env);
        const size_t vk =
          GMMetrics_coord_global_from_index(metrics, index, 2, env);
        if (vi >= nva || vj >= nva || vk >= nva) {
          continue;
        }
        for (int i0 = 0; i0 < 2; ++i0) {
          for (int i1 = 0; i1 < 2; ++i1) {
            for (int i2 = 0; i2 < 2; ++i2) {
              const GMFloat value =
               GMMetrics_ccc_get_from_index_3( metrics, index, i0, i1, i2, env);

              GMTally1 rijk = 0;
              GMTally1 si = 0;
              GMTally1 sj = 0;
              GMTally1 sk = 0;
              GMTally1 ci = 0;
              GMTally1 cj = 0;
              GMTally1 ck = 0;
              GMTally1 cijk = 0;

              for (size_t g=0; g<num_group; ++g) {

                const size_t pf_min = g * group_size_max;
                const size_t pf_max = gm_min_i8((g+1) * group_size_max, nfa);
                const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

                const size_t pvi = perm(g, vi, nva);
                const size_t pvj = perm(g, vj, nva);
                const size_t pvk = perm(g, vk, nva);

                const size_t value_i = value_min + ( pvi * value_max ) /
                                                   (nva+value_min);
                const size_t value_j = value_min + ( pvj * value_max ) /
                                                   (nva+value_min);
                const size_t value_k = value_min + ( pvk * value_max ) /
                                                   (nva+value_min);

                const GMBits2 bval_i = ((size_t)3) & (value_i - value_min);
                const GMBits2 bval_j = ((size_t)3) & (value_j - value_min);
                const GMBits2 bval_k = ((size_t)3) & (value_k - value_min);

                const int bval_i_0 = !!(bval_i&1);
                const int bval_i_1 = !!(bval_i&2);
                const int bval_j_0 = !!(bval_j&1);
                const int bval_j_1 = !!(bval_j&2);
                const int bval_k_0 = !!(bval_k&1);
                const int bval_k_1 = !!(bval_k&2);


                const bool unknown_i = env->sparse && bval_i == GM_2BIT_UNKNOWN;
                const bool unknown_j = env->sparse && bval_j == GM_2BIT_UNKNOWN;
                const bool unknown_k = env->sparse && bval_k == GM_2BIT_UNKNOWN;
                const bool unknown_ijk = unknown_i || unknown_j || unknown_k;

                if (! unknown_i) {
                  ci += gs_this;
                  si += ((bval_i_0 == i0) + (bval_i_1 == i0)) * gs_this;
                }

                if (! unknown_j) {
                  cj += gs_this;
                  sj += ((bval_j_0 == i1) + (bval_j_1 == i1)) * gs_this;
                }

                if (! unknown_k) {
                  ck += gs_this;
                  sk += ((bval_k_0 == i2) + (bval_k_1 == i2)) * gs_this;
                }

                if (! unknown_ijk) {
                  cijk += 8 * gs_this;
                  rijk += (((bval_i_0==i0) && (bval_j_0==i1) && (bval_k_0==i2))+
                           ((bval_i_1==i0) && (bval_j_0==i1) && (bval_k_0==i2))+
                           ((bval_i_0==i0) && (bval_j_1==i1) && (bval_k_0==i2))+
                           ((bval_i_1==i0) && (bval_j_1==i1) && (bval_k_0==i2))+
                           ((bval_i_0==i0) && (bval_j_0==i1) && (bval_k_1==i2))+
                           ((bval_i_1==i0) && (bval_j_0==i1) && (bval_k_1==i2))+
                           ((bval_i_0==i0) && (bval_j_1==i1) && (bval_k_1==i2))+
                           ((bval_i_1==i0) && (bval_j_1==i1) && (bval_k_1==i2)))
                          * gs_this;
                }
              } //---g

              GMFloat value_expected_floatcalc = 0;
              if (!(ci == 0 || cj == 0 || ck == 0 || cijk == 0)) {
                const GMFloat f_one = 1;
  
                const GMFloat recip_ci = env->sparse ? f_one/ci
                                                     : metrics->recip_m;
                const GMFloat recip_cj = env->sparse ? f_one/cj
                                                     : metrics->recip_m;
                const GMFloat recip_ck = env->sparse ? f_one/ck
                                                     : metrics->recip_m;
  
                const GMFloat recip_sumcijk = env->sparse ? f_one/cijk :
                                               (f_one / 8) * metrics->recip_m;
  
                value_expected_floatcalc =
                  GMMetrics_ccc_value_3(metrics, rijk, si, sj, sk, recip_ci,
                                        recip_cj, recip_ck, recip_sumcijk, env);
              }

              GMFloat value_expected = value_expected_floatcalc;

#if 0
//#ifdef HAVE_INT128
              if (env->are_ccc_params_default) {
              if (!(0 == ci || 0 == cj || 0 == ck || 0 == cijk)) {
                value_expected = GMMetrics_ccc_value_nofp_3(metrics,
                  rijk, si, sj, sk, ci, cj, ck, cijk, env); 
              }
              }
#endif

              const bool is_incorrect = value_expected != value;
              if (is_incorrect) {
                const double diff = fabs(value - value_expected);
                max_incorrect_diff = diff > max_incorrect_diff ? diff : max_incorrect_diff;
                if (num_incorrect < max_to_print) {
                  printf("Error: incorrect result detected.  coords %zu %zu %zu  "
                         "expected %.20e  actual %.20e  diff %.20e\n", vi, vj, vk,
                         (double)value_expected, (double)value,
                         (double)value-(double)value_expected);
                }
              }

              num_incorrect += is_incorrect;
            } //---k
          } //---j
        } //---i
      } //---for index
    } break;
    /*--------------------*/
    default:
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
  do_->num_incorrect += num_incorrect;
  do_->max_incorrect_diff = max_incorrect_diff > do_->max_incorrect_diff ?
                            max_incorrect_diff : do_->max_incorrect_diff;
}

//=============================================================================

void check_metrics(GMMetrics* metrics, DriverOptions* do_, GMEnv* env) {
  GMInsist(metrics && do_ && env);

  if (NULL != do_->input_file_path) {
    return;
  }

  if (GM_PROBLEM_TYPE_ANALYTIC == do_->problem_type) {
    check_metrics_analytic_(metrics, do_, env);
  }
}

//=============================================================================
