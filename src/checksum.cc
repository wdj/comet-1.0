//-----------------------------------------------------------------------------
/*!
 * \file   checksum.cc
 * \author Wayne Joubert
 * \date   Mon Aug  7 14:47:01 EDT 2017
 * \brief  Utility to compute checksum of data in a metrics object.
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

#include "cstdlib"
#include "cstdint"
#include "string.h"
#include "math.h"

#include "mpi.h"

#include "env.hh"
#include "metrics.hh"
#include "checksum.hh"

//-----------------------------------------------------------------------------

namespace CoMet {

//=============================================================================
/// \brief Checksum default constructor.

Checksum::Checksum(bool computing_checksum)
  : is_overflowed_(false)
  , value_max_(-DBL_MAX)
  , sum_d_(0)
  , is_started_(false)
  , computing_checksum_(computing_checksum) {

  for (int i=0; i<SIZE; ++i) {
    data_[i] = 0;
  }
}

//-----------------------------------------------------------------------------
/// \brief Manual copy of checksum entries.

void Checksum::copy(const Checksum& cksum) {
  for (int i=0; i<SIZE; ++i) {
    data_[i] = cksum.data_[i];
  }

  is_overflowed_ = cksum.is_overflowed_;
  value_max_ = cksum.value_max_;
  sum_ = cksum.sum_;
  sum_d_ = cksum.sum_d_;
  is_started_ = cksum.is_started_;
  computing_checksum_ = cksum.computing_checksum_;
}

//-----------------------------------------------------------------------------
/// \brief Check whether two checksums are equal.

bool Checksum::is_equal(const Checksum& cksum2) const {
  bool result = true;

  // Don't perform this test if not computing both checksums.
  if (this->computing_checksum_ && cksum2.computing_checksum_) {
    for (int i = 0; i < SIZE; ++i) {
      result = result && this->data_[i] == cksum2.data_[i];
    }
    // ISSUE: for now, check overflow like this because
    // overflow is not necessarily an error (e.g., setting of
    // ccc_multiplier).
    // TODO: fix this better.
    result = result && this->is_overflowed_ == cksum2.is_overflowed_;
    result = result && this->value_max_ == cksum2.value_max_;
  }

  return result;
}

//-----------------------------------------------------------------------------
/// \brief Print ckecksum to stdout.

void Checksum::print(GMEnv& env) {

  for (int i = 0; i < SIZE; ++i) {
    printf("%s%li", i == 0 ? "" : "-",
           this->data_[SIZE - 1 - i]);
  }
  if (this->is_overflowed_) {
    printf("-OVFL");
    printf("-%e", this->value_max_);
  }
}

//-----------------------------------------------------------------------------
/// \brief Checksum helper function: perform one bubble sort step.

inline static void makegreater(size_t& i, size_t& j, int& ind_i, int& ind_j) {
  if (i < j) {
    const size_t tmp = i;
    i = j;
    j = tmp;
    const int tmp2 = ind_i;
    ind_i = ind_j;
    ind_j = tmp2;
  }
}

//-----------------------------------------------------------------------------
/// \brief Checksum helper function: left shift that works for any shift amount.

inline static size_t lshift(size_t a, int j) {
  if (j >= 64 || j <= -64) {
    return 0;
  }
  return j > 0 ? a << j : a >> (-j);
}

//-----------------------------------------------------------------------------
/// \brief Checksum helper: return largest value in metrics object.
///
///        Note: this computes the max value on proc, not across procs,
///        so each proc's result in general can have a different value.

double Checksum::metrics_max_value(GMMetrics& metrics, GMEnv& env) {

  double result = -DBL_MAX;

  // TODO: make this unnecessary.
  if (! GMEnv_is_proc_active(&env)) {
    return result;
  }

  // Loop over metrics indices to find max.
  #pragma omp parallel for reduction(max:result)
  for (size_t index = 0; index < metrics.num_elts_local; ++index) {
    // Determine whether this cell is active.
    bool is_active = true;
    for (int i = 0; i < GMEnv_num_way(&env); ++i) {
      const size_t coord = GMMetrics_coord_global_from_index(&metrics, index,
                                                             i, &env);
      is_active = is_active && coord < metrics.num_vector_active;
    }
    double value_max = -DBL_MAX;
    if (is_active) {
      // Loop over data values at this index
      for (int i_value = 0; i_value < metrics.data_type_num_values; ++i_value) {
        // Pick up value of this metrics elt
        double value = 0;
        switch (metrics.data_type_id) {
          // --------------
          case GM_DATA_TYPE_FLOAT: {
            value = GMMetrics_czek_get_from_index(&metrics, index, &env);
          } break;
          // --------------
          case GM_DATA_TYPE_TALLY2X2: {
            const int i0 = i_value / 2;
            const int i1 = i_value % 2;
            value = GMEnv_metric_type(&env) == GM_METRIC_TYPE_CCC ?
              GMMetrics_ccc_get_from_index_2(&metrics, index, i0, i1, &env) :
              GMMetrics_duo_get_from_index_2(&metrics, index, i0, i1, &env);
          } break;
          // --------------
          case GM_DATA_TYPE_TALLY4X2: {
            const int i0 = i_value / 4;
            const int i1 = (i_value / 2) % 2;
            const int i2 = i_value % 2;
            value =
              GMMetrics_ccc_get_from_index_3(&metrics, index, i0, i1, i2, &env);
          } break;
          // --------------
          default:
          GMInsist(false && "Invalid metrics data type. metrics.data_type_id.");
        } // switch
        // value_max is the largest of the values at this index.
        value_max = value > value_max ? value : value_max;
      } // for i_value
    } // if is_active
    result = value_max > result ? value_max : result;
  } // for index

  return result;
} // Checksum::metrics_max_value

//-----------------------------------------------------------------------------
/// \brief compute (global and local) checksum of metrics object.
///
///        The global checksum (across all procs) is generally of most
///        interest.  The local checksum is to debug difficult cases
///        when checksum errors need to be narrowed down to the
///        specific proc.
///
///        Note that cksum and cksum_local are input/output
///        variables; they are added to for multiple stages or phases
///        of the calculation.

void Checksum::compute(Checksum& cksum, Checksum& cksum_local,
                       GMMetrics& metrics, GMEnv& env){
  // TODO: make this check unnecessary.
  GMInsist(metrics.data || ! GMEnv_is_proc_active(&env));

  // TODO: make this unnecessary.
  if (! GMEnv_is_proc_active(&env)) {
    return;
  }

  enum { NUM_WAY_MAX = GM_NUM_NUM_WAY + 1 };
  GMInsist(GMEnv_num_way(&env) <= NUM_WAY_MAX && "This num_way not supported.");

  // Check for NaNs if appropriate

  // TODO: put this in metrics class - a heavyweight validity check function
  switch (metrics.data_type_id) {
    case GM_DATA_TYPE_FLOAT: {
      GMFloat_check((GMFloat*)(metrics.data), metrics.num_elts_local);
    } break;
  }

  // Get max metrics value.

  // The only impact this has on the output is to determine whether
  // an "overflow" in size of the metric value beyond what the checksum
  // functionality can compute has occurred.
  // It might also be useful for debugging.
  // Note there is no effort to strictly differentiate the
  // local (per proc) vs. global value.

  double value_max_tmp = Checksum::metrics_max_value(metrics, env);
  value_max_tmp = value_max_tmp > cksum.value_max_ ?
                  value_max_tmp : cksum.value_max_;

  int mpi_code = MPI_Allreduce(&value_max_tmp, &cksum.value_max_, 1,
                        MPI_DOUBLE, MPI_MAX, GMEnv_mpi_comm_repl_vector(&env));
  GMInsist(mpi_code == MPI_SUCCESS && "Faiure in call to MPI_Allreduce.");
  cksum_local.value_max_ = cksum.value_max_;

  // Check whether values are within a range for which we can compute
  // the checksum with this code.

  // The largest we expect any value to be if using "special" inputs.
  //
  // This size constraint may be violated if there is an error in the
  // calculation.  It may also be violated by setting of ccc_multiplier large.
  // TODO: examine whether the bound here could be made tighter.
  // TODO: consider polling the metrics object for what the max value
  // upper bound should be - e.g., ccc_multiplier * (1 + roundoff_fuzz)
  const int log2_value_max_allowed = 4;
  const double value_max_allowed = 1 << log2_value_max_allowed;

  cksum.is_overflowed_ = cksum_local.is_overflowed_ =
    cksum.is_overflowed_ && cksum.value_max_ > value_max_allowed;

  // Scaling factor for values - so that after scaling, value is <= 1.
  const double scaling = value_max_allowed;

  //--------------------
  // Calculate checksum
  //--------------------

  typedef uint64_t UI64_t;

  const int w = 30; // 2*w is the integer size
  GMInsist(64 - 2 * w >= 4); // fits into uint64, with some headroom
  const UI64_t one64 = 1; // the constant "1"

  const UI64_t lomask = (one64 << w) - 1; // masks for lo and hi parts of int
  const UI64_t lohimask = (one64 << (2 * w)) - 1;

  MultiprecInt sum_local; // = 0 // checksum valune on this proc
  double sum_d_local = 0; // floating point representation of the same, as check

  #pragma omp parallel
  {
    MultiprecInt sum_local_private; // = 0
    double sum_d_local_private = 0;
    // Loop over metrics indices to get checksum contribution.
    #pragma omp for collapse(2)
    for (size_t index = 0; index < metrics.num_elts_local; ++index) {
      // Loop over data values at this index
      for (int i_value = 0; i_value < metrics.data_type_num_values; ++i_value) {
        // Obtain global coords of metrics elt
        size_t coords[NUM_WAY_MAX];
        int ind_coords[NUM_WAY_MAX]; // permutation index
        for (int i = 0; i < NUM_WAY_MAX; ++i) {
          coords[i] = 0;
          ind_coords[i] = i;
        }
        bool is_active = true;
        for (int i = 0; i < GMEnv_num_way(&env); ++i) {
          const size_t coord =
            GMMetrics_coord_global_from_index(&metrics, index, i, &env);
          // Ignore padding vectors.
          is_active = is_active && coord < metrics.num_vector_active;
          coords[i] = coord;
        }
        // Reflect coords by symmetry to get uniform result -
        //   sort into descending order
        //
        // The idea here is that, because the tensor has reflective
        // symmetries, a different equivlent reflected tensor value may be
        // computed based on the parallel decomposition.
        // This permutation puts the indices into a uniform order
        // so that this is not viewed as a difference in the results.
        // Note also below we will permute i0 / i1 / i2 as needed.
        makegreater(coords[1], coords[2], ind_coords[1], ind_coords[2]);
        makegreater(coords[0], coords[1], ind_coords[0], ind_coords[1]);
        makegreater(coords[1], coords[2], ind_coords[1], ind_coords[2]);

        // Pick up value of this metrics elt
        double value = 0;
        switch (metrics.data_type_id) {
          // --------------
          case GM_DATA_TYPE_FLOAT: {
            value = GMMetrics_czek_get_from_index(&metrics, index, &env);
          } break;
          // --------------
          case GM_DATA_TYPE_TALLY2X2: {
            const int i0_unpermuted = i_value / 2;
            const int i1_unpermuted = i_value % 2;
            const int i0 = ind_coords[0] == 0 ? i0_unpermuted : i1_unpermuted;
            const int i1 = ind_coords[0] == 0 ? i1_unpermuted : i0_unpermuted;
            value = GMEnv_metric_type(&env) == GM_METRIC_TYPE_CCC ?
              GMMetrics_ccc_get_from_index_2(&metrics, index, i0, i1, &env) :
              GMMetrics_duo_get_from_index_2(&metrics, index, i0, i1, &env);
          } break;
          // --------------
          case GM_DATA_TYPE_TALLY4X2: {
            const int i0_unpermuted = i_value / 4;
            const int i1_unpermuted = (i_value / 2) % 2;
            const int i2_unpermuted = i_value % 2;
            const int i0 = ind_coords[0] == 0 ? i0_unpermuted :
                           ind_coords[1] == 0 ? i1_unpermuted :
                                                i2_unpermuted;
            const int i1 = ind_coords[0] == 1 ? i0_unpermuted :
                           ind_coords[1] == 1 ? i1_unpermuted :
                                                i2_unpermuted;
            const int i2 = ind_coords[0] == 2 ? i0_unpermuted :
                           ind_coords[1] == 2 ? i1_unpermuted :
                                                i2_unpermuted;
            value =
              GMMetrics_ccc_get_from_index_3(&metrics, index, i0, i1, i2, &env);
          } break;
          // --------------
          default:
            GMInsist(false && "Invalid data type. metrics.data_type_id.");
        } // switch

        // Convert to uint64.  Store only 2*w+1 bits, at most -
        // if (value / scaling) <= 1, which it should be if
        // floating point arithmetic works as expected,
        // must have ivalue <= (1<<(2*w)).
        // Note: it would have been better to set this to be 2*w bits max
        // instead - would need to subtract 1 from the second multiplicand
        // and make sure floating point arith works as expected.
        // HOWEVER, see note below.
        UI64_t ivalue = (UI64_t)( (value / scaling) * (one64 << (2 * w)) );
        // Construct an id that is a single number representing the coord
        // and value number.
        UI64_t uid = coords[0];
        for (int i = 1; i < GMEnv_num_way(&env); ++i) {
          uid = uid * metrics.num_vector_active + coords[i];
        }
        uid = uid * metrics.data_type_num_values + i_value;
        // Randomize this id
        const UI64_t rand1 = gm_randomize(uid + 956158765);
        const UI64_t rand2 = gm_randomize(uid + 842467637);
        UI64_t rand_value = rand1 + gm_randomize_max() * rand2;
        // Truncate to 2*w bits.
        rand_value &= lohimask;
        // Multiply the two values.
        const UI64_t a = rand_value;
        const UI64_t alo = a & lomask;
        const UI64_t ahi = a >> w;
        const UI64_t b = ivalue;
        const UI64_t blo = b & lomask;
        const UI64_t bhi = b >> w;
        const UI64_t cx = alo * bhi + ahi * blo;
        // Note: since a < (1<<(2*w)) and b <= (1<<(2*w)),
        // it is guaranteed that c < (1<<(4*w)),
        // so the result is in the right range of bits.
        UI64_t clo = alo * blo + ((cx & lomask) << w);
        UI64_t chi = ahi * bhi + (cx >> w);
        // (move the carry bits)
        chi += clo >> (2 * w);
        clo &= lohimask;
        // The checksum is the product of the metric value and the
        // global coord / value number id, this summed across all.
        const double value_d =
            ivalue * (double)rand_value / ((double)(one64 << (2 * w)));
        // Accumulate
        if (is_active) {
          sum_d_local_private += value_d; // (private) reduction
          // Split the product into one-char chunks, accumulate to sums
          for (int i = 0; i < 8; ++i) {
            const UI64_t value0 = (clo << (64 - 8 - 8 * i)) >> (64 - 8);
            const UI64_t value1 = (chi << (64 - 8 - 8 * i)) >> (64 - 8);
            sum_local_private.data_[0 + i] += value0; // (private) reduction
            sum_local_private.data_[8 + i] += value1; // (private) reduction
          }
        }
      } // for i_value
    } // for index
    // omp for collapse

    #pragma omp critical
    {
        sum_d_local += sum_d_local_private; // critical reduction
        for (int i = 0; i < 8; ++i) {
          sum_local.data_[0 + i] += sum_local_private.data_[0 + i]; // critical
          sum_local.data_[8 + i] += sum_local_private.data_[8 + i]; // reduction
        }
    }
  } // omp parallel

  // Global sum of multiprecision int

  MultiprecInt sum; // = 0
  mpi_code = MPI_Allreduce(sum_local.data_, sum.data_,
                           MultiprecInt::SIZE,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(&env));
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Allreduce.");
  // Add to multiprec int data we have so far.
  for (int i = 0; i < MultiprecInt::SIZE; ++i) {
    cksum.sum_.data_[i] += sum.data_[i];
    cksum_local.sum_.data_[i] += sum_local.data_[i];
  }

  // Condense results into smaller number of 64 bit ints.

  for (int i = 0; i < SIZE; ++i) {
    cksum.data_[i] = 0;
    cksum_local.data_[i] = 0;
    for (int j = 0; j < 8; ++j) {
      cksum.data_[i] +=
        lshift(cksum.sum_.data_[0 + j], 8*j - 2*w*i) & lohimask;
      cksum.data_[i] +=
        lshift(cksum.sum_.data_[8 + j], 8*j - 2*w*(i - 1)) & lohimask;
      cksum_local.data_[i] +=
        lshift(cksum_local.sum_.data_[0 + j], 8*j - 2*w*i) & lohimask;
      cksum_local.data_[i] +=
        lshift(cksum_local.sum_.data_[8 + j], 8*j - 2*w*(i - 1)) & lohimask;
    }
  }
  // Adjustments: move the carry bits
  cksum.data_[1] += cksum.data_[0] >> (2 * w);
  cksum.data_[0] &= lohimask;
  cksum.data_[2] += cksum.data_[1] >> (2 * w);
  cksum.data_[1] &= lohimask;
  cksum_local.data_[1] += cksum_local.data_[0] >> (2 * w);
  cksum_local.data_[0] &= lohimask;
  cksum_local.data_[2] += cksum_local.data_[1] >> (2 * w);
  cksum_local.data_[1] &= lohimask;

  // Validate by checking against floating point result

  double sum_d = 0;
  mpi_code = MPI_Allreduce(&sum_d_local, &sum_d, 1, MPI_DOUBLE, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(&env));
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Allreduce.");
  cksum.sum_d_ += sum_d;
  cksum_local.sum_d_ += sum_d_local;

  double result_d = cksum.data_[0] / ((double)(one64 << (2 * w))) +
                    cksum.data_[1] +
                    cksum.data_[2] * ((double)(one64 << (2 * w)));
  GMInsist(fabs(cksum.sum_d_ - result_d) <= cksum.sum_d_ * 1.e-10 &&
           "Error in checksum calculation");
} // Checksum::compute

//=============================================================================

} // namespace CoMet

//-----------------------------------------------------------------------------
