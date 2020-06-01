//-----------------------------------------------------------------------------
/*!
 * \file   decomp_mgr.hh
 * \author Wayne Joubert
 * \date   Tue Aug  8 19:58:57 EDT 2017
 * \brief  Define distribution of vectors to MPI ranks, padding needed, etc.
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
#include "string.h"

#include "mpi.h"

#include "env.hh"
#include "linalg_accel.hh"
#include "decomp_mgr.hh"

//=============================================================================

size_t gm_num_vector_local_required(size_t num_vector_active_local,
                                    GMEnv* const env) {
  GMInsist(env);
  // NOTE: this function should receive the same num_vector_active_local
  // and give the same result independent of MPI rank.

  const size_t factor_4 = gm_gemm_divisibility_required(env);

  const bool need_divisible_by_6 = GMEnv_num_way(env) == GM_NUM_WAY_3 &&
                                   GMEnv_all2all(env) &&
                                   GMEnv_num_proc_vector_i(env) > 2;

  const size_t lcm = (! need_divisible_by_6) ? factor_4 :
                     factor_4 % 2 == 0 ? 3 * factor_4 : 6 * factor_4;

  return gm_ceil_i8(num_vector_active_local, lcm)*lcm;
}

//-----------------------------------------------------------------------------
// Set to null

GMDecompMgr GMDecompMgr_null() {
  GMDecompMgr result;
  memset((void*)&result, 0, sizeof(result));
  return result;
}

//-----------------------------------------------------------------------------
// (Pseudo) constructor

void GMDecompMgr_create(GMDecompMgr* dm,
                        bool fields_by_local,
                        bool vectors_by_local,
                        size_t num_field_specifier,
                        size_t num_vector_specifier,
                        int vectors_data_type_id,
                        GMEnv* env) {
  GMInsist(dm && env);

  if (! GMEnv_is_proc_active(env)) {
    *dm = GMDecompMgr_null();
    return;
  }

  //const int vectors_data_type_id = GMEnv_data_type_vectors(env);

  //--------------------
  // Vector counts
  //--------------------

  if (vectors_by_local) {
    dm->num_vector_local = num_vector_specifier;
    const size_t num_vector_local_required = gm_num_vector_local_required(
                                              dm->num_vector_local, env);
    GMInsistInterface(env, dm->num_vector_local == num_vector_local_required &&
         "Manual selection of nvl requires divisibility condition");
    // All vectors active on every proc.
    dm->num_vector_active_local = dm->num_vector_local;
#if 0
    dm->num_vector_active_local = num_vector_specifier;
    dm->num_vector_local = gm_num_vector_local_required(
                               dm->num_vector_active_local, env);
#endif
    dm->num_vector_active = dm->num_vector_active_local *
                            GMEnv_num_proc_vector_i(env);
    dm->num_vector = dm->num_vector_local * GMEnv_num_proc_vector_i(env);
  } else { // ! vectors_by_local
    dm->num_vector_active = num_vector_specifier;
    // Pad up as needed, require every proc has same number
    const int num_proc = GMEnv_num_proc_vector_i(env);
    const int proc_num = GMEnv_proc_num_vector_i(env);
    //dm->num_vector_local = gm_ceil_i8(dm->num_vector_active, num_proc);
    dm->num_vector_local = gm_num_vector_local_required(
      gm_ceil_i8(dm->num_vector_active, GMEnv_num_proc_vector_i(env)), env);
    dm->num_vector = dm->num_vector_local * num_proc;
    // Lower procs fully packed with active values
    // Upper procs fully inactive
    // One proc in between may be mixed
    const size_t nvl = dm->num_vector_local;
    const size_t nva = dm->num_vector_active;
    // Pack in lower procs with no gaps to ensure indexing of actives works
    // right independent of decomposition
    dm->num_vector_active_local = nva <= nvl * proc_num ? 0 :
                                  nva >= nvl * (proc_num + 1) ? nvl :
                                  nva - nvl * proc_num;
  } // if vectors_by_local

  //--------------------
  // Check the sizes
  //--------------------

  int mpi_code = 0;
  size_t sum = 0;

#if 0
printf("%i %i %i %i %i\n",
(int)GMEnv_num_proc_vector_i(env),
(int)dm->num_vector_active,
(int)dm->num_vector,
(int)dm->num_vector_active_local,
(int)dm->num_vector_local
);
#endif

  GMInsist(dm->num_vector_active >= 0 &&
           dm->num_vector_active <= dm->num_vector);
  GMInsist(dm->num_vector_active_local >= 0 &&
           dm->num_vector_active_local <= dm->num_vector_local);

  mpi_code = MPI_Allreduce(&dm->num_vector_local, &sum, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Allreduce.");
  GMInsist(sum == dm->num_vector_local * GMEnv_num_proc_repl_vector(env) &&
           "Every process must have the same number of vectors.");
  GMInsist(sum == dm->num_vector * GMEnv_num_proc_repl(env) &&
           "Error in local/global sizes computation.");

  mpi_code = MPI_Allreduce(&dm->num_vector_active_local, &sum, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Allreduce.");
  GMInsist(sum == dm->num_vector_active * GMEnv_num_proc_repl(env) &&
           "Error in local/global sizes computation.");

  //--------------------
  // Field counts
  //--------------------

  if (fields_by_local) {
    dm->num_field_active_local = num_field_specifier;
    dm->num_field_local = dm->num_field_active_local;
    dm->num_field = dm->num_field_local * GMEnv_num_proc_field(env);
    dm->num_field_active = dm->num_field;
    dm->field_base = dm->num_field_local * GMEnv_proc_num_field(env);
  } else { // ! fields_by_local
    dm->num_field_active = num_field_specifier;
    // Pad up as needed so that every proc has same number
    const int num_proc = GMEnv_num_proc_field(env);
    const int proc_num = GMEnv_proc_num_field(env);
    dm->num_field_local = gm_ceil_i8(dm->num_field_active, num_proc);
    dm->num_field = dm->num_field_local * num_proc;
    // Lower procs fully packed with active values
    // Upper procs fully inactive
    // One proc in between may be mixed
    // NOTE: see below for possible more inactive fields in packedfield
    const size_t nfl = dm->num_field_local;
    const size_t nfa = dm->num_field_active;
    dm->num_field_active_local =
      nfa <= nfl * proc_num       ? 0 :
      nfa >= nfl * (proc_num + 1) ? nfl :
                                    nfa - nfl * proc_num;
    // Pack in lower procs with no gaps to ensure indexing of actives works
    // right independent of decomposition
    dm->field_base = nfl * proc_num <= nfa ? nfl * proc_num : nfa;

  } // if fields_by_local

  //--------------------
  // Check the sizes
  //--------------------

  GMInsist(dm->num_field_active >= 0 && dm->num_field_active <= dm->num_field);
  GMInsist(dm->num_field_active_local >= 0 &&
           dm->num_field_active_local <= dm->num_field_local);

  mpi_code = MPI_Allreduce(&dm->num_field_local, &sum, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_field(env));
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Allreduce.");
  GMInsist(sum == dm->num_field_local * GMEnv_num_proc_field(env) &&
           "Every process must have the same number of fields.");
  GMInsist(sum == dm->num_field &&
           "Error in local/global sizes computation.");

  mpi_code = MPI_Allreduce(&dm->num_field_active_local, &sum, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_field(env));
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Allreduce.");
  GMInsist(sum == dm->num_field_active &&
           "Error in local/global sizes computation.");

  //--------------------
  // Element sizes
  //--------------------

  const int bits_per_byte = 8;

  switch (vectors_data_type_id) {
    //--------------------
    case GM_DATA_TYPE_FLOAT: {
      dm->num_bits_per_field = bits_per_byte * sizeof(GMFloat);
      dm->num_bits_per_packedfield = bits_per_byte * sizeof(GMFloat);
    } break;
    //--------------------
    case GM_DATA_TYPE_BITS2: {
      dm->num_bits_per_field = GM_BITS2_MAX_VALUE_BITS;
      dm->num_bits_per_packedfield = bits_per_byte * sizeof(GMBits2x64);
      // By design can only store this number of fields for this metric
      // TODO: later may be able to permit higher via rounding -
      // have 2-way method on-proc be exact, then for 3-way combining
      // or num_proc_field>1 drop low order bits to allow to fit.
      const int table_entry_limit =
        GMEnv_metric_type(env) == GM_METRIC_TYPE_DUO ?
        1 : 1 << GMEnv_num_way(env);
      GMInsistInterface(env,
               ((uint64_t)(table_entry_limit * dm->num_field)) <
                       (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS)
                && "Number of fields requested is too large for this metric");
    } break;
    //--------------------
    default:
      GMInsist(false && "Invalid vectors_data_type_id.");
  } //---switch---

  GMInsist(dm->num_bits_per_packedfield % bits_per_byte == 0 &&
           "Error in size computation.");

  dm->num_field_per_packedfield = dm->num_bits_per_packedfield /
                                  dm->num_bits_per_field;

  //--------------------
  // Packedfield counts
  //--------------------

  dm->num_packedfield_local =
      gm_ceil_i8(dm->num_field_local * dm->num_bits_per_field,
                 dm->num_bits_per_packedfield);

  //--------------------
  // Number of non-active fields on proc.
  //--------------------

  dm->num_pad_field_local =
    dm->num_packedfield_local *
    dm->num_field_per_packedfield -
    dm->num_field_active_local;

  //--------------------
  // tc memory
  //--------------------

  gm_tc_bufs_malloc(dm->num_vector_local, dm->num_field_local,
                    dm->num_packedfield_local, dm->tc_bufs, env);
}

//-----------------------------------------------------------------------------
// (Pseudo) destructor

void GMDecompMgr_destroy(GMDecompMgr* dm, GMEnv* env) {
  GMInsist(dm && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  //--------------------
  // tc memory
  //--------------------

  gm_tc_bufs_free(dm->tc_bufs, env);
}

//-----------------------------------------------------------------------------
