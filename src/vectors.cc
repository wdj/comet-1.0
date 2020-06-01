//-----------------------------------------------------------------------------
/*!
 * \file   vectors.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Class to manage the set of vectors taken as input to the methods.
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

#include "mpi.h"

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"

//=============================================================================
/*---Null object---*/

GMVectors GMVectors_null() {
  GMVectors result;
  memset((void*)&result, 0, sizeof(GMVectors));
  return result;
}

//=============================================================================
/*---Set unused (pad) vector entries to zero---*/

void GMVectors_initialize_pad(GMVectors* vectors, GMEnv* env) {
  GMInsist(vectors && env);

  /*---Ensure final pad words/bits of each vector are set to zero so that
       word-wise summations of bits aren't corrupted with bad trailing data---*/

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  const size_t nfal = vectors->dm->num_field_active_local;

  const size_t pfl_min = nfal / vectors->dm->num_field_per_packedfield;
  const size_t pfl_max = vectors->dm->num_packedfield_local;

  switch (vectors->data_type_id) {
    case GM_DATA_TYPE_FLOAT: {
      GMFloat zero = 0;
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          GMVectors_float_set(vectors, pfl, vl, zero, env);
        }
      }
    } break;
    case GM_DATA_TYPE_BITS2: {
      const GMBits2x64 zero = GMBits2x64_null();
      const uint64_t allbits = 0xffffffffffffffff;
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          const size_t fl = 64 * pfl;

          if (fl >= nfal) {

            GMVectors_bits2x64_set(vectors, pfl, vl, zero, env);

          } else if (fl + 32 >= nfal) {

            GMBits2x64 val = GMVectors_bits2x64_get(vectors, pfl, vl, env);
            const int shift_dist = 64 - 2*(nfal-fl);
            GMAssert(shift_dist >= 0 && shift_dist < 64);
            val.data[0] &= allbits >> shift_dist;
            val.data[1] = 0;
            GMVectors_bits2x64_set(vectors, pfl, vl, val, env);

          } else if (fl + 64 >= nfal) {

            GMBits2x64 val = GMVectors_bits2x64_get(vectors, pfl, vl, env);
            const int shift_dist = 64 - 2*(nfal-fl-32);
            GMAssert(shift_dist >= 0 && shift_dist < 64);
            val.data[1] &= allbits >> shift_dist;
            GMVectors_bits2x64_set(vectors, pfl, vl, val, env);

          } // if
        } // for pfl
      }
    } break;
    default:
      GMInsist(false && "Invalid vectors data_type_id.");
  } /*---switch---*/
}

//=============================================================================
/*---Vectors pseudo-constructor---*/

void GMVectors_create_imp_(GMVectors* vectors,
                           int data_type_id,
                           GMDecompMgr* dm,
                           GMEnv* env) {
  GMInsist(vectors && dm && env);

  vectors->data_type_id = data_type_id;
  vectors->dm = dm;
  vectors->num_field = dm->num_field;
  vectors->num_field_active = dm->num_field_active;
  vectors->num_field_local = dm->num_field_local;
  vectors->num_vector = dm->num_vector;
  vectors->num_vector_local = dm->num_vector_local;

  // Set element sizes

  vectors->num_bits_per_val = dm->num_bits_per_field;
  vectors->num_bits_per_packedval = dm->num_bits_per_packedfield;
  vectors->num_val_per_packedval = dm->num_field_per_packedfield;

  const int bits_per_byte = 8;

  // Allocation size for vector storage

  vectors->num_packedval_field_local = dm->num_packedfield_local;

  vectors->num_packedval_local =
      vectors->num_packedval_field_local * dm->num_vector_local;

  vectors->data_size = vectors->num_packedval_local *
                       (vectors->num_bits_per_packedval / bits_per_byte);

  // Set up vector storage, mirrored buffer

  vectors->buf = GMMirroredBuf_null();

  if (vectors->has_buf) {
    GMMirroredBuf_create(&(vectors->buf),vectors->num_packedval_field_local,
                         vectors->num_vector_local, env);
    vectors->data = vectors->buf.h; // alias vector storage to buf
  } else {
    vectors->data = gm_malloc(vectors->data_size, env);
  }

  // Set pad entries to zero

  GMVectors_initialize_pad(vectors, env);

//  gm_tc_bufs_malloc(env, vectors->num_vector_local,
//                    vectors->num_packedval_field_local);
}

//-----------------------------------------------------------------------------

void GMVectors_create(GMVectors* vectors,
                      int data_type_id,
                      GMDecompMgr* dm,
                      GMEnv* env) {
  GMInsist(vectors && dm && env);

  *vectors = GMVectors_null();

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  vectors->has_buf = false;

  GMVectors_create_imp_(vectors, data_type_id, dm, env);
}

//-----------------------------------------------------------------------------

void GMVectors_create_with_buf(GMVectors* vectors,
                               int data_type_id,
                               GMDecompMgr* dm,
                               GMEnv* env) {
  GMInsist(vectors && dm && env);

  *vectors = GMVectors_null();

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  vectors->has_buf = GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU;

  GMVectors_create_imp_(vectors, data_type_id, dm, env);
}

//=============================================================================
/*---Vectors pseudo-destructor---*/

void GMVectors_destroy(GMVectors* vectors, GMEnv* env) {
  GMInsist(vectors && env);
  GMInsist(vectors->data || ! GMEnv_is_proc_active(env));

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  if (vectors->has_buf) {
    GMMirroredBuf_destroy(&vectors->buf, env);
  } else {
    gm_free(vectors->data, vectors->data_size, env);
  }

  *vectors = GMVectors_null();
}


//=============================================================================
// Copy vectors to mirrored buffer

void gm_vectors_to_buf(GMMirroredBuf* vectors_buf,
                       GMVectors* vectors,
                       GMEnv* env) {
  GMInsist(vectors && vectors_buf && env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*---Copy vectors into GPU buffers if needed---*/

  switch (GMEnv_metric_type(env)) {
    case GM_METRIC_TYPE_CZEK: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_field_local; ++fl) {
          GMMirroredBuf_elt<GMFloat>(vectors_buf, fl, i) =
            GMVectors_float_get(vectors, fl, i, env);
        }
      }
    } break;
    case GM_METRIC_TYPE_CCC: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_packedval_field_local; ++fl) {
          GMMirroredBuf_elt<GMBits2x64>(vectors_buf, fl, i) =
            GMVectors_bits2x64_get(vectors, fl, i, env);
        }
      }
    } break;
    case GM_METRIC_TYPE_DUO: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_packedval_field_local; ++fl) {
          GMMirroredBuf_elt<GMBits2x64>(vectors_buf, fl, i) =
            GMVectors_bits2x64_get(vectors, fl, i, env);
        }
      }
    } break;
    default:
      GMInsistInterface(env, false && "Unimplemented metric_type.");
  } /*---case---*/
}

//=============================================================================
// Print entries of vectors

void GMVectors_print(GMVectors* vectors, GMEnv* env) {
  GMInsist(vectors && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  const int nval = vectors->dm->num_vector_active_local;
  const int nfal = vectors->dm->num_field_active_local;
  
  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
      for (int vl = 0; vl < nval; ++vl) {
        for (int fl = 0; fl < nfal; ++fl) {
          const GMFloat float_value = GMVectors_float_get(vectors, fl, vl, env);
            printf("vec_proc %i vec %i field_proc %i field %i value %e\n",
                   GMEnv_proc_num_vector_i(env), vl,
                   GMEnv_proc_num_field(env), fl, float_value);
        } /*---fl---*/
      }   /*---vl---*/
    } break;       
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/
      for (int vl = 0; vl < nval; ++vl) {
        for (int fl = 0; fl < nfal; ++fl) {
          const GMBits2 value = GMVectors_bits2_get(vectors, fl, vl, env);
            printf("vec_proc %i vec %i "
                   "field_proc %i field %i value %.1i%.1i\n",
                   GMEnv_proc_num_vector_i(env), vl,
                   GMEnv_proc_num_field(env), fl, value / 2, value % 2);
        } /*---fl---*/
      }   /*---vl---*/
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMInsist(false && "Invalid data_type_vectors.");
  } /*---switch---*/
}

//=============================================================================
// hecksum of entries.

size_t GMVectors_cksum(GMVectors* vectors, GMEnv* env) {
  GMInsist(vectors && env);

  if (! GMEnv_is_proc_active(env)) {
    return 0;
  }

  return gm_array_cksum((unsigned char*)(vectors->data), vectors->data_size);
}

//-----------------------------------------------------------------------------
