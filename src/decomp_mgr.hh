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

#ifndef _gm_decomp_mgr_hh_
#define _gm_decomp_mgr_hh_

#include "env.hh"
#include "linalg_accel.hh"

//=============================================================================
// DecompMgr struct

typedef struct {
  // Field counts
  size_t num_field;
  size_t num_field_local;
  size_t num_field_active;
  size_t num_field_active_local;
  size_t field_base;
  // Vector counts
  size_t num_vector;
  size_t num_vector_local;
  size_t num_vector_active;
  size_t num_vector_active_local;
  // Packed field info
  int num_bits_per_field;
  int num_bits_per_packedfield;
  int num_field_per_packedfield;
  int num_pad_field_local;
  size_t num_packedfield_local;
  // Other
  TCBufs tc_bufs;
} GMDecompMgr;

//-----------------------------------------------------------------------------

size_t gm_num_vector_local_required(size_t num_vector_active,
                                    GMEnv* const env);

//-----------------------------------------------------------------------------
// Set to null

GMDecompMgr GMDecompMgr_null();

//-----------------------------------------------------------------------------
// (Pseudo) constructor

void GMDecompMgr_create(GMDecompMgr* dm,
                        bool fields_by_local,
                        bool vectors_by_local,
                        size_t num_field_specifier,
                        size_t num_vector_specifier,
                        int vectors_data_type_id,
                        GMEnv* env);

//-----------------------------------------------------------------------------
// (Pseudo) destructor

void GMDecompMgr_destroy(GMDecompMgr* dm, GMEnv* env);

//=============================================================================

#endif // _gm_decomp_mgr_hh_

//-----------------------------------------------------------------------------
