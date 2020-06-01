//-----------------------------------------------------------------------------
/*!
 * \file   vector_sums.hh
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Compute the denominators needed by the methods.
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

#ifndef _gm_vector_sums_hh_
#define _gm_vector_sums_hh_

#include "env.hh"
#include "vectors.hh"

//=============================================================================
/*---Struct declaration---*/

typedef struct {
  GMFloat* __restrict__ sums;
  GMFloat* __restrict__ counts;
  GMFloat* __restrict__ sums_tmp_;
  GMFloat* __restrict__ counts_tmp_;
  size_t size_;
  size_t num_field_;
} GMVectorSums;

//=============================================================================
/*---Null object---*/

GMVectorSums GMVectorSums_null(void);

//=============================================================================
/*---Pseudo-constructor---*/

void GMVectorSums_create(GMVectorSums* this_,
                         int num_vector_local,
                         GMEnv* env);

//=============================================================================
/*---Pseudo-destructor---*/

void GMVectorSums_destroy(GMVectorSums* this_, GMEnv* env);

//=============================================================================
/*---Compute---*/

void GMVectorSums_compute(GMVectorSums* this_, GMVectors* vectors, GMEnv* env);

//=============================================================================
/*---Accessors---*/

GMFloat GMVectorSums_sum(const GMVectorSums* this_, int i,  GMEnv* env);

GMFloat GMVectorSums_count(const GMVectorSums* this_, int i,  GMEnv* env);

//-----------------------------------------------------------------------------

#endif // _gm_vector_sums_hh_

//-----------------------------------------------------------------------------
