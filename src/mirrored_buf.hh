//-----------------------------------------------------------------------------
/*!
 * \file   mirrored_buf.hh
 * \author Wayne Joubert
 * \date   Thu Aug  3 15:04:05 EDT 2017
 * \brief  Struct/code to manage dual CPU/GPU reflected arrays.
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

#ifndef _gm_mirrored_buf_hh_
#define _gm_mirrored_buf_hh_

#include "env.hh"

//=============================================================================

typedef struct {
  void* __restrict__ h;
  void* __restrict__ d;
  size_t size;
  size_t dim0;
  size_t dim1;
  bool is_alias;
} GMMirroredBuf;

// TODO: is it appropriate to put copy to host / copy to device fns here

//-----------------------------------------------------------------------------

GMMirroredBuf GMMirroredBuf_null(void);

void GMMirroredBuf_create(GMMirroredBuf* p, size_t dim0, size_t dim1,
                          GMEnv* env);

void GMMirroredBuf_create(GMMirroredBuf* p, GMMirroredBuf* p_old, size_t dim0,
                          GMEnv* env);

void GMMirroredBuf_destroy(GMMirroredBuf* p, GMEnv* env);

//-----------------------------------------------------------------------------
/// \brief Mirrored buf element accessor.

template<typename T>
static T& GMMirroredBuf_elt(GMMirroredBuf* p, int i0, int i1) {
  GMAssert(p);
  GMAssert(i0 >= 0 && (size_t)i0 < p->dim0);
  GMAssert(i1 >= 0 && (size_t)i1 < p->dim1);

  return ((T*)(p->h))[i0 + p->dim0 * i1];
}

//-----------------------------------------------------------------------------
/// \brief Mirrored buf const element accessor.

template<typename T>
static T GMMirroredBuf_elt_const(const GMMirroredBuf* p, int i0, int i1) {
  GMAssert(p);
  GMAssert(i0 >= 0 && (size_t)i0 < p->dim0);
  GMAssert(i1 >= 0 && (size_t)i1 < p->dim1);

  return ((T*)(p->h))[i0 + p->dim0 * i1];
}

//=============================================================================

#endif // _gm_mirrored_buf_hh_

//-----------------------------------------------------------------------------
