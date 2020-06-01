//-----------------------------------------------------------------------------
/*!
 * \file   mirrored_buf.cc
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

#include "linalg.hh"
#include "mirrored_buf.hh"

//=============================================================================

GMMirroredBuf GMMirroredBuf_null(void) {
  GMMirroredBuf p;
  p.h = NULL;
  p.d = NULL;
  p.size = 0;
  p.dim0 = 0;
  p.dim1 = 0;
  p.is_alias = false;
  return p;
}

//-----------------------------------------------------------------------------
/// \brief Mirrored buf pseudo-constructor: allocate CPU and GPU arrays.

void GMMirroredBuf_create(GMMirroredBuf* p, size_t dim0, size_t dim1, 
                          GMEnv* env) {
  GMInsist(p && env);
  GMInsist(dim0 + 1 >= 1 && dim1 + 1 >= 1);

  gm_linalg_malloc(p, dim0, dim1, env);
}

//-----------------------------------------------------------------------------
/// \brief Mirrored buf pseudo-constructor: create alias to existing mirror.

void GMMirroredBuf_create(GMMirroredBuf* p, GMMirroredBuf* p_old, size_t dim0,
                          GMEnv* env) {

  GMInsist(p && p_old && env);
  GMInsist(dim0 <= p_old->dim0);

  p->h = p_old->h;
  p->d = p_old->d;
  p->size = p_old->size;
  p->dim0 = dim0;
  p->dim1 = p_old->dim1;
  p->is_alias = true;
}

//-----------------------------------------------------------------------------
/// \brief Mirrored buf pseudo-destructor

void GMMirroredBuf_destroy(GMMirroredBuf* p, GMEnv* env) {
  GMInsist(p && env);

  if (! p->is_alias) {
    gm_linalg_free(p, env);
  }

  *p = GMMirroredBuf_null();
}

//=============================================================================

//-----------------------------------------------------------------------------
