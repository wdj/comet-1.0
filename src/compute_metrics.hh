//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Top-level function to calculate metrics.
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

#ifndef _gm_compute_metrics_hh_
#define _gm_compute_metrics_hh_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "compute_metrics_2way.hh"
#include "compute_metrics_3way.hh"

//=============================================================================

typedef struct {
  GMComputeMetrics2Way compute_metrics_2way;
  GMComputeMetrics3Way compute_metrics_3way;
} GMComputeMetrics;

//=============================================================================

void GMComputeMetrics_create(
    GMComputeMetrics* this_,
    GMDecompMgr* dm,
    GMEnv* env);

void GMComputeMetrics_destroy(
    GMComputeMetrics* this_,
    GMEnv* env);

//-----------------------------------------------------------------------------

void gm_compute_metrics(GMComputeMetrics* compute_metrics, GMMetrics* metrics,
                        GMVectors* vectors, GMEnv* env);

void gm_compute_metrics(GMMetrics* metrics, GMVectors* vectors, GMEnv* env);

//=============================================================================

#endif // _gm_compute_metrics_hh_

//-----------------------------------------------------------------------------
