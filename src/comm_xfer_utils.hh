//-----------------------------------------------------------------------------
/*!
 * \file   comm_xfer_utils.hh
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Utilities for communication and CPU/GPU transfer of vectors, metrics.
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

#ifndef _gm_comm_xfer_utils_hh_
#define _gm_comm_xfer_utils_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"

//=============================================================================

MPI_Request gm_send_vectors_start(GMVectors* vectors, int proc_num, int mpi_tag,
                                  GMEnv* env);

MPI_Request gm_recv_vectors_start(GMVectors* vectors, int proc_num, int mpi_tag,
                                  GMEnv* env);

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

//--------------------

void gm_reduce_metrics(GMMetrics* metrics,
                       GMMirroredBuf* metrics_buf_target,
                       GMMirroredBuf* metrics_buf_source,
                       GMEnv* env);

MPI_Request gm_reduce_metrics_start(GMMetrics* metrics,
                                    GMMirroredBuf* metrics_buf_target,
                                    GMMirroredBuf* metrics_buf_source,
                                    GMEnv* env);

void gm_reduce_metrics_wait(MPI_Request* mpi_request, GMEnv* env);

//--------------------

void gm_set_vectors_start(GMVectors* vectors,
                          GMMirroredBuf* vectors_buf,
                          GMEnv* env);

void gm_set_vectors_wait(GMEnv* env);

void gm_get_metrics_start(GMMetrics* metrics,
                          GMMirroredBuf* metrics_buf,
                          GMEnv* env);

void gm_get_metrics_wait(GMMetrics* metrics,
                         GMMirroredBuf* metrics_buf,
                         GMEnv* env);

//=============================================================================

#endif // _gm_comm_xfer_utils_hh_

//-----------------------------------------------------------------------------
