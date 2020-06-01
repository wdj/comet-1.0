//-----------------------------------------------------------------------------
/*!
 * \file   linalg.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Interface to generalized linear algebra functions, e.g. MAGMA.
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

#ifdef USE_CUDA
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "magma_mgemm2.h"
#include "magma_mgemm2_lapack.h"
#include "magma_mgemm3.h"
#include "magma_mgemm3_lapack.h"
#include "magma_mgemm4.h"
#include "magma_mgemm4_lapack.h"
#include "magma_mgemm5.h"
#include "magma_mgemm5_lapack.h"
#endif

#include "env.hh"
#include "assertions.hh"
#include "linalg_accel.hh"
#include "decomp_mgr.hh"
#include "linalg.hh"

//=============================================================================
// Helpers

static bool use_minproduct(GMEnv* env) {
  return GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK;
}

static bool use_mgemm2(GMEnv* env) {
  return GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
         GMEnv_num_way(env) == GM_NUM_WAY_2 && ! env->sparse;
}

static bool use_mgemm3(GMEnv* env) {
  return GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
         GMEnv_num_way(env) == GM_NUM_WAY_3 && ! env->sparse;
}

static bool use_mgemm4(GMEnv* env) {
  return GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse;
}

static bool use_mgemm5(GMEnv* env) {
  return GMEnv_metric_type(env) == GM_METRIC_TYPE_DUO;
}

//=============================================================================
/*---Magma setup, teardown---*/

void gm_linalg_initialize(GMEnv* env) {
  GMInsist(env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  // need magma blasSetKernelStream -- see
  // http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
  // page 14

#ifdef USE_CUDA
  if (use_minproduct(env)) { //--------------------

    magma_minproduct_int_t magma_code = magma_minproduct_init();
    GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_init.");
    magma_code = magma_minproductblasSetKernelStream(GMEnv_stream_compute(env));
    GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproductblasSetKernelStream.");

  } else if (use_mgemm4(env)) { //--------------------

    magma_mgemm4_int_t magma_code = magma_mgemm4_init();
    GMInsist(magma_code == MAGMA_mgemm4_SUCCESS &&
                   "Error in call to magma_mgemm4_init.");
    magma_code = magma_mgemm4blasSetKernelStream(GMEnv_stream_compute(env));
    GMInsist(magma_code == MAGMA_mgemm4_SUCCESS &&
                   "Error in call to magma_mgemm4blasSetKernelStream.");

  } else if (use_mgemm2(env)) { //--------------------

    magma_mgemm2_int_t magma_code = magma_mgemm2_init();
    GMInsist(magma_code == MAGMA_mgemm2_SUCCESS &&
                   "Error in call to magma_mgemm2_init.");
    magma_code = magma_mgemm2blasSetKernelStream(GMEnv_stream_compute(env));
    GMInsist(magma_code == MAGMA_mgemm2_SUCCESS &&
                   "Error in call to magma_mgemm2blasSetKernelStream.");

  } else if (use_mgemm3(env)) { //--------------------

    magma_mgemm3_int_t magma_code = magma_mgemm3_init();
    GMInsist(magma_code == MAGMA_mgemm3_SUCCESS &&
                   "Error in call to magma_mgemm3_init.");
    magma_code = magma_mgemm3blasSetKernelStream(GMEnv_stream_compute(env));
    GMInsist(magma_code == MAGMA_mgemm3_SUCCESS &&
                   "Error in call to magma_mgemm3blasSetKernelStream.");

  } else if (use_mgemm5(env)) { //--------------------

    magma_mgemm5_int_t magma_code = magma_mgemm5_init();
    GMInsist(magma_code == MAGMA_mgemm5_SUCCESS &&
                   "Error in call to magma_mgemm5_init.");
    magma_code = magma_mgemm5blasSetKernelStream(GMEnv_stream_compute(env));
    GMInsist(magma_code == MAGMA_mgemm5_SUCCESS &&
                   "Error in call to magma_mgemm5blasSetKernelStream.");

  } else { //--------------------

      GMInsistInterface(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // USE_CUDA
}

//-----------------------------------------------------------------------------

void gm_linalg_finalize(GMEnv* env) {
  GMInsist(env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  // TODO: (maybe) reset kernel stream (probably not really needed)

#ifdef USE_CUDA
  if (use_minproduct(env)) { //--------------------

    magma_minproduct_int_t magma_code = magma_minproduct_finalize();
    GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_finalize.");

  } else if (use_mgemm4(env)) { //--------------------

    magma_mgemm4_int_t magma_code = magma_mgemm4_finalize();
    GMInsist(magma_code == MAGMA_mgemm4_SUCCESS &&
                   "Error in call to magma_mgemm4_finalize.");

  } else if (use_mgemm2(env)) { //--------------------

    magma_mgemm2_int_t magma_code = magma_mgemm2_finalize();
    GMInsist(magma_code == MAGMA_mgemm2_SUCCESS &&
                   "Error in call to magma_mgemm2_finalize.");

  } else if (use_mgemm3(env)) { //--------------------

    magma_mgemm3_int_t magma_code = magma_mgemm3_finalize();
    GMInsist(magma_code == MAGMA_mgemm3_SUCCESS &&
                   "Error in call to magma_mgemm3_finalize.");

  } else if (use_mgemm5(env)) { //--------------------

    magma_mgemm5_int_t magma_code = magma_mgemm5_finalize();
    GMInsist(magma_code == MAGMA_mgemm5_SUCCESS &&
                   "Error in call to magma_mgemm5_finalize.");

  } else { //--------------------

      GMInsistInterface(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // USE_CUDA
}

//=============================================================================
/*---Allocate/free host and device memory---*/

void gm_linalg_malloc(GMMirroredBuf* p, size_t dim0, size_t dim1, GMEnv* env) {
  GMInsist(p && env);
  GMInsist(dim0 + 1 >= 1 && dim1 + 1 >= 1);

  *p = GMMirroredBuf_null();

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

#ifdef USE_CUDA
  p->dim0 = dim0;
  p->dim1 = dim1;

  const size_t n = dim0 * dim1;

  if (use_minproduct(env)) { //--------------------

    magma_minproduct_int_t magma_code = 0;

    if (GM_FP_PRECISION_DOUBLE) {
      magma_code = magma_minproduct_dmalloc_pinned((double**)&p->h, n);
      GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_dmalloc_pinned,"
                   " possibly due to insufficient memory.");
    } else {
      magma_code = magma_minproduct_smalloc_pinned((float**)&p->h, n);
      GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_smalloc_pinned,"
                   " possibly due to insufficient memory.");
    }
    GMFloat_fill_nan((GMFloat*)p->h, n);

    if (GM_FP_PRECISION_DOUBLE) {
      magma_code = magma_minproduct_dmalloc((double**)&p->d, n);
      GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_dmalloc,"
                   " possibly due to insufficient memory.");
    } else {
      magma_code = magma_minproduct_smalloc((float**)&p->d, n);
      GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_smalloc,"
                   " possibly due to insufficient memory.");
    }
    // TODO: ? fill GPU memory with NaNs

    p->size = n*sizeof(GMFloat);
    gm_cpu_mem_inc(p->size, env);
    gm_gpu_mem_inc(p->size, env);

  } else if (use_mgemm4(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_int_t magma_code = 0;

    magma_code = magma_mgemm4_zmalloc_pinned((Float_t**)&p->h, n);
    GMInsist(magma_code == MAGMA_mgemm4_SUCCESS &&
                   "Error in call to magma_mgemm4_zmalloc_pinned,"
                   " possibly due to insufficient memory.");

    magma_code = magma_mgemm4_zmalloc((Float_t**)&p->d, n);
    GMInsist(magma_code == MAGMA_mgemm4_SUCCESS &&
                   "Error in call to magma_mgemm4_zmalloc,"
                   " possibly due to insufficient memory.");

    p->size = n*sizeof(Float_t);
    gm_cpu_mem_inc(p->size, env);
    gm_gpu_mem_inc(p->size, env);

  } else if (use_mgemm2(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_int_t magma_code = 0;

    magma_code = magma_mgemm2_zmalloc_pinned((Float_t**)&p->h, n);
    GMInsist(magma_code == MAGMA_mgemm2_SUCCESS &&
                   "Error in call to magma_mgemm2_zmalloc_pinned,"
                   " possibly due to insufficient memory.");

    magma_code = magma_mgemm2_zmalloc((Float_t**)&p->d, n);
    GMInsist(magma_code == MAGMA_mgemm2_SUCCESS &&
                   "Error in call to magma_mgemm2_zmalloc,"
                   " possibly due to insufficient memory.");

    p->size = n*sizeof(Float_t);
    gm_cpu_mem_inc(p->size, env);
    gm_gpu_mem_inc(p->size, env);

  } else if (use_mgemm3(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_int_t magma_code = 0;

    magma_code = magma_mgemm3_zmalloc_pinned((Float_t**)&p->h, n);
    GMInsist(magma_code == MAGMA_mgemm3_SUCCESS &&
                   "Error in call to magma_mgemm3_zmalloc_pinned,"
                   " possibly due to insufficient memory.");

    magma_code = magma_mgemm3_zmalloc((Float_t**)&p->d, n);
    GMInsist(magma_code == MAGMA_mgemm3_SUCCESS &&
                   "Error in call to magma_mgemm3_zmalloc,"
                   " possibly due to insufficient memory.");

    p->size = n*sizeof(Float_t);
    gm_cpu_mem_inc(p->size, env);
    gm_gpu_mem_inc(p->size, env);

  } else if (use_mgemm5(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_int_t magma_code = 0;

    magma_code = magma_mgemm5_zmalloc_pinned((Float_t**)&p->h, n);
    GMInsist(magma_code == MAGMA_mgemm5_SUCCESS &&
                   "Error in call to magma_mgemm5_zmalloc_pinned,"
                   " possibly due to insufficient memory.");

    magma_code = magma_mgemm5_zmalloc((Float_t**)&p->d, n);
    GMInsist(magma_code == MAGMA_mgemm5_SUCCESS &&
                   "Error in call to magma_mgemm5_zmalloc,"
                   " possibly due to insufficient memory.");

    p->size = n*sizeof(Float_t);
    gm_cpu_mem_inc(p->size, env);
    gm_gpu_mem_inc(p->size, env);

  } else { //--------------------

      GMInsistInterface(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // USE_CUDA

  GMInsist(p->h && "Invalid host pointer created in gm_linalg_malloc,"
                   " possibly due to insufficient memory.");
  GMInsist(p->d && "Invalid device pointer created in gm_linalg_malloc,"
                   " possibly due to insufficient memory.");
  p->is_alias = false;
}

//-----------------------------------------------------------------------------

void gm_linalg_free(GMMirroredBuf* p, GMEnv* env) {
  GMInsist(p && env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  GMInsist(! p->is_alias);

#ifdef USE_CUDA
  const size_t size = p->size;

  if (use_minproduct(env)) { //--------------------

    magma_minproduct_int_t magma_code = magma_minproduct_free_pinned(p->h);
    GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
             "Failure in call to magma_minproduct_free_pinned.");
    magma_code = magma_minproduct_free(p->d);
    GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
             "Failure in call to magma_minproduct_free.");

    gm_cpu_mem_dec(size, env);
    gm_gpu_mem_dec(size, env);

  } else if (use_mgemm4(env)) { //--------------------

    magma_mgemm4_int_t magma_code = magma_mgemm4_free_pinned(p->h);
    GMInsist(magma_code == MAGMA_mgemm4_SUCCESS &&
             "Failure in call to magma_mgemm4_free_pinned.");
    magma_code = magma_mgemm4_free(p->d);
    GMInsist(magma_code == MAGMA_mgemm4_SUCCESS &&
             "Failure in call to magma_mgemm4_free.");

    gm_cpu_mem_dec(size, env);
    gm_gpu_mem_dec(size, env);

  } else if (use_mgemm2(env)) { //--------------------

    magma_mgemm2_int_t magma_code = magma_mgemm2_free_pinned(p->h);
    GMInsist(magma_code == MAGMA_mgemm2_SUCCESS &&
             "Failure in call to magma_mgemm2_free_pinned.");
    magma_code = magma_mgemm2_free(p->d);
    GMInsist(magma_code == MAGMA_mgemm2_SUCCESS &&
             "Failure in call to magma_mgemm2_free.");

    gm_cpu_mem_dec(size, env);
    gm_gpu_mem_dec(size, env);

  } else if (use_mgemm3(env)) { //--------------------

    magma_mgemm3_int_t magma_code = magma_mgemm3_free_pinned(p->h);
    GMInsist(magma_code == MAGMA_mgemm3_SUCCESS &&
             "Failure in call to magma_mgemm3_free_pinned.");
    magma_code = magma_mgemm3_free(p->d);
    GMInsist(magma_code == MAGMA_mgemm3_SUCCESS &&
             "Failure in call to magma_mgemm3_free.");

    gm_cpu_mem_dec(size, env);
    gm_gpu_mem_dec(size, env);

  } else if (use_mgemm5(env)) { //--------------------

    magma_mgemm5_int_t magma_code = magma_mgemm5_free_pinned(p->h);
    GMInsist(magma_code == MAGMA_mgemm5_SUCCESS &&
             "Failure in call to magma_mgemm5_free_pinned.");
    magma_code = magma_mgemm5_free(p->d);
    GMInsist(magma_code == MAGMA_mgemm5_SUCCESS &&
             "Failure in call to magma_mgemm5_free.");

    gm_cpu_mem_dec(size, env);
    gm_gpu_mem_dec(size, env);

  } else { //--------------------

      GMInsistInterface(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // USE_CUDA
}

//-----------------------------------------------------------------------------

void gm_linalg_set_matrix_zero_start(GMMirroredBuf* matrix_buf,
                                     GMEnv* env) {
  GMInsist(matrix_buf && env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

#ifdef USE_CUDA
  const size_t mat_dim1 = matrix_buf->dim0;
  const size_t mat_dim2 = matrix_buf->dim1;

  // ISSUE: these MAGMA routines don't return an error code.

  if (use_minproduct(env)) { //--------------------

    if (GM_FP_PRECISION_DOUBLE) {
      magma_minproductblas_dlaset
        (Magma_minproductFull, mat_dim1, mat_dim2, (double)0, (double)0,
         (double*)matrix_buf->d, mat_dim1);
    } else {
      magma_minproductblas_slaset
        (Magma_minproductFull, mat_dim1, mat_dim2, (float)0, (float)0,
         (float*)matrix_buf->d, mat_dim1);
    }

  } else if (use_mgemm4(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_mgemm4blas_zlaset(Magma_mgemm4Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  } else if (use_mgemm2(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_mgemm2blas_zlaset(Magma_mgemm2Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  } else if (use_mgemm3(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_mgemm3blas_zlaset(Magma_mgemm3Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  } else if (use_mgemm5(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_mgemm5blas_zlaset(Magma_mgemm5Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  } else { //--------------------

      GMInsistInterface(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // USE_CUDA
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_block_start(magma_minproduct_int_t m,
                                magma_minproduct_int_t n,
                                magma_minproduct_int_t k,
                                void* dA,
                                magma_minproduct_int_t ldda,
                                void* dB,
                                magma_minproduct_int_t lddb,
                                void* dC,
                                magma_minproduct_int_t lddc,
                                bool is_beta_one,
                                GMEnv* env) {
  GMInsist(dA && dB && dC && env);
  GMInsist(m >= 0 && n >= 0 && k >= 0);
  GMInsist(ldda >= 0 && lddb >= 0);
  GMInsist(GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU);

#ifdef USE_CUDA
  {
    int TransA = 1;
    int TransB = 0;

    magma_minproduct_int_t Am = ( ! TransA ? m : k);
    magma_minproduct_int_t An = ( ! TransA ? k : m);
    magma_minproduct_int_t Bm = ( ! TransB ? k : n);
    magma_minproduct_int_t Bn = ( ! TransB ? n : k);
    size_t sizeA = (size_t) ldda * (An - 1) + Am;
    size_t sizeB = (size_t) lddb * (Bn - 1) + Bm;

    size_t CUBLAS_MAX_1DBUF_SIZE = ((1 << 27) - 512);
    GMInsist((! (sizeA >= CUBLAS_MAX_1DBUF_SIZE ||
                 sizeB >= CUBLAS_MAX_1DBUF_SIZE )) &&
             "Error in MAGMA block sizes.");
  }

  // ISSUE: these MAGMA routines don't return an error code.

  if (use_minproduct(env)) { //--------------------

    const GMFloat alpha = 1;
    const GMFloat beta = is_beta_one ? 1 : 0;

    if (GM_FP_PRECISION_DOUBLE) {
      magma_minproductblas_dgemm
        (Magma_minproductTrans, Magma_minproductNoTrans, m, n, k, alpha,
         (double*)dA, ldda, (double*)dB, lddb, beta, (double*)dC, lddc);
      GMInsist(GMEnv_accel_last_call_succeeded(env) &&
               "Failure in call to magma_minproductblas_dgemm.");
    } else {
      magma_minproductblas_sgemm
        (Magma_minproductTrans, Magma_minproductNoTrans, m, n, k, alpha,
         (float*)dA, ldda, (float*)dB, lddb, beta, (float*)dC, lddc);
      GMInsist(GMEnv_accel_last_call_succeeded(env) &&
               "Failure in call to magma_minproductblas_sgemm.");
    }

    env->ops_local += 2 * m * (double)n * (double)k;

  } else if (use_mgemm4(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    magma_mgemm4blas_zgemm(Magma_mgemm4Trans, Magma_mgemm4NoTrans, m, n, k,
                           alpha, (Float_t*)dA, ldda, (Float_t*)dB, lddb,
                           beta, (Float_t*)dC, lddc);
    GMInsist(GMEnv_accel_last_call_succeeded(env) &&
             "Failure in call to magma_mgemm4blas_zgemm.");

  } else if (use_mgemm2(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    magma_mgemm2blas_zgemm(Magma_mgemm2Trans, Magma_mgemm2NoTrans, m, n, k,
                           alpha, (Float_t*)dA, ldda, (Float_t*)dB, lddb,
                           beta, (Float_t*)dC, lddc);
    GMInsist(GMEnv_accel_last_call_succeeded(env) &&
             "Failure in call to magma_mgemm2blas_zgemm.");

  } else if (use_mgemm3(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    magma_mgemm3blas_zgemm(Magma_mgemm3Trans, Magma_mgemm3NoTrans, m, n, k,
                           alpha, (Float_t*)dA, ldda, (Float_t*)dB, lddb,
                           beta, (Float_t*)dC, lddc);
    GMInsist(GMEnv_accel_last_call_succeeded(env) &&
             "Failure in call to magma_mgemm3blas_zgemm.");

  } else if (use_mgemm5(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

//printf("%i %i %i %i %i %i\n", (int)m, (int)n, (int)k, (int)ldda, (int)lddb, (int)lddc);
    magma_mgemm5blas_zgemm(Magma_mgemm5Trans, Magma_mgemm5NoTrans, m, n, k,
                           alpha, (Float_t*)dA, ldda, (Float_t*)dB, lddb,
                           beta, (Float_t*)dC, lddc);
    GMInsist(GMEnv_accel_last_call_succeeded(env) &&
             "Failure in call to magma_mgemm5blas_zgemm.");

  } else { //--------------------

      GMInsistInterface(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // USE_CUDA
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_start(magma_minproduct_int_t m,
                          magma_minproduct_int_t n,
                          magma_minproduct_int_t k,
                          void* dA,
                          magma_minproduct_int_t ldda,
                          void* dB,
                          magma_minproduct_int_t lddb,
                          void* dC,
                          magma_minproduct_int_t lddc,
                          GMDecompMgr* dm,
                          GMEnv* env) {
  GMInsist(dA && dB && dC && env);
  GMInsist(m >= 0 && n >= 0 && k >= 0);
  GMInsist(ldda >= 0 && lddb >= 0);
  GMInsist(GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU);

#ifdef USE_CUDA
  if (m==0 || n==0 || k==0) {
    return;
  }

  if ((GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC ||
       GMEnv_metric_type(env) == GM_METRIC_TYPE_DUO) && env->tc) {
    gm_tc_gemm_start(m, n, k, dA, ldda, dB, lddb, dC, lddc, dm->tc_bufs, env);
    return;
  }

  const size_t rows = k;
  const size_t cols_A = m;
  const size_t cols_B = n;

  const size_t elt_size =
    GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK ? sizeof(GMFloat) :
   (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse) ?
                                         sizeof(magma_mgemm4DoubleComplex) :
   (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
    GMEnv_num_way(env) == GM_NUM_WAY_2) ? sizeof(magma_mgemm2DoubleComplex) :
   (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
    GMEnv_num_way(env) == GM_NUM_WAY_3) ? sizeof(magma_mgemm3DoubleComplex) :
   (GMEnv_metric_type(env) == GM_METRIC_TYPE_DUO) ?
                                         sizeof(magma_mgemm5DoubleComplex) : 0;
  GMInsist(elt_size > 0 && "Error in gemm block calculation.");

  const size_t align_factor = 128 / elt_size;
  const size_t max_elts = (1 << 27) - 512;

  /*---TODO: can we improve aspect ratios of submatrices---*/
//  const size_t max_rows_per_block_raw = (1 << 14);
//  const size_t max_cols_per_block_raw = max_elts / max_rows_per_block_raw;

  const size_t max_rows_per_block_raw = rows + align_factor;
  const size_t max_cols_per_block_raw = max_elts / rows;

  const size_t max_rows_per_block = (max_rows_per_block_raw / align_factor)
                                                            * align_factor;
  const size_t max_cols_per_block = (max_cols_per_block_raw / align_factor)
                                                            * align_factor;

  GMInsist(max_rows_per_block != 0 && "Error in gemm block calculation.");
  GMInsist(max_cols_per_block != 0 && "Error in gemm block calculation.");

  const size_t cols_per_block_A = gm_min_i8(cols_A, max_cols_per_block);
  const size_t cols_per_block_B = gm_min_i8(cols_B, max_cols_per_block);

  const size_t rows_per_block = gm_min_i8(rows, max_rows_per_block);

  for (size_t row_base=0; row_base<rows; row_base+=rows_per_block) {
    const size_t rows_remaining = rows - row_base;
    const size_t rows_this = gm_min_i8(rows_remaining, rows_per_block);

    for (size_t col_A_base=0; col_A_base<cols_A; col_A_base+=cols_per_block_A) {
      const size_t cols_A_remaining = cols_A - col_A_base;
      const size_t cols_A_this = gm_min_i8(cols_A_remaining, cols_per_block_A);

      void* dA_this = (char*)dA + (row_base + ldda*col_A_base)*elt_size;

      for (size_t col_B_base=0; col_B_base<cols_B;
           col_B_base+=cols_per_block_B) {

        const size_t cols_B_remaining = cols_B - col_B_base;
        const size_t cols_B_this = gm_min_i8(cols_B_remaining,
                                             cols_per_block_B);

        void* dB_this = (char*)dB + (row_base + ldda*col_B_base)*elt_size;

        void* dC_this = (char*)dC + (col_A_base + lddc*col_B_base)*elt_size;

        gm_linalg_gemm_block_start(cols_A_this, cols_B_this, rows_this,
          dA_this, ldda, dB_this, lddb, dC_this, lddc, row_base > 0,  env);
      }
    }
  }
#endif // USE_CUDA
}

//-----------------------------------------------------------------------------
/*---Wait for any computation on the GPU to complete---*/

void gm_compute_wait(GMEnv* env) {
  GMInsist(env);

  GMEnv_stream_synchronize(GMEnv_stream_compute(env), env);
}

//=============================================================================
/*---Start/end transfer of generic matrix to GPU---*/

void gm_linalg_set_matrix_start(GMMirroredBuf* matrix_buf, GMEnv* env) {
  GMInsist(matrix_buf && env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

#ifdef USE_CUDA
  const size_t mat_dim1 = matrix_buf->dim0;
  const size_t mat_dim2 = matrix_buf->dim1;

  /*---Send vectors to GPU---*/

  // ISSUE: these MAGMA routines don't return an error code.

  if (use_minproduct(env)) { //--------------------

    if (GM_FP_PRECISION_DOUBLE) {
      magma_minproduct_dsetmatrix_async(
        mat_dim1, mat_dim2, (double*)matrix_buf->h, mat_dim1,
        (double*)matrix_buf->d, mat_dim1, GMEnv_stream_togpu(env));
    } else {
      magma_minproduct_ssetmatrix_async(
        mat_dim1, mat_dim2, (float*)matrix_buf->h, mat_dim1,
        (float*)matrix_buf->d, mat_dim1, GMEnv_stream_togpu(env));
    }

  } else if (use_mgemm4(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_zsetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->h,
                                  mat_dim1, (Float_t*)matrix_buf->d, mat_dim1,
                                  GMEnv_stream_togpu(env));

  } else if (use_mgemm2(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_zsetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->h,
                                  mat_dim1, (Float_t*)matrix_buf->d, mat_dim1,
                                  GMEnv_stream_togpu(env));

  } else if (use_mgemm3(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_zsetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->h,
                                  mat_dim1, (Float_t*)matrix_buf->d, mat_dim1,
                                  GMEnv_stream_togpu(env));

  } else if (use_mgemm5(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_zsetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->h,
                                  mat_dim1, (Float_t*)matrix_buf->d, mat_dim1,
                                  GMEnv_stream_togpu(env));

  } else { //--------------------

      GMInsistInterface(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // USE_CUDA
}

//-----------------------------------------------------------------------------

void gm_linalg_set_matrix_wait(GMEnv* env) {
  GMInsist(env);

  GMEnv_stream_synchronize(GMEnv_stream_togpu(env), env);
}

//=============================================================================
/*---Start/end transfer of generic matrix from GPU---*/

void gm_linalg_get_matrix_start(GMMirroredBuf* matrix_buf,
                                GMEnv* env) {
  GMInsist(matrix_buf && env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

#ifdef USE_CUDA
  const size_t mat_dim1 = matrix_buf->dim0;
  const size_t mat_dim2 = matrix_buf->dim1;

  /*---Get vectors from GPU---*/

  // ISSUE: these MAGMA routines don't return an error code.

  if (use_minproduct(env)) { //--------------------

    if (GM_FP_PRECISION_DOUBLE) {
      magma_minproduct_dgetmatrix_async(
        mat_dim1, mat_dim2, (double*)matrix_buf->d, mat_dim1,
        (double*)matrix_buf->h, mat_dim1, GMEnv_stream_fromgpu(env));
    } else {
      magma_minproduct_sgetmatrix_async(
        mat_dim1, mat_dim2, (float*)matrix_buf->d, mat_dim1,
        (float*)matrix_buf->h, mat_dim1, GMEnv_stream_fromgpu(env));
    }

  } else if (use_mgemm4(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_zgetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->d,
                                  mat_dim1, (Float_t*)matrix_buf->h, mat_dim1,
                                  GMEnv_stream_fromgpu(env));

  } else if (use_mgemm2(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_zgetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->d,
                                  mat_dim1, (Float_t*)matrix_buf->h, mat_dim1,
                                  GMEnv_stream_fromgpu(env));

  } else if (use_mgemm3(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_zgetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->d,
                                  mat_dim1, (Float_t*)matrix_buf->h, mat_dim1,
                                  GMEnv_stream_fromgpu(env));

  } else if (use_mgemm5(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_zgetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->d,
                                  mat_dim1, (Float_t*)matrix_buf->h, mat_dim1,
                                  GMEnv_stream_fromgpu(env));

  } else { //--------------------

      GMInsistInterface(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // USE_CUDA
}

//-----------------------------------------------------------------------------

void gm_linalg_get_matrix_wait(GMEnv* env) {
  GMInsist(env);

  GMEnv_stream_synchronize(GMEnv_stream_fromgpu(env), env);
}

//-----------------------------------------------------------------------------
