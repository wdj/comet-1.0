//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way.cc
 * \author Wayne Joubert
 * \date   Thu Jan  7 10:21:09 EST 2016
 * \brief  Calculate metrics, 2-way.
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

#include "string.h"

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "comm_xfer_utils.hh"
#include "compute_metrics_2way_block_nums.hh"
#include "compute_metrics_2way_block_combine.hh"
#include "compute_metrics_2way.hh"

//=============================================================================

void GMComputeMetrics2Way_create(
    GMComputeMetrics2Way* this_,
    GMDecompMgr* dm,
    GMEnv* env) {
  GMInsist(this_ && dm && env);

  if (!(GMEnv_num_way(env) == 2 && GMEnv_all2all(env))) {
    return;
  }

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  *this_ = {0};

  GMVectorSums_create(&this_->vector_sums_onproc, dm->num_vector_local, env);
  GMVectorSums_create(&this_->vector_sums_offproc, dm->num_vector_local, env);

  for (int i = 0; i < 2; ++i) {
    GMVectors_create_with_buf(&this_->vectors_01[i],
                              GMEnv_data_type_vectors(env), dm, env);
  }

  for (int i = 0; i < 2; ++i) {
    GMMirroredBuf_create(&this_->metrics_buf_01[i],
                         dm->num_vector_local, dm->num_vector_local, env);
  }

  GMMirroredBuf_create(&this_->vectors_buf, dm->num_packedfield_local,
                       dm->num_vector_local, env);

  if (env->do_reduce) {
    GMMirroredBuf_create(&this_->metrics_tmp_buf,
                         dm->num_vector_local, dm->num_vector_local, env);
  }
}

//-----------------------------------------------------------------------------

void GMComputeMetrics2Way_destroy(
    GMComputeMetrics2Way* this_,
    GMEnv* env) {
  GMInsist(this_ && env);

  if (!(GMEnv_num_way(env) == 2 && GMEnv_all2all(env))) {
    return;
  }

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  GMVectorSums_destroy(&this_->vector_sums_onproc, env);
  GMVectorSums_destroy(&this_->vector_sums_offproc, env);

  for (int i = 0; i < 2; ++i) {
    GMVectors_destroy(&this_->vectors_01[i], env);
  }

  for (int i = 0; i < 2; ++i) {
    GMMirroredBuf_destroy(&this_->metrics_buf_01[i], env);
  }

  GMMirroredBuf_destroy(&this_->vectors_buf, env);

  if (env->do_reduce) {
    GMMirroredBuf_destroy(&this_->metrics_tmp_buf, env);
  }
}

//=============================================================================

void gm_compute_metrics_2way_notall2all(
  GMComputeMetrics2Way* this_,
  GMMetrics* metrics,
  GMVectors* vectors,
  GMEnv* env) {

  GMInsist(metrics && vectors && env);
  GMInsist(! GMEnv_all2all(env));

  // Denominator

  GMVectorSums vector_sums = GMVectorSums_null();
  GMVectorSums_create(&vector_sums, vectors->num_vector_local, env);

  GMVectorSums_compute(&vector_sums, vectors, env);

  // Numerator

  gm_linalg_initialize(env);

  const int nvl = vectors->num_vector_local;
  const int npvfl = vectors->num_packedval_field_local;

  // Allocate memory for vectors and for result 

  GMMirroredBuf vectors_buf = GMMirroredBuf_null();
  GMMirroredBuf_create(&vectors_buf, npvfl, nvl, env);

  GMMirroredBuf metrics_buf = GMMirroredBuf_null();
  GMMirroredBuf_create(&metrics_buf, nvl, nvl, env);

  GMMirroredBuf metrics_tmp_buf = GMMirroredBuf_null();
  if (env->do_reduce) {
    GMMirroredBuf_create(&metrics_tmp_buf, nvl, nvl, env);
  }

  GMMirroredBuf* metrics_buf_ptr =
      env->do_reduce ?  &metrics_tmp_buf : &metrics_buf;

  // Copy in vectors

  gm_vectors_to_buf(&vectors_buf, vectors, env);

  // Send vectors to GPU

  gm_set_vectors_start(vectors, &vectors_buf, env);
  gm_set_vectors_wait(env);

  gm_compute_2way_proc_nums_start(vectors, vectors, metrics, &vectors_buf,
                                  &vectors_buf, metrics_buf_ptr,
                                  GMEnv_proc_num_vector_i(env),
                                  true, env);
  gm_compute_wait(env);

  // Copy result from GPU

  gm_get_metrics_start(metrics, metrics_buf_ptr, env);
  gm_get_metrics_wait(metrics, metrics_buf_ptr, env);
  gm_metrics_pad_adjust(metrics, metrics_buf_ptr, env);

  // Do reduction across field procs if needed

  if (env->do_reduce) {
    gm_reduce_metrics(metrics, &metrics_buf, metrics_buf_ptr, env);
  }

  // Combine

  gm_compute_2way_proc_combine(metrics, &metrics_buf,
                               &vector_sums, &vector_sums,
                               GMEnv_proc_num_vector_i(env), true, env);

  // Terminations

  GMVectorSums_destroy(&vector_sums, env);

  GMMirroredBuf_destroy(&vectors_buf, env);
  GMMirroredBuf_destroy(&metrics_buf, env);

  if (env->do_reduce) {
    GMMirroredBuf_destroy(&metrics_tmp_buf, env);
  }

  gm_linalg_finalize(env);
}

//=============================================================================

static void lock(bool& lock_val) {
  GMInsist(! lock_val);
  lock_val = true;
};

static void unlock(bool& lock_val) {
  GMInsist(lock_val);
  lock_val = false;
};

//=============================================================================

void gm_compute_metrics_2way_all2all(
  GMComputeMetrics2Way* this_,
  GMMetrics* metrics,
  GMVectors* vectors,
  GMEnv* env) {

  GMInsist(metrics && vectors && env);
  GMInsist(GMEnv_all2all(env));

  // Initializations

  const int num_block = GMEnv_num_block_vector(env);
  const int i_block = GMEnv_proc_num_vector_i(env);

  GMVectorSums vector_sums_onproc = this_->vector_sums_onproc;
  GMVectorSums vector_sums_offproc = this_->vector_sums_offproc;

  gm_linalg_initialize(env);

  // Create double buffer of vectors objects for send/recv

  GMVectors* const & vectors_01 = this_->vectors_01;

  // Allocate GPU buffers
  // To overlap transfers with compute, set up double buffers for the
  // vectors sent to the GPU and the metrics received from the GPU.

  GMMirroredBuf* const & metrics_buf_01 = this_->metrics_buf_01;
  GMMirroredBuf& vectors_buf = this_->vectors_buf;
  GMMirroredBuf& metrics_tmp_buf = this_->metrics_tmp_buf;

  // Result matrix is diagonal block and half the blocks to the right
  // (including wraparound to left side of matrix when appropriate).
  // For even number of vector blocks, block rows of lower half of matrix
  //  have one less block to make correct count.

  const int num_proc_r = GMEnv_num_proc_repl(env);
  const int proc_num_r = GMEnv_proc_num_repl(env);

  // Flatten the proc_vector and proc_repl indices into a single index.

  const int num_proc_rv = num_block * num_proc_r;
  const int proc_num_rv = proc_num_r + num_proc_r * i_block;

  MPI_Request mpi_requests[2];

  // Prepare for loop over blocks of result.

  /*---Summary of the opertions in this loop:

    send VECTORS next step start
    recv VECTORS next step start
    set VECTORS this step wait
    compute numerators start
    get METRICS prev step wait
    combine prev step wait (GPU case)
    recv VECTORS next step wait
    set VECTORS next step start
    compute numerators wait
    get METRICS this step start
    compute denominators
    combine this step (CPU case)
    send VECTORS next step wait

  ---*/

  // Add extra step at begin/end to fill/drain pipeline.

  const int extra_step = 1;

  // Lowest/highest (block) diag to be computed for this phase,
  // measured from (block) main diag.
  // For all repl procs.

  const int j_i_offset_min = gm_bdiag_computed_min(env);
  const int j_i_offset_max = gm_bdiag_computed_max(env);
  const int j_i_offset_this_row_max = gm_block_computed_this_row_max(env);

  const int num_bdiag_computed = j_i_offset_max - j_i_offset_min;

  // Num steps to take to compute blocks
  // (note: at each step, num_proc_r processors each compute a block)
  // NOTE: num_step should be consistent within same proc_r.

  const int num_step = gm_ceil_i(num_bdiag_computed, num_proc_r);

  typedef struct {
    GMVectors* vectors_right;
    GMMirroredBuf* vectors_right_buf;
    GMMirroredBuf* metrics_buf;
    bool is_compute_step;
    bool is_first_compute_step;
    bool do_compute_block;
    bool is_main_diag;
    bool is_right_aliased;
    int step_num;
    int index_01;
    int j_i_offset;
    int j_block;
  } LoopVars;

  LoopVars vars = {0};
  LoopVars vars_prev = {0};
  LoopVars vars_next = {0};

  // Use locks to verify no race condition on a buffer.
  // Lock buffer when in use for read or write, unlock when done.

  bool lock_vectors_01_buf_h[2] = {false, false};
  bool lock_vectors_01_buf_d[2] = {false, false};
  bool lock_metrics_buf_01_h[2] = {false, false};
  bool lock_metrics_buf_01_d[2] = {false, false};
  bool lock_vectors_buf_h = false;
  bool lock_vectors_buf_d = false;
  //bool lock_metrics_tmp_buf_d = false; // Not needed
  bool lock_metrics_tmp_buf_h = false;

  //========================================
  for (int step_num = 0-extra_step; step_num < num_step+extra_step; ++step_num){
  //========================================
//double t0 = GMEnv_get_time();

    // Set per-step variables

    vars_prev = vars;
    vars = vars_next;

    vars_next.step_num = step_num + 1;
    vars_next.is_compute_step = vars_next.step_num >= 0 &&
                                vars_next.step_num < num_step;
    vars_next.is_first_compute_step = vars_next.step_num == 0;
    vars_next.index_01 = gm_mod_i(vars_next.step_num, 2);
    vars_next.j_i_offset = j_i_offset_min + vars_next.step_num * num_proc_r
                           + proc_num_r;
    vars_next.is_main_diag = vars_next.j_i_offset == 0;
    vars_next.j_block = gm_mod_i(i_block + vars_next.j_i_offset, num_block);
    vars_next.do_compute_block = vars_next.is_compute_step &&
                   vars_next.j_i_offset < j_i_offset_this_row_max;

    // Pointers to left/right-side vecs.
    // Here we are computing V^T W, for V, W containing column vectors.

    vars_next.is_right_aliased = vars_next.is_main_diag;

    GMVectors* vectors_left = vectors;
    vars_next.vectors_right = vars_next.is_right_aliased ?
      vectors_left : &vectors_01[vars_next.index_01];

    GMMirroredBuf* vectors_left_buf = &vectors_buf;
    vars_next.vectors_right_buf = vars_next.is_right_aliased ?
      vectors_left_buf : &vectors_01[vars_next.index_01].buf;

    // Pointer to metrics buffer

    vars_next.metrics_buf = &metrics_buf_01[vars_next.index_01];

    // Set up lock aliases

    bool& lock_metrics_buf_ptr_h_prev
                                  = lock_metrics_buf_01_h[vars_prev.index_01];
    bool& lock_metrics_buf_ptr_d_prev
                                  = lock_metrics_buf_01_d[vars_prev.index_01];

    bool& lock_metrics_buf_ptr_h = lock_metrics_buf_01_h[vars.index_01];
    bool& lock_metrics_buf_ptr_d = lock_metrics_buf_01_d[vars.index_01];

    bool& lock_vectors_left_buf_h = lock_vectors_buf_h;
    bool& lock_vectors_left_buf_d = lock_vectors_buf_d;

    bool& lock_vectors_right_buf_h_next = vars_next.is_right_aliased ?
      lock_vectors_left_buf_h : lock_vectors_01_buf_h[vars_next.index_01];

    bool& lock_vectors_right_buf_d_next = vars_next.is_right_aliased ?
      lock_vectors_left_buf_d : lock_vectors_01_buf_d[vars_next.index_01];

    bool& lock_vectors_right_buf_h = vars.is_right_aliased ?
      lock_vectors_left_buf_h : lock_vectors_01_buf_h[vars.index_01];

    bool& lock_vectors_right_buf_d = vars.is_right_aliased ?
      lock_vectors_left_buf_d : lock_vectors_01_buf_d[vars.index_01];

    // Prepare for sends/recvs: procs for communication

    const int proc_send = gm_mod_i(proc_num_rv
        - vars_next.j_i_offset*num_proc_r, num_proc_rv);

    const int proc_recv = gm_mod_i(proc_num_rv
        + vars_next.j_i_offset*num_proc_r, num_proc_rv);

    const bool comm_with_self = vars_next.is_main_diag;

    // Initiate sends/recvs for vecs needed on next step

    if (vars_next.is_compute_step && ! comm_with_self) {
      const int mpi_tag = step_num + 1;
      // NOTE: the following order helps performance
      GMInsist((!vars_next.is_right_aliased) &&
               "Next step should always compute off-diag block.");
      lock(lock_vectors_right_buf_h_next);
      mpi_requests[1] = gm_recv_vectors_start(vars_next.vectors_right,
                                              proc_recv, mpi_tag, env);
      mpi_requests[0] = gm_send_vectors_start(vectors_left,
                                              proc_send, mpi_tag, env);
    }

    // Send right vectors to GPU end

    if (vars.is_compute_step && vars.do_compute_block &&
        ! vars.is_right_aliased) {
      gm_set_vectors_wait(env);
      unlock(lock_vectors_right_buf_h);
      unlock(lock_vectors_right_buf_d);
    }

    // First step (for any repl or phase): send (left) vecs to GPU

    if (vars_next.is_first_compute_step) {
      lock(lock_vectors_left_buf_h);
      gm_vectors_to_buf(vectors_left_buf, vectors_left, env);
      lock(lock_vectors_left_buf_d);
      gm_set_vectors_start(vectors_left, vectors_left_buf, env);
      // TODO: examine whether overlap possible.
      // May not be possible for general repl and phase (??).
      gm_set_vectors_wait(env);
      unlock(lock_vectors_left_buf_h);
      unlock(lock_vectors_left_buf_d);
    }

    // Commence numerators computation

    if (vars.is_compute_step && vars.do_compute_block) {
      lock(lock_vectors_left_buf_d);
      if (! vars.is_right_aliased) {
        lock(lock_vectors_right_buf_d);
      }
      lock(lock_metrics_buf_ptr_d);
      gm_compute_2way_proc_nums_start(
        vectors_left, vars.vectors_right, metrics,
        vectors_left_buf, vars.vectors_right_buf, vars.metrics_buf,
        vars.j_block, vars.is_main_diag, env);
    }

    // GPU case: wait for prev step get metrics to complete, then combine.
    // Note this is hidden under GPU computation

    if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
      if (vars_prev.is_compute_step && vars_prev.do_compute_block) {
        gm_get_metrics_wait(metrics, vars_prev.metrics_buf, env);
        unlock(lock_metrics_buf_ptr_d_prev);
        unlock(lock_metrics_buf_ptr_h_prev);
        lock(lock_metrics_buf_ptr_h_prev);
        gm_metrics_pad_adjust(metrics, vars_prev.metrics_buf, env);
        unlock(lock_metrics_buf_ptr_h_prev);

        GMVectorSums* vector_sums_left = &vector_sums_onproc;
        GMVectorSums* vector_sums_right =
          vars_prev.is_main_diag
          ? &vector_sums_onproc : &vector_sums_offproc;

        //TODO: remove need to allocate metrics_tmp_buf device array
        GMMirroredBuf* metrics_buf_prev_ptr =
            env->do_reduce ?  &metrics_tmp_buf : vars_prev.metrics_buf;

        lock(lock_metrics_buf_ptr_h_prev); // semantics not perfect but ok

        if (env->do_reduce) {
          lock(lock_metrics_tmp_buf_h);
          gm_reduce_metrics(metrics, metrics_buf_prev_ptr,
                            vars_prev.metrics_buf, env);
        }

        gm_compute_2way_proc_combine(
          metrics, metrics_buf_prev_ptr,
          vector_sums_left, vector_sums_right,
          vars_prev.j_block,
          vars_prev.is_main_diag, env);

        unlock(lock_metrics_buf_ptr_h_prev); // semantics not perfect but ok

        if (env->do_reduce) {
          unlock(lock_metrics_tmp_buf_h);
        }
      }
    }

    // ISSUE: it may be possible to increase performance by swapping the
    // some code blocks below and the one code block above.  It depends
    // on the relative speeds.  If these would be put in two different
    // CPU threads, then it wouldn't matter.

    // Compute sums for denominators

    //const bool compute_sums_early = GMEnv_is_ppc64();
    const bool compute_sums_early = true;

    if (compute_sums_early) { // put it here for speed on this arch
      if (vars.is_compute_step && vars.do_compute_block) {
        //TODO: possibly move this
        if (vars.is_first_compute_step) {
          GMVectorSums_compute(&vector_sums_onproc, vectors_left, env);
        }
        if (! vars.is_main_diag) {
          GMVectorSums_compute(&vector_sums_offproc, vars.vectors_right, env);
        }
      }
    }

    // Wait for recvs to complete

    if (vars_next.is_compute_step && ! comm_with_self) {
      gm_recv_vectors_wait(&(mpi_requests[1]), env);
      GMInsist((!vars_next.is_right_aliased) &&
               "Next step should always compute off-diag block.");
      unlock(lock_vectors_right_buf_h_next);
    }

    // Send right vectors for next step to GPU start

    if (vars_next.is_compute_step && vars_next.do_compute_block &&
        ! vars_next.is_right_aliased) {
      // ISSUE: make sure not necessary if vars_next.is_right_aliased
      lock(lock_vectors_right_buf_h_next);
      lock(lock_vectors_right_buf_d_next);
      gm_set_vectors_start(vars_next.vectors_right,
                           vars_next.vectors_right_buf, env);
    }

    // Wait for numerators computation to complete

    if (vars.is_compute_step && vars.do_compute_block) {
      gm_compute_wait(env);
      unlock(lock_vectors_left_buf_d);
      if (! vars.is_right_aliased) {
        unlock(lock_vectors_right_buf_d);
      }
      unlock(lock_metrics_buf_ptr_d);
    }

    // Commence copy of completed numerators back from GPU

    if (vars.is_compute_step && vars.do_compute_block) {
      lock(lock_metrics_buf_ptr_h);
      lock(lock_metrics_buf_ptr_d);
      gm_get_metrics_start(metrics, vars.metrics_buf, env);
    }

    // Compute sums for denominators

    if (! compute_sums_early) { // put it here for speed on this arch
      if (vars.is_compute_step && vars.do_compute_block) {
        //TODO: possibly move this
        if (vars.is_first_compute_step) {
          GMVectorSums_compute(&vector_sums_onproc, vectors_left, env);
        }
        if (! vars.is_main_diag) {
          GMVectorSums_compute(&vector_sums_offproc, vars.vectors_right, env);
        }
      }
    }

    // CPU case: combine numerators, denominators to obtain final result

    if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
      if (vars.is_compute_step && vars.do_compute_block) {
        GMVectorSums* vector_sums_left = &vector_sums_onproc;
        GMVectorSums* vector_sums_right =
            vars.is_main_diag
            ? &vector_sums_onproc : &vector_sums_offproc;
        gm_get_metrics_wait(metrics, vars.metrics_buf, env); // NO-OP
        unlock(lock_metrics_buf_ptr_d);
        unlock(lock_metrics_buf_ptr_h);
        lock(lock_metrics_buf_ptr_h);
        gm_compute_2way_proc_combine(
          metrics, vars.metrics_buf, vector_sums_left,
          vector_sums_right, vars.j_block,
          vars.is_main_diag, env);
        unlock(lock_metrics_buf_ptr_h);
      }
    }

    // Wait for sends to complete

    if (vars_next.is_compute_step && ! comm_with_self) {
      gm_send_vectors_wait(&(mpi_requests[0]), env);
    }

//double t1 = GMEnv_get_time();
//printf("%i %f\n", step_num, t1-t0);
  //========================================
  } // step_num
  //========================================

  // Terminations

  for (int i=0; i<2; ++i) {
    GMInsist(!lock_vectors_01_buf_h[i]);
    GMInsist(!lock_vectors_01_buf_d[i]);
    GMInsist(!lock_metrics_buf_01_h[i]);
    GMInsist(!lock_metrics_buf_01_d[i]);
  }
  GMInsist(!lock_vectors_buf_h);
  GMInsist(!lock_vectors_buf_d);
  GMInsist(!lock_metrics_tmp_buf_h);

  gm_linalg_finalize(env);
}

//-----------------------------------------------------------------------------
