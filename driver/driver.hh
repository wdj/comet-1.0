//-----------------------------------------------------------------------------
/*!
 * \file   driver.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions, header.
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

#ifndef _gm_driver_hh_
#define _gm_driver_hh_

#include "env.hh"
#include "checksum.hh"

//=============================================================================
/*---Struct to hold driver options (options not in GMEnv)---*/

typedef struct {
  int num_field_local;
  int num_vector_local;
  size_t num_field;
  size_t num_vector;
  size_t num_field_active;
  size_t num_vector_active;
  bool num_field_local_initialized;
  bool num_field_active_initialized;
  bool num_vector_local_initialized;
  bool num_vector_active_initialized;
  int verbosity;
  int stage_min_0based;
  int stage_max_0based;
  int phase_min_0based;
  int phase_max_0based;
  char* input_file_path;
  char* metrics_file_path_stub;
  int problem_type;
  size_t num_incorrect;
  double max_incorrect_diff;
  double threshold;
  bool checksum;
} DriverOptions;

enum {
  GM_PROBLEM_TYPE_RANDOM = 1,
  GM_PROBLEM_TYPE_ANALYTIC = 2
};

//=============================================================================

//void finish_parsing(int argc, char** argv, DriverOptions* do_, GMEnv* env);

void perform_run(int argc, char** argv, const char* const description,
                 MPI_Comm base_comm = MPI_COMM_WORLD, GMEnv* env = 0);

void perform_run(const char* const options,
                 MPI_Comm base_comm = MPI_COMM_WORLD, GMEnv* env = 0);

void perform_run(CoMet::Checksum& cksum, int argc, char** argv,
                 const char* const description,
                 MPI_Comm base_comm = MPI_COMM_WORLD, GMEnv* env = 0);

void perform_run(CoMet::Checksum& cksum, const char* const options,
                 MPI_Comm base_comm = MPI_COMM_WORLD, GMEnv* env = 0);

//=============================================================================

#endif // _gm_driver_hh_

//-----------------------------------------------------------------------------
