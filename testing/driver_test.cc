//-----------------------------------------------------------------------------
/*!
 * \file   driver_test.cc
 * \author Wayne Joubert
 * \date   Fri Nov  6 18:18:21 EST 2015
 * \brief  Tester for driver.
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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "gtest/gtest.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "checksum.hh"
#include "compute_metrics.hh"

#include "driver.hh"
#include "input_output.hh"
#include "test_problems.hh"

enum {PROCS_MAX = TEST_PROCS_MAX};

//=============================================================================

bool compare_2runs(const char* options1, const char* options2) {
  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

  /*---Do runs---*/

  if (proc_num == 0) {
    printf("%s\n", options1);
  }

  CoMet::Checksum checksum1;
  perform_run(checksum1, options1);

  if (proc_num == 0) {
    printf("%s\n", options2);
  }

  CoMet::Checksum checksum2;
  perform_run(checksum2, options2);

  /*---Need test result only on proc 0---*/

  const bool is_passed = proc_num != 0 ? true :
                         checksum1.is_equal(checksum2);

  return is_passed;
}

//=============================================================================

bool compare_3runs(const char* options1,
                   const char* options2,
                   const char* options3) {
  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

  /*---Do runs---*/

  if (proc_num == 0) {
    printf("%s\n", options1);
  }
  CoMet::Checksum checksum1;
  perform_run(checksum1, options1);

  if (proc_num == 0) {
    printf("%s\n", options2);
  }
  CoMet::Checksum checksum2;
  perform_run(checksum2, options2);

  if (proc_num == 0) {
    printf("%s\n", options3);
  }
  CoMet::Checksum checksum3;
  perform_run(checksum3, options3);

  /*---Need test result only on proc 0---*/

  const bool is_passed = proc_num != 0 ? true :
                         checksum1.is_equal(checksum2) &&
                         checksum1.is_equal(checksum3);
  return is_passed;
}

//=============================================================================

void test_2runs(const char* options1,
                const char* options2) {
  EXPECT_EQ(true, compare_2runs(options1, options2));
}

//=============================================================================

void create_vectors_file(const char* file_path, int num_field, int num_vector,
                         int metric_type, int num_way, int problem_type,
                         int verbosity) {

  GMEnv env_value = GMEnv_null(), *env = &env_value;
  GMEnv_create(env, MPI_COMM_WORLD, 0, NULL, NULL);
  env->metric_type_ = metric_type;
  env->num_way_ = num_way;
  GMEnv_set_num_proc(env, 1, 1, 1);

  GMDecompMgr dm_value = GMDecompMgr_null(), *dm = &dm_value;
  GMDecompMgr_create(dm, false, false, num_field, num_vector,
                     GMEnv_data_type_vectors(env), env);

  GMVectors vectors_value = GMVectors_null(), *vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), dm, env);
  set_vectors_synthetic(vectors, problem_type, verbosity, env);

  write_vectors_to_file(vectors, file_path, env);

  GMVectors_destroy(vectors, env);
  GMDecompMgr_destroy(dm, env);
  GMEnv_destroy(env);
}

//=============================================================================

void DriverTest_czek_() {

  //----------
  //---2-way, all2all no
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method GPU"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU"));

  //----------
  //---2-way, all2all yes, small
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU ",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU ",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method GPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 4 "
                    "--compute_method CPU ",
                    "--num_proc_vector 2 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 4 "
                    "--compute_method CPU ",
                    "--num_proc_vector 2 --num_field 1 --num_vector_local 2 "
                    "--compute_method GPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_3runs("--num_proc_vector 1 --num_field 2 --num_vector 5 "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 2 --num_field 2 --num_vector 5 "
                    "--compute_method CPU --all2all yes",
                    "--num_proc_vector 2 --num_field 2 --num_vector 5 "
                    "--compute_method GPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_proc_field 1 "
                    "--num_field 7 --num_vector 5 "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 1 --num_proc_field 3 "
                    "--num_field 7 --num_vector 5 "
                    "--compute_method GPU --all2all yes"));

  //----------
  //---2-way, all2all yes, large
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU"
                    " --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 2 --num_field 100 --num_vector_local 24 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 2 --num_field 100 --num_vector_local 24 "
                    "--compute_method GPU --all2all yes"));

  //----------
  //---3-way, all2all no
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU "
                    "--num_way 3",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method GPU --num_way 3"));
  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3"));

  //----------
  //---3-way, all2all yes, small
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 18 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 3 --num_field 1 --num_vector_local 6 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 18 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 3 --num_field 1 --num_vector_local 6 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_3runs("--num_proc_vector 1 --num_field 2 --num_vector 16 "
                    "--compute_method REF --all2all yes --num_way 3",
                    "--num_proc_vector 3 --num_field 2 --num_vector 16 "
                    "--compute_method CPU --all2all yes --num_way 3",
                    "--num_proc_vector 3 --num_field 2 --num_vector 16 "
                    "--compute_method GPU --all2all yes --num_way 3"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_proc_field 1 "
                    "--num_field 13 --num_vector 6 "
                    "--compute_method REF --all2all yes --num_way 3",
                    "--num_proc_vector 1 --num_proc_field 5 "
                    "--num_field 13 --num_vector 6 "
                    "--compute_method GPU --all2all yes --num_way 3"));

  //----------
  //---3-way, all2all yes, large
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3 --all2all yes",
                    "--num_proc_vector 4 --num_field 100 --num_vector_local 12 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3 --all2all yes",
                    "--num_proc_vector 4 --num_field 100 --num_vector_local 12 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  //----------
  //---num_proc_field
  //----------

  EXPECT_EQ(true, compare_2runs("--num_proc_vector 1 --num_proc_field "
                                "1 --num_field 2 --num_vector_local 2 "
                                "--compute_method CPU",
                                "--num_proc_vector 1 --num_proc_field "
                                "2 --num_field 2 --num_vector_local 2 "
                                "--compute_method GPU"));

  EXPECT_EQ(true, compare_2runs("--num_proc_vector 1 --num_proc_field "
                                "1 --num_field 2 --num_vector_local 4 "
                                "--compute_method CPU",
                                "--num_proc_vector 2 --num_proc_field "
                                "2 --num_field 2 --num_vector_local 2 "
                                "--compute_method GPU --all2all yes"));

  EXPECT_EQ(true, compare_2runs("--num_proc_vector 1 --num_proc_field "
                                "1 --num_field 2 --num_vector_local 3 "
                                "--compute_method CPU --num_way 3",
                                "--num_proc_vector 1 --num_proc_field "
                                "2 --num_field 2 --num_vector_local 3 "
                                "--compute_method GPU --num_way 3"));

  EXPECT_EQ(true,
            compare_2runs("--num_proc_vector 1 --num_proc_field 1 --num_field "
                          "2 --num_vector_local 18"
                          " --compute_method CPU --num_way 3",
                          "--num_proc_vector 3 --num_proc_field 2 --num_field "
                          "2 --num_vector_local 6 "
                          " --compute_method GPU --num_way 3 --all2all yes"));

  //----------
  //---num_repl, 2-way
  //----------

  char options1[1024];
  char options2[1024];

  char options_template_1[] =
      "--metric_type czekanowski "
      "--num_field 4 --num_vector_local %i --compute_method %s --all2all yes "
      "--num_proc_vector %i --num_proc_repl %i "
      "--num_proc_field %i --num_way %i --num_stage %i";

  for (int gpu=0; gpu<=1; ++gpu) {
    for (int num_vector_local=4; num_vector_local<=5; ++num_vector_local) {
      for (int num_proc_vector=1; num_proc_vector<=6; ++num_proc_vector) {
        for (int num_proc_repl=2; num_proc_repl<=6; ++num_proc_repl) {
          const int num_proc_field = gpu ? 2 : 1;
          if (num_proc_vector * num_proc_field * num_proc_repl > PROCS_MAX) {
            continue;
          }
          const int num_way = 2;
          sprintf(options1, options_template_1,
                  num_vector_local*num_proc_vector,
                  "GPU", 1, 1,
                  1, num_way, 1);
          sprintf(options2, options_template_1, num_vector_local,
                  gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                  num_proc_field, num_way, 1);
          test_2runs(options1, options2);
        }
      }
    }
  }

  //----------
  //---num_repl, num_stage, 3-way
  //----------

  for (int gpu=0; gpu<=1; ++gpu) {
    for (int num_vector_local=6; num_vector_local<=18; num_vector_local+=12) {
      for (int num_proc_vector=1; num_proc_vector<=6; ++num_proc_vector) {
        for (int num_proc_repl=2; num_proc_repl<=6; ++num_proc_repl) {
          for (int num_stage=1; num_stage<=6; num_stage+=4) {
          const int num_proc_field = gpu ? 2 : 1;
            if (num_proc_vector * num_proc_field * num_proc_repl > PROCS_MAX) {
              continue;
            }
            const int num_way = 3;
            sprintf(options1, options_template_1,
                    num_vector_local*num_proc_vector,
                    "GPU", 1, 1,
                   1, num_way, 1);
            sprintf(options2, options_template_1, num_vector_local,
                    gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                    num_proc_field, num_way, num_stage);
            test_2runs(options1, options2);
          }
        }
      }
    }
  }

  //----------
  //---num_phase, 2-way
  //----------

  for (int num_proc_vector=1; num_proc_vector<=8; ++num_proc_vector) {
    for (int num_proc_repl=1; num_proc_repl<=8; ++num_proc_repl) {
      for (int num_phase=2; num_phase<=8; ++num_phase) {
        if (!(num_phase <= 1 + num_proc_vector/2)) {
          continue;
        }
        char options_template[] =
          "--metric_type czekanowski "
          "--num_field 7 --num_vector 12 --compute_method GPU --all2all yes "
          "--num_proc_vector %i --num_proc_repl %i --num_phase %i "
          "--num_way 2";
        sprintf(options1, options_template, num_proc_vector, num_proc_repl,
                num_phase);
        sprintf(options2, options_template, 1, 1, 1);
        test_2runs(options1, options2);
      }
    }
  }

  char options_template_2[] =
      "--metric_type czekanowski "
      "--num_field 4 --num_vector_local %i --compute_method %s --all2all yes "
      "--num_proc_vector %i --num_proc_repl %i "
      "--num_proc_field %i --num_way %i --num_phase %i";

  //----------
  //---num_repl, num_phase, 3-way
  //----------

  for (int gpu=0; gpu<=1; ++gpu) {
    for (int num_vector_local=6; num_vector_local<=6; num_vector_local+=12) {
      for (int num_proc_vector=1; num_proc_vector<=4; ++num_proc_vector) {
        for (int num_proc_repl=1; num_proc_repl<=4; ++num_proc_repl) {
          const int npv = num_proc_vector;
          const int num_phase_max = num_proc_repl==1 ? npv*npv - 2*npv + 2 :
                                                      (npv+1)*(npv+2);
          const int num_phase_min = num_phase_max / 2;
          for (int num_phase=num_phase_min; num_phase<=num_phase_max;
               num_phase+=(num_phase_max-num_phase_min)) {
            if (num_phase < 1) {
              continue;
            }
            const int num_proc_field = gpu ? 2 : 1;
            if (num_proc_vector * num_proc_field * num_proc_repl > PROCS_MAX) {
              continue;
            }
            const int num_way = 3;
            sprintf(options1, options_template_2,
                    num_vector_local*num_proc_vector,
                    "GPU", 1, 1,
                   1, num_way, 1);
            sprintf(options2, options_template_2, num_vector_local,
                    gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                    num_proc_field, num_way, num_phase);
            test_2runs(options1, options2);
          }
        }
      }
    }
  }

  //----------
  //---file output, 2-way
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_proc_vector 1 "
                    "--num_field 7 --num_vector 10 "
                    "--num_way 2 --metric_type czekanowski "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 1 --num_proc_vector 3 "
                    "--num_field 7 --num_vector 10 "
                    "--num_way 2 --metric_type czekanowski "
                    "--compute_method GPU --all2all yes "
                    "--verbosity 1 "
                    "--output_file_stub test_czek_2way"));

  //----------
  //---file output, 3-way
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_proc_vector 1 "
                    "--num_field 7 --num_vector 18 "
                    "--num_way 3 --metric_type czekanowski "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 1 --num_proc_vector 3 "
                    "--num_field 7 --num_vector 18 "
                    "--num_way 3 --metric_type czekanowski "
                    "--compute_method GPU --all2all yes "
                    "--verbosity 1 "
                    "--output_file_stub test_czek_3way"));

  //----------
  //---Misc options
  //----------

  EXPECT_EQ(
      true,
      compare_2runs("--num_proc_vector 1 --num_field 30 --num_vector 3 "
                    "--verbosity 1 --all2all yes "
                    "--compute_method GPU",
                    "--num_proc_vector 1 --num_field 30 --num_vector 3 "
                    "--verbosity 1 --all2all yes "
                    "--compute_method GPU --threshold .65"));

  // TODO: set up better test
  //EXPECT_EQ(
  //    true,
  //    compare_2runs("--num_proc_vector 1 --num_field 3 --num_vector 3 "
  //                  "--compute_method GPU",
  //                  "--num_proc_vector 1 --num_field 3 --num_vector 3 "
  //                  "--compute_method GPU --checksum no"));

  //----------
  //---file input
  //----------

  for (int num_vector = 13; num_vector <= 13; ++num_vector) {
    for (int num_field = 1; num_field <= 10; ++num_field) {
      create_vectors_file("czek_2way_in.bin", num_field, num_vector,
                          GM_METRIC_TYPE_CZEK, 2, problem_type_default(), 1);
      for (int num_proc_vector=1; num_proc_vector<=4; ++num_proc_vector) {
        for (int num_proc_field=1; num_proc_field<=5; ++num_proc_field) {

          char options1[1024];
          char options2[1024];

          char options_template[] =
                 "--num_vector %i --num_field %i "
                 "--num_proc_vector %i --num_proc_field %i "
                 "--compute_method GPU --all2all yes %s --verbosity 1";

          sprintf(options1, options_template,
                  num_vector, num_field, num_proc_vector, num_proc_field, "");

          sprintf(options2, options_template,
                  num_vector, num_field, num_proc_vector, num_proc_field,
                  "--input_file czek_2way_in.bin");

          test_2runs(options1, options2);
        }
      }
    }
  }
} // DriverTest_czek_

//=============================================================================

void DriverTest_ccc2_simple_compute_method(int compute_method) {
  const int num_field = 5;
  const int num_vector_local = 2;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env, MPI_COMM_WORLD, 0, NULL, NULL);
  env->metric_type_ = GM_METRIC_TYPE_CCC;
  env->num_way_ = 2;
  env->all2all_ = false;
  GMEnv_set_compute_method(env, compute_method);
  GMEnv_set_num_proc(env, 1, 1, 1);

  GMDecompMgr dm_value = GMDecompMgr_null(), *dm = &dm_value;
  GMDecompMgr_create(dm, true, true, num_field, num_vector_local,
                     GMEnv_data_type_vectors(env), env);

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), dm, env);

  if (GMEnv_is_proc_active(env)) {
    {
      const int G = 0;
      const int T = 1;
      int f = 0;
      const int i = 0;
      GMVectors_bits2_set(vectors, f++, i, 2 * G + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
    }
    {
      const int G = 0;
      const int A = 1;
      int f = 0;
      const int i = 1;
      GMVectors_bits2_set(vectors, f++, i, 2 * G + 1 * G, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * G, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * G + 1 * G, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * G + 1 * G, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * G + 1 * A, env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  GMMetricsMem metrics_mem(env);
  GMMetrics_create(metrics, GMEnv_data_type_metrics(env), dm,
                   &metrics_mem, env);

  gm_compute_metrics(metrics, vectors, env);

  if (GMEnv_is_proc_active(env)) {
    const double result00 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 0, 0, env);
    const double result01 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 0, 1, env);
    const double result10 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 1, 0, env);
    const double result11 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 1, 1, env);

    printf("G G  %.5f\n", result00);
    printf("G A  %.5f\n", result01);
    printf("T G  %.5f\n", result10);
    printf("T A  %.5f\n", result11);
    printf("\n");

    const double ref00 = .196;
    const double ref01 = .000;
    const double ref10 = .588;
    const double ref11 = .312;

    printf("G G  %.5f\n", ref00);
    printf("G A  %.5f\n", ref01);
    printf("T G  %.5f\n", ref10);
    printf("T A  %.5f\n", ref11);
    printf("\n");

    const double eps = 1.e-5;

    EXPECT_EQ(true, fabs(result00 - ref00) < eps);
    EXPECT_EQ(true, fabs(result01 - ref01) < eps);
    EXPECT_EQ(true, fabs(result10 - ref10) < eps);
    EXPECT_EQ(true, fabs(result11 - ref11) < eps);
  }

  GMMetrics_destroy(metrics, env);
  GMVectors_destroy(vectors, env);
  GMDecompMgr_destroy(dm, env);
  GMEnv_destroy(env);
} // DriverTest_ccc2_simple_compute_method

//=============================================================================

void DriverTest_ccc2_simple_() {
  DriverTest_ccc2_simple_compute_method(GM_COMPUTE_METHOD_REF);
  DriverTest_ccc2_simple_compute_method(GM_COMPUTE_METHOD_CPU);
  DriverTest_ccc2_simple_compute_method(GM_COMPUTE_METHOD_GPU);
}

//=============================================================================

void DriverTest_ccc2_simple_sparse_compute_method(int compute_method) {
  const int num_field = 5;
  const int num_vector_local = 2;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env, MPI_COMM_WORLD, 0, NULL, NULL);
  env->metric_type_ = GM_METRIC_TYPE_CCC;
  env->num_way_ = 2;
  env->all2all_ = false;
  GMEnv_set_compute_method(env, compute_method);
  GMEnv_set_num_proc(env, 1, 1, 1);
  env->sparse = true;

  GMDecompMgr dm_value = GMDecompMgr_null(), *dm = &dm_value;
  GMDecompMgr_create(dm, true, true, num_field, num_vector_local,
                     GMEnv_data_type_vectors(env), env);

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), dm, env);

  if (GMEnv_is_proc_active(env)) {
    const int UN = 2 * 1 + 1 * 0;
    {
      const int G = 0;
      const int T = 1;
    //const int GG =  2 * G + 1 * G;
      const int GT =  2 * G + 1 * T;
    //const int TG =  2 * G + 1 * T;
      const int TT =  2 * T + 1 * T;
      int f = 0;
      const int i = 0;
      GMVectors_bits2_set(vectors, f++, i, GT, env);
      GMVectors_bits2_set(vectors, f++, i, TT, env);
      GMVectors_bits2_set(vectors, f++, i, TT, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
    }
    {
      const int G = 0;
      const int A = 1;
      const int GG =  2 * G + 1 * G;
      const int GA =  2 * G + 1 * A;
      const int AG =  2 * G + 1 * A;
    //const int AA =  2 * Q + 1 * A;
      int f = 0;
      const int i = 1;
      GMVectors_bits2_set(vectors, f++, i, GG, env);
      GMVectors_bits2_set(vectors, f++, i, AG, env);
      GMVectors_bits2_set(vectors, f++, i, GG, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, GA, env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  GMMetricsMem metrics_mem(env);
  GMMetrics_create(metrics, GMEnv_data_type_metrics(env), dm,
                   &metrics_mem, env);

  gm_compute_metrics(metrics, vectors, env);

  if (GMEnv_is_proc_active(env)) {
    const double result00 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 0, 0, env);
    const double result01 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 0, 1, env);
    const double result10 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 1, 0, env);
    const double result11 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 1, 1, env);

    printf("G G  %.5f\n", result00);
    printf("G A  %.5f\n", result01);
    printf("T G  %.5f\n", result10);
    printf("T A  %.5f\n", result11);
    printf("\n");

    const double s0_0 = 1;
    const double s0_1 = 5;

    const double s1_0 = 6;
    const double s1_1 = 2;

    const double c0 = s0_0 + s0_1;
    const double c1 = s1_0 + s1_1;

    const double f0_0 = s0_0 / c0;
    const double f0_1 = s0_1 / c0;

    const double f1_0 = s1_0 / c1;
    const double f1_1 = s1_1 / c1;

    const double r0_00 = 2;
    const double r0_01 = 0;
    const double r0_10 = 2;
    const double r0_11 = 0;

    const double r1_00 = 0;
    const double r1_01 = 0;
    const double r1_10 = 2;
    const double r1_11 = 2;

    const double r2_00 = 0;
    const double r2_01 = 0;
    const double r2_10 = 4;
    const double r2_11 = 0;

    const double r_00 = r0_00 + r1_00 + r2_00;
    const double r_01 = r0_01 + r1_01 + r2_01;
    const double r_10 = r0_10 + r1_10 + r2_10;
    const double r_11 = r0_11 + r1_11 + r2_11;

    const double c = r_00 + r_01 + r_10 + r_11;

    const double f_00 = r_00 / c;
    const double f_01 = r_01 / c;
    const double f_10 = r_10 / c;
    const double f_11 = r_11 / c;

    const double fm = 9 / (double) 2;
    const double cp = 2 / (double) 3;

    const double ref00 = fm * f_00 * ( 1 - cp * f0_0 ) * ( 1 - cp * f1_0 );
    const double ref01 = fm * f_01 * ( 1 - cp * f0_0 ) * ( 1 - cp * f1_1 );
    const double ref10 = fm * f_10 * ( 1 - cp * f0_1 ) * ( 1 - cp * f1_0 );
    const double ref11 = fm * f_11 * ( 1 - cp * f0_1 ) * ( 1 - cp * f1_1 );

    printf("G G  %.5f\n", ref00);
    printf("G A  %.5f\n", ref01);
    printf("T G  %.5f\n", ref10);
    printf("T A  %.5f\n", ref11);
    printf("\n");

    const double eps = 1.e-5;

    EXPECT_EQ(true, fabs(result00 - ref00) < eps);
    EXPECT_EQ(true, fabs(result01 - ref01) < eps);
    EXPECT_EQ(true, fabs(result10 - ref10) < eps);
    EXPECT_EQ(true, fabs(result11 - ref11) < eps);
  }

  GMMetrics_destroy(metrics, env);
  GMVectors_destroy(vectors, env);
  GMDecompMgr_destroy(dm, env);
  GMEnv_destroy(env);
} // DriverTest_ccc2_simple_sparse_compute_method

//=============================================================================

void DriverTest_ccc2_simple_sparse_() {
  DriverTest_ccc2_simple_sparse_compute_method(GM_COMPUTE_METHOD_REF);
  DriverTest_ccc2_simple_sparse_compute_method(GM_COMPUTE_METHOD_CPU);
  DriverTest_ccc2_simple_sparse_compute_method(GM_COMPUTE_METHOD_GPU);
}

//=============================================================================

void DriverTest_duo2_simple_sparse_compute_method(int compute_method) {
  const int num_field = 5;
  const int num_vector_local = 2;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env, MPI_COMM_WORLD, 0, NULL, NULL);
  env->metric_type_ = GM_METRIC_TYPE_DUO;
  env->num_way_ = 2;
  env->all2all_ = false;
  GMEnv_set_compute_method(env, compute_method);
  GMEnv_set_num_proc(env, 1, 1, 1);
  env->sparse = true;

  GMDecompMgr dm_value = GMDecompMgr_null(), *dm = &dm_value;
  GMDecompMgr_create(dm, true, true, num_field, num_vector_local,
                     GMEnv_data_type_vectors(env), env);

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), dm, env);

  if (GMEnv_is_proc_active(env)) {
    // entry choices, binary representation
    const int MIN = 2 * (0) + 1 * (0);
    const int MAX = 2 * (1) + 1 * (1);
    const int UNK = 2 * (1) + 1 * (0);
    // define first vector
    {
      int f = 0;
      const int i = 0;
      GMVectors_bits2_set(vectors, f++, i, MIN, env);
      GMVectors_bits2_set(vectors, f++, i, MAX, env);
      GMVectors_bits2_set(vectors, f++, i, MIN, env);
      GMVectors_bits2_set(vectors, f++, i, UNK, env);
      GMVectors_bits2_set(vectors, f++, i, UNK, env);
    }
    // define second vector
    {
      int f = 0;
      const int i = 1;
      GMVectors_bits2_set(vectors, f++, i, UNK, env);
      GMVectors_bits2_set(vectors, f++, i, MAX, env);
      GMVectors_bits2_set(vectors, f++, i, MAX, env);
      GMVectors_bits2_set(vectors, f++, i, MAX, env);
      GMVectors_bits2_set(vectors, f++, i, UNK, env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  GMMetricsMem metrics_mem(env);
  GMMetrics_create(metrics, GMEnv_data_type_metrics(env), dm,
                   &metrics_mem, env);

  gm_compute_metrics(metrics, vectors, env);

  if (GMEnv_is_proc_active(env)) {
    const double result00 =
        GMMetrics_duo_get_from_index_2(metrics, 0, 0, 0, env);
    const double result01 =
        GMMetrics_duo_get_from_index_2(metrics, 0, 0, 1, env);
    const double result10 =
        GMMetrics_duo_get_from_index_2(metrics, 0, 1, 0, env);
    const double result11 =
        GMMetrics_duo_get_from_index_2(metrics, 0, 1, 1, env);

    printf("COMPUTED: MIN MIN  %.5f\n", result00);
    printf("COMPUTED: MIN MAX  %.5f\n", result01);
    printf("COMPUTED: MAX MIN  %.5f\n", result10);
    printf("COMPUTED: MAX MAX  %.5f\n", result11);
    printf("\n");

    // calculate by hand the expected DUO value.

    // recall that num_field = 5;

    const double s0_0 = 2; // vector 0, number of MINs
    const double s0_1 = 1; // vector 0, number of MAXs

    const double s1_0 = 0;
    const double s1_1 = 3;

    const double c0 = s0_0 + s0_1; // vector 0, number of MINs and MAXs
    const double c1 = s1_0 + s1_1;

    // const double unk_0 = num_field - c0; // = 2 // for vector 0, number of UNK
    // const double unk_1 = num_field - c1; // = 2 // for vector 1, number of UNK

    const double f0_0 = s0_0 / c0; // vector 0, f_0(MIN)
    const double f0_1 = s0_1 / c0; // vector 0, f_0(MAX)

    const double f1_0 = s1_0 / c1;
    const double f1_1 = s1_1 / c1;

    // numerators from table computation
    const double r_00 = 0;
    const double r_01 = 1;
    const double r_10 = 0;
    const double r_11 = 1;

    // sum of all table entries
    const double c = r_00 + r_01 + r_10 + r_11;

    // Dij of DUO
    const double d_00 = r_00 / c;
    const double d_01 = r_01 / c;
    const double d_10 = r_10 / c;
    const double d_11 = r_11 / c;

    // Constants needed by method
    const double fm = 4;
    const double cp = 2 / (double) 3;

    // DUO values
    const double ref00 = fm * d_00 * ( 1 - cp * f0_0 ) * ( 1 - cp * f1_0 );
    const double ref01 = fm * d_01 * ( 1 - cp * f0_0 ) * ( 1 - cp * f1_1 );
    const double ref10 = fm * d_10 * ( 1 - cp * f0_1 ) * ( 1 - cp * f1_0 );
    const double ref11 = fm * d_11 * ( 1 - cp * f0_1 ) * ( 1 - cp * f1_1 );

    printf("EXPECTED: MIN MIN  %.5f\n", ref00);
    printf("EXPECTED: MIN MAX  %.5f\n", ref01);
    printf("EXPECTED: MAX MIN  %.5f\n", ref10);
    printf("EXPECTED: MAX MAX  %.5f\n", ref11);
    printf("\n");

    const double eps = 1.e-5;

    EXPECT_EQ(true, fabs(result00 - ref00) < eps);
    EXPECT_EQ(true, fabs(result01 - ref01) < eps);
    EXPECT_EQ(true, fabs(result10 - ref10) < eps);
    EXPECT_EQ(true, fabs(result11 - ref11) < eps);
  }

  GMMetrics_destroy(metrics, env);
  GMVectors_destroy(vectors, env);
  GMDecompMgr_destroy(dm, env);
  GMEnv_destroy(env);
} // DriverTest_duo2_simple_sparse_compute_method

//=============================================================================

void DriverTest_duo2_simple_sparse_() {
  DriverTest_duo2_simple_sparse_compute_method(GM_COMPUTE_METHOD_REF);
  DriverTest_duo2_simple_sparse_compute_method(GM_COMPUTE_METHOD_CPU);
  DriverTest_duo2_simple_sparse_compute_method(GM_COMPUTE_METHOD_GPU);
}

//=============================================================================

void DriverTest_ccc3_simple_compute_method(int compute_method) {
  const int num_field = 10;
  const int num_vector_local = 3;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env, MPI_COMM_WORLD, 0, NULL, NULL);
  env->metric_type_ = GM_METRIC_TYPE_CCC;
  env->num_way_ = 3;
  env->all2all_ = true;
  GMEnv_set_compute_method(env, compute_method);
  GMEnv_set_num_proc(env, 1, 1, 1);

  GMDecompMgr dm_value = GMDecompMgr_null(), *dm = &dm_value;
  GMDecompMgr_create(dm, true, true, num_field, num_vector_local,
                     GMEnv_data_type_vectors(env), env);

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), dm, env);

  if (GMEnv_is_proc_active(env)) {
    {
      const int A = 0;
      const int T = 1;
      int f = 0;
      const int i = 0;
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
    }
    {
      const int A = 0;
      const int T = 1;
      int f = 0;
      const int i = 1;
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
    }
    {
      const int A = 0;
      const int T = 1;
      int f = 0;
      const int i = 2;
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  GMMetricsMem metrics_mem(env);
  GMMetrics_create(metrics, GMEnv_data_type_metrics(env), dm,
                   &metrics_mem, env);

  gm_compute_metrics(metrics, vectors, env);

  if (GMEnv_is_proc_active(env)) {
    const double result000 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 0, 0, env);
    const double result001 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 0, 1, env);
    const double result010 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 1, 0, env);
    const double result011 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 1, 1, env);
    const double result100 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 0, 0, env);
    const double result101 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 0, 1, env);
    const double result110 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 1, 0, env);
    const double result111 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 1, 1, env);

    printf("A A A  %.8f\n", result000);
    printf("A A T  %.5f\n", result001);
    printf("A T A  %.8f\n", result010);
    printf("A T T  %.8f\n", result011);
    printf("T A A  %.8f\n", result100);
    printf("T A T  %.8f\n", result101);
    printf("T T A  %.8f\n", result110);
    printf("T T T  %.8f\n", result111);
    printf("\n");

    const double fm = 9 / (double) 2;
    const double ref000 = fm * .055;
    const double ref001 = fm * .016;
    const double ref010 = fm * .016;
    const double ref011 = fm * .030;
    const double ref100 = fm * .039;
    const double ref101 = fm * .008;
    const double ref110 = fm * .008;
    const double ref111 = fm * .015;

    printf("A A A  %.8f\n", ref000);
    printf("A A T  %.5f\n", ref001);
    printf("A T A  %.8f\n", ref010);
    printf("A T T  %.8f\n", ref011);
    printf("T A A  %.8f\n", ref100);
    printf("T A T  %.8f\n", ref101);
    printf("T T A  %.8f\n", ref110);
    printf("T T T  %.8f\n", ref111);
    printf("\n");

    const double eps = 1.e-3;

    EXPECT_EQ(true, fabs(result000 - ref000) < eps);
    EXPECT_EQ(true, fabs(result001 - ref001) < eps);
    EXPECT_EQ(true, fabs(result010 - ref010) < eps);
    EXPECT_EQ(true, fabs(result011 - ref011) < eps);
    EXPECT_EQ(true, fabs(result100 - ref100) < eps);
    EXPECT_EQ(true, fabs(result101 - ref101) < eps);
    EXPECT_EQ(true, fabs(result110 - ref110) < eps);
    EXPECT_EQ(true, fabs(result111 - ref111) < eps);
  }

  GMMetrics_destroy(metrics, env);
  GMVectors_destroy(vectors, env);
  GMDecompMgr_destroy(dm, env);
  GMEnv_destroy(env);
} // DriverTest_ccc3_simple_compute_method

//=============================================================================

void DriverTest_ccc3_simple_() {
  DriverTest_ccc3_simple_compute_method(GM_COMPUTE_METHOD_REF);
  DriverTest_ccc3_simple_compute_method(GM_COMPUTE_METHOD_CPU);
  DriverTest_ccc3_simple_compute_method(GM_COMPUTE_METHOD_GPU);
}

//=============================================================================

void DriverTest_ccc3_simple_sparse_compute_method(int compute_method) {
  const int num_field = 10;
  const int num_vector_local = 3;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env, MPI_COMM_WORLD, 0, NULL, NULL);
  env->metric_type_ = GM_METRIC_TYPE_CCC;
  env->num_way_ = 3;
  env->all2all_ = true;
  GMEnv_set_compute_method(env, compute_method);
  GMEnv_set_num_proc(env, 1, 1, 1);
  env->sparse = true;

  GMDecompMgr dm_value = GMDecompMgr_null(), *dm = &dm_value;
  GMDecompMgr_create(dm, true, true, num_field, num_vector_local,
                     GMEnv_data_type_vectors(env), env);

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), dm, env);

  if (GMEnv_is_proc_active(env)) {
    const int UN = 2 * 1 + 1 * 0;
    {
      const int A = 0;
      const int T = 1;
      const int AA =  2 * A + 1 * A;
      const int AT =  2 * A + 1 * T;
    //const int TA =  2 * A + 1 * T;
      const int TT =  2 * T + 1 * T;
      int f = 0;
      const int i = 0;
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, TT, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
    }
    {
      const int A = 0;
      const int T = 1;
      const int AA =  2 * A + 1 * A;
      const int AT =  2 * A + 1 * T;
    //const int TA =  2 * A + 1 * T;
    //const int TT =  2 * T + 1 * T;
      int f = 0;
      const int i = 1;
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
    }
    {
      const int A = 0;
      const int T = 1;
      const int AA =  2 * A + 1 * A;
      const int AT =  2 * A + 1 * T;
    //const int TA =  2 * A + 1 * T;
    //const int TT =  2 * T + 1 * T;
      int f = 0;
      const int i = 2;
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  GMMetricsMem metrics_mem(env);
  GMMetrics_create(metrics, GMEnv_data_type_metrics(env), dm,
                   &metrics_mem, env);

  gm_compute_metrics(metrics, vectors, env);

  if (GMEnv_is_proc_active(env)) {
    const double result000 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 0, 0, env);
    const double result001 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 0, 1, env);
    const double result010 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 1, 0, env);
    const double result011 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 1, 1, env);
    const double result100 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 0, 0, env);
    const double result101 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 0, 1, env);
    const double result110 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 1, 0, env);
    const double result111 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 1, 1, env);

    printf("A A A  %.8f\n", result000);
    printf("A A T  %.5f\n", result001);
    printf("A T A  %.8f\n", result010);
    printf("A T T  %.8f\n", result011);
    printf("T A A  %.8f\n", result100);
    printf("T A T  %.8f\n", result101);
    printf("T T A  %.8f\n", result110);
    printf("T T T  %.8f\n", result111);
    printf("\n");

    const double s0_0 = 9;
    const double s0_1 = 5;
    const double c0 = s0_0 + s0_1;

    const double s1_0 = 12;
    const double s1_1 = 4;
    const double c1 = s1_0 + s1_1;

    const double s2_0 = 14;
    const double s2_1 = 2;
    const double c2 = s2_0 + s2_1;

    const double f0_0 = s0_0 / c0;
    const double f0_1 = s0_1 / c0;
    const double f1_0 = s1_0 / c1;
    const double f1_1 = s1_1 / c1;
    const double f2_0 = s2_0 / c2;
    const double f2_1 = s2_1 / c2;

    const double r0_000 = 4;
    const double r0_001 = 0;
    const double r0_010 = 4;
    const double r0_011 = 0;
    const double r0_100 = 0;
    const double r0_101 = 0;
    const double r0_110 = 0;
    const double r0_111 = 0;

    const double r1_000 = 4;
    const double r1_001 = 0;
    const double r1_010 = 0;
    const double r1_011 = 0;
    const double r1_100 = 4;
    const double r1_101 = 0;
    const double r1_110 = 0;
    const double r1_111 = 0;

    const double r2_000 = 0;
    const double r2_001 = 0;
    const double r2_010 = 0;
    const double r2_011 = 0;
    const double r2_100 = 8;
    const double r2_101 = 0;
    const double r2_110 = 0;
    const double r2_111 = 0;

    const double r4_000 = 1;
    const double r4_001 = 1;
    const double r4_010 = 1;
    const double r4_011 = 1;
    const double r4_100 = 1;
    const double r4_101 = 1;
    const double r4_110 = 1;
    const double r4_111 = 1;

    const double r5_000 = 4;
    const double r5_001 = 0;
    const double r5_010 = 0;
    const double r5_011 = 0;
    const double r5_100 = 4;
    const double r5_101 = 0;
    const double r5_110 = 0;
    const double r5_111 = 0;

    const double r_000 = r0_000 + r1_000 + r2_000 + r4_000 + r5_000;
    const double r_001 = r0_001 + r1_001 + r2_001 + r4_001 + r5_001;
    const double r_010 = r0_010 + r1_010 + r2_010 + r4_010 + r5_010;
    const double r_011 = r0_011 + r1_011 + r2_011 + r4_011 + r5_011;
    const double r_100 = r0_100 + r1_100 + r2_100 + r4_100 + r5_100;
    const double r_101 = r0_101 + r1_101 + r2_101 + r4_101 + r5_101;
    const double r_110 = r0_110 + r1_110 + r2_110 + r4_110 + r5_110;
    const double r_111 = r0_111 + r1_111 + r2_111 + r4_111 + r5_111;

    const double c = r_000 + r_001 + r_010 + r_011 +
                     r_100 + r_101 + r_110 + r_111;

    const double f_000 = r_000 / c;
    const double f_001 = r_001 / c;
    const double f_010 = r_010 / c;
    const double f_011 = r_011 / c;
    const double f_100 = r_100 / c;
    const double f_101 = r_101 / c;
    const double f_110 = r_110 / c;
    const double f_111 = r_111 / c;

    const double fm = 9 / (double) 2;
    const double cp = 2 / (double) 3;

    const double ref000 = fm * f_000 * (1-cp*f0_0) * (1-cp*f1_0) * (1-cp*f2_0);
    const double ref001 = fm * f_001 * (1-cp*f0_0) * (1-cp*f1_0) * (1-cp*f2_1);
    const double ref010 = fm * f_010 * (1-cp*f0_0) * (1-cp*f1_1) * (1-cp*f2_0);
    const double ref011 = fm * f_011 * (1-cp*f0_0) * (1-cp*f1_1) * (1-cp*f2_1);
    const double ref100 = fm * f_100 * (1-cp*f0_1) * (1-cp*f1_0) * (1-cp*f2_0);
    const double ref101 = fm * f_101 * (1-cp*f0_1) * (1-cp*f1_0) * (1-cp*f2_1);
    const double ref110 = fm * f_110 * (1-cp*f0_1) * (1-cp*f1_1) * (1-cp*f2_0);
    const double ref111 = fm * f_111 * (1-cp*f0_1) * (1-cp*f1_1) * (1-cp*f2_1);

    printf("A A A  %.8f\n", ref000);
    printf("A A T  %.5f\n", ref001);
    printf("A T A  %.8f\n", ref010);
    printf("A T T  %.8f\n", ref011);
    printf("T A A  %.8f\n", ref100);
    printf("T A T  %.8f\n", ref101);
    printf("T T A  %.8f\n", ref110);
    printf("T T T  %.8f\n", ref111);
    printf("\n");

    //const double ref000 = .055;
    //const double ref001 = .016;
    //const double ref010 = .016;
    //const double ref011 = .030;
    //const double ref100 = .039;
    //const double ref101 = .008;
    //const double ref110 = .008;
    //const double ref111 = .015;

    const double eps = 1.e-3;

    EXPECT_EQ(true, fabs(result000 - ref000) < eps);
    EXPECT_EQ(true, fabs(result001 - ref001) < eps);
    EXPECT_EQ(true, fabs(result010 - ref010) < eps);
    EXPECT_EQ(true, fabs(result011 - ref011) < eps);
    EXPECT_EQ(true, fabs(result100 - ref100) < eps);
    EXPECT_EQ(true, fabs(result101 - ref101) < eps);
    EXPECT_EQ(true, fabs(result110 - ref110) < eps);
    EXPECT_EQ(true, fabs(result111 - ref111) < eps);
  }

  GMMetrics_destroy(metrics, env);
  GMVectors_destroy(vectors, env);
  GMDecompMgr_destroy(dm, env);
  GMEnv_destroy(env);
} // DriverTest_ccc3_simple_sparse_compute_method

//=============================================================================

void DriverTest_ccc3_simple_sparse_() {
  DriverTest_ccc3_simple_sparse_compute_method(GM_COMPUTE_METHOD_REF);
  DriverTest_ccc3_simple_sparse_compute_method(GM_COMPUTE_METHOD_CPU);
  DriverTest_ccc3_simple_sparse_compute_method(GM_COMPUTE_METHOD_GPU);
}

//=============================================================================

void DriverTest_ccc_duo_(const char* const metric_type) {
  const bool is_duo = metric_type[0] == 'd';

  {
    char options1[1024];
    char options2[1024];

    char options_template_tc[] =
        //"--num_proc_vector 1 --num_field 1 --num_vector_local 5 "
        "--metric_type %s "
        "--num_proc_vector 2 --num_field 100 --num_vector %i "
        "--compute_method %s --sparse %s "
        "--problem_type random --verbosity %i --tc %i --num_way %i "
        "--num_tc_steps %i --all2all yes" ;

    for (int num_tc_steps=1; num_tc_steps<=3; ++num_tc_steps) {
    for (int nvl=3; nvl<=10; ++nvl) {
    for (int num_way=2; num_way<=3; ++num_way) {
      if (is_duo && num_way == 3) continue;
    for (int sparse=0; sparse<=1; ++sparse) {
      if (is_duo && sparse == 0) continue;
    for (int tc=1; tc<=2; ++tc) {
      if (!gm_is_tc_valid(tc)) continue;
      if (nvl/2 < num_way) continue;

      sprintf(options1, options_template_tc, metric_type, nvl, "REF",
              sparse==0 ? "no" : "yes", 1, 0, num_way, 1);
      sprintf(options2, options_template_tc, metric_type, nvl, "GPU",
              sparse==0 ? "no" : "yes", 1, tc, num_way, num_tc_steps);
      EXPECT_EQ(true, compare_2runs(options1, options2));
    }
    }
    }
    }
    }
  }

//FIX
#if 1
  char options1[1024];
  char options2[1024];
  char options3[1024];

  //----------
  //---2-way, all2all no
  //----------

  char options_template_1[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      //"--problem_type analytic " 
      "--num_proc_vector 1 --num_field %i --num_vector_local %i "
      "--compute_method %s";

  sprintf(options1, options_template_1, metric_type, 2, 1, 2, "REF");
  sprintf(options2, options_template_1, metric_type, 2, 1, 2, "CPU");
  sprintf(options3, options_template_1, metric_type, 2, 1, 2, "GPU");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_1, metric_type, 1, 100, 2, "REF");
  sprintf(options2, options_template_1, metric_type, 1, 100, 2, "CPU");
  sprintf(options3, options_template_1, metric_type, 1, 100, 2, "GPU");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_1, metric_type, 1, 100, 48, "REF");
  sprintf(options2, options_template_1, metric_type, 1, 100, 48, "CPU");
  sprintf(options3, options_template_1, metric_type, 1, 100, 48, "GPU");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  for (int i = 1; i <= 100; ++i) {
    if ((i+3) % 32 > 6) {
      continue;
    }
    sprintf(options1, options_template_1, metric_type, 1, i, 48, "REF");
    sprintf(options2, options_template_1, metric_type, 1, i, 48, "CPU");
    sprintf(options3, options_template_1, metric_type, 1, i, 48, "GPU");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
  }

  //----------
  //---2-way, all2all yes, small
  //----------

  char options_template_2[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      //"--problem_type analytic "
      "--num_proc_vector %i --num_field %i --num_vector_local %i "
      "--compute_method %s --all2all %s";

  sprintf(options1, options_template_2, metric_type, 2, 1, 1, 2,
          "REF", "no");
  sprintf(options2, options_template_2, metric_type, 2, 1, 1, 2,
          "CPU", "yes");
  sprintf(options3, options_template_2, metric_type, 2, 1, 1, 2,
          "GPU", "yes");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_2, metric_type, 1, 1, 1, 4,
          "REF", "no");
  sprintf(options2, options_template_2, metric_type, 1, 2, 1, 2,
          "CPU", "yes");
  sprintf(options3, options_template_2, metric_type, 1, 2, 1, 2,
          "GPU", "yes");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  char options_template_2a[] =
      "--sparse yes "
      "--num_proc_vector 1 --num_field 2 --num_vector 5 "
      //"--problem_type analytic "
      "--compute_method REF --all2all yes --metric_type %s";

  char options_template_2b[] =
      "--sparse yes "
      "--num_proc_vector 2 --num_field 2 --num_vector 5 "
      //"--problem_type analytic "
      "--compute_method CPU --all2all yes --metric_type %s";

  char options_template_2c[] =
      "--sparse yes "
      "--num_proc_vector 2 --num_field 2 --num_vector 5 "
      //"--problem_type analytic "
      "--compute_method GPU --all2all yes --metric_type %s";

  sprintf(options1, options_template_2a, metric_type);
  sprintf(options2, options_template_2b, metric_type);
  sprintf(options3, options_template_2c, metric_type);
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  char options_template_2a2[] =
      "--sparse yes "
      "--num_proc_vector 1 --num_proc_field 1 "
      "--num_field 7 --num_vector 2 "
      //"--problem_type analytic "
      "--compute_method REF --all2all yes --metric_type %s";

  char options_template_2b2[] =
      "--sparse yes "
      "--num_proc_vector 1 --num_proc_field 3 "
      "--num_field 7 --num_vector 2 "
      //"--problem_type analytic "
      "--compute_method GPU --all2all yes --metric_type %s";

  sprintf(options1, options_template_2a2, metric_type);
  sprintf(options2, options_template_2b2, metric_type);
  EXPECT_EQ(true, compare_2runs(options1, options2));

  //----------
  //---2-way, all2all yes, large
  //----------

  sprintf(options1, options_template_2, metric_type, 1, 1, 100, 48,
          "REF", "no");
  sprintf(options2, options_template_2, metric_type, 1, 1, 100, 48,
          "CPU", "yes");
  sprintf(options3, options_template_2, metric_type, 1, 1, 100, 48,
          "GPU", "yes");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_2, metric_type, 1, 1, 100, 48,
          "REF", "no");
  sprintf(options2, options_template_2, metric_type, 1, 2, 100, 24,
          "CPU", "yes");
  sprintf(options3, options_template_2, metric_type, 1, 2, 100, 24,
          "GPU", "yes");
  EXPECT_EQ(true, compare_3runs(options1, options2, options3));

//FIX
//#endif
//#if 1

  //----------
  //---3-way, all2all no
  //----------

  char options_template_3[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      "--num_proc_vector 1 --num_field %i --num_vector_local %i "
      //"--problem_type analytic "
      "--compute_method %s --num_way 3";

  if (!is_duo) {
    sprintf(options1, options_template_3, metric_type, 2, 1, 3, "REF");
    sprintf(options2, options_template_3, metric_type, 2, 1, 3, "CPU");
    sprintf(options3, options_template_3, metric_type, 2, 1, 3, "GPU");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    sprintf(options1, options_template_3, metric_type, 1, 100, 48, "REF");
    sprintf(options2, options_template_3, metric_type, 1, 100, 48, "CPU");
    sprintf(options3, options_template_3, metric_type, 1, 100, 48, "GPU");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    for (int i = 1; i <= 100; ++i) {
      if ((i+3) % 32 > 6) {
        continue;
      }
      sprintf(options1, options_template_3, metric_type, 0, i, 24, "REF");
      sprintf(options2, options_template_3, metric_type, 0, i, 24, "CPU");
      sprintf(options3, options_template_3, metric_type, 0, i, 24, "GPU");
      EXPECT_EQ(true, compare_3runs(options1, options2, options3));
    }
  }

  //----------
  //---3-way, all2all yes, small
  //----------

  char options_template_4[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      "--num_proc_vector %i --num_field %i --num_vector_local %i "
      //"--problem_type analytic "
      "--compute_method %s --all2all %s --num_way 3";

  if (!is_duo) {
    sprintf(options1, options_template_4, metric_type, 2, 1, 1, 3, "REF", "no");
    sprintf(options2, options_template_4, metric_type, 2, 1, 1, 3, "CPU", "yes");
    sprintf(options3, options_template_4, metric_type, 2, 1, 1, 3, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    sprintf(options1, options_template_4, metric_type, 1, 1, 1, 6, "REF", "no");
    sprintf(options2, options_template_4, metric_type, 1, 2, 1, 3, "CPU", "yes");
    sprintf(options3, options_template_4, metric_type, 1, 2, 1, 3, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
  }

  //----------
  //---3-way, all2all yes, large
  //----------

  if (!is_duo) {
    sprintf(options1, options_template_4, metric_type, 1, 1, 100, 48, "REF", "no");
    sprintf(options2, options_template_4, metric_type, 1, 1, 100, 48, "CPU", "yes");
    sprintf(options3, options_template_4, metric_type, 1, 1, 100, 48, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    sprintf(options1, options_template_4, metric_type, 1, 1, 100, 48, "REF", "no");
    sprintf(options2, options_template_4, metric_type, 1, 4, 100, 12, "CPU", "yes");
    sprintf(options3, options_template_4, metric_type, 1, 4, 100, 12, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

  //  sprintf(options1, options_template_4, metric_type, 1, 3, 1, 6, "REF", "yes");
  //  sprintf(options2, options_template_4, metric_type, 1, 1, 1, 18, "GPU", "yes");
  //  sprintf(options3, options_template_4, metric_type, 1, 3, 1, 6, "GPU", "yes");
  //  EXPECT_EQ(true, compare_3runs(options1, options2, options3));
    //EXPECT_EQ(true, compare_2runs(options1, options2));

    char options_template_4a[] =
        "--sparse yes "
        "--num_proc_vector 1 --num_field 2 --num_vector 16 "
        //"--problem_type analytic "
        "--compute_method REF --all2all yes --num_way 3 "
        "--metric_type %s";

    char options_template_4b[] =
        "--sparse yes "
        "--num_proc_vector 3 --num_field 2 --num_vector 16 "
        //"--problem_type analytic "
        "--compute_method CPU --all2all yes --num_way 3 "
        "--metric_type %s";

    char options_template_4c[] =
        "--sparse yes "
        "--num_proc_vector 3 --num_field 2 --num_vector 16 "
        //"--problem_type analytic "
        "--compute_method GPU --all2all yes --num_way 3 "
        "--metric_type %s";

    sprintf(options1, options_template_4a, metric_type);
    sprintf(options2, options_template_4b, metric_type);
    sprintf(options3, options_template_4c, metric_type);
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    char options_template_4a2[] =
        "--sparse yes "
        "--num_proc_vector 1 --num_proc_field 1 "
        //"--problem_type analytic "
        "--num_field 13 --num_vector 3 "
        "--compute_method REF --all2all yes --num_way 3 "
        "--metric_type %s";

    char options_template_4b2[] =
        "--sparse yes "
        "--num_proc_vector 1 --num_proc_field 5 "
        //"--problem_type analytic "
        "--num_field 13 --num_vector 3 "
        "--compute_method GPU --all2all yes --num_way 3 "
        "--metric_type %s";

    sprintf(options1, options_template_4a2, metric_type);
    sprintf(options2, options_template_4b2, metric_type);
    EXPECT_EQ(true, compare_2runs(options1, options2));
  }

  //----------
  //---num_proc_field
  //----------

  char options_template_5[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      "--num_proc_vector %i --num_proc_field %i "
      " --num_field %i --num_vector_local %i "
      //"--problem_type analytic "
      "--compute_method %s --all2all %s";

  for (int i = 1; i <= 3; ++i) {
    sprintf(options1, options_template_5, metric_type, 1, 1, 1, 60, 48, "REF", "no");
    sprintf(options2, options_template_5, metric_type, 1, 1, 2 * i - 1, 60, 48, "GPU",
            "yes");
    sprintf(options3, options_template_5, metric_type, 1, 1, 2 * i, 60, 48, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
    sprintf(options1, options_template_5, metric_type, 1, 1, 1, 60, 48, "REF", "no");
    sprintf(options2, options_template_5, metric_type, 1, 1, i, 60, 48, "GPU", "yes");
    sprintf(options3, options_template_5, metric_type, 1, 2, i, 60, 24, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
  }

  //----------
  //---num_repl, 2-way
  //----------

  char options_template_10[] =
      "--sparse yes "
      "--metric_type %s "
      //"--problem_type analytic "
      "--num_field 4 --num_vector_local %i --compute_method %s --all2all yes "
      "--num_proc_vector %i --num_proc_repl %i "
      "--num_proc_field %i --num_way %i --num_stage %i";

  int num_vector_local = 0;
  int num_proc_vector = 0;
  int num_proc_repl = 0;
  int gpu = 0;
  int num_stage = 1;

  for (gpu=0; gpu<=1; ++gpu) {
    for (num_vector_local=4; num_vector_local<=5; ++num_vector_local) {
      for (num_proc_vector=1; num_proc_vector<=4; ++num_proc_vector) {
        for (num_proc_repl=2; num_proc_repl<=5; ++num_proc_repl) {
          const int num_proc_field = gpu ? 2 : 1;
          const int num_way = 2;
          sprintf(options1, options_template_10, metric_type,
                  num_vector_local*num_proc_vector,
                  "GPU", 1, 1,
                  1, num_way, 1);
          sprintf(options2, options_template_10, metric_type, num_vector_local,
                  gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                  num_proc_field, num_way, num_stage);
          test_2runs(options1, options2);
        }
      }
    }
  }

  //----------
  //---num_repl, num_stage, 3-way
  //----------

  if (!is_duo) {
    for (gpu=0; gpu<=1; ++gpu) {
      for (num_vector_local=6; num_vector_local<=18; num_vector_local+=12) {
        for (num_proc_vector=1; num_proc_vector<=4; ++num_proc_vector) {
          for (num_proc_repl=2; num_proc_repl<=5; ++num_proc_repl) {
            for (num_stage=1; num_stage<=6; num_stage+=4) {
              const int num_proc_field = gpu ? 2 : 1;
              const int num_way = 3;
              sprintf(options1, options_template_10, metric_type,
                      num_vector_local*num_proc_vector,
                      "GPU", 1, 1,
                      1, num_way, 1);
              sprintf(options2, options_template_10, metric_type,
                      num_vector_local,
                      gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                      num_proc_field, num_way, num_stage);
              test_2runs(options1, options2);
            }
          }
        }
      }
    }
  }

  //----------
  //---ccc_param
  //----------

  char options_template_11a[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      //"--problem_type analytic "
      "--num_proc_vector 1 --num_field 30 --num_vector_local 40 "
      "--compute_method GPU";
  char options_template_11b[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      //"--problem_type analytic "
      "--num_proc_vector 1 --num_field 30 --num_vector_local 40 "
      "--compute_method GPU --ccc_param %.20e --problem_type analytic";

  sprintf(options1, options_template_11a, metric_type, 1);
  sprintf(options2, options_template_11b, metric_type, 1, ((double)2) / ((double)3));
  EXPECT_EQ(true, compare_2runs(options1, options2));

  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

  sprintf(options1, options_template_11a, metric_type, 1);
  sprintf(options2, options_template_11b, metric_type, 1, ((double)1) / ((double)2));
  const int result11 = compare_2runs(options1, options2);
  EXPECT_EQ(true, proc_num==0 ? ! result11 : true);

  //----------
  //---num_phase, 2-way
  //----------

  int num_phase = 0;

  for (num_proc_vector=1; num_proc_vector<=8; ++num_proc_vector) {
    for (num_proc_repl=1; num_proc_repl<=8; ++num_proc_repl) {
      for (num_phase=2; num_phase<=8; ++num_phase) {
        if (!(num_phase <= 1 + num_proc_vector/2)) {
          continue;
        }
        char options_template[] =
          "--sparse yes "
          "--metric_type %s "
          //"--problem_type analytic "
          "--num_field 7 --num_vector 12 --compute_method GPU --all2all yes "
          "--num_proc_vector %i --num_proc_repl %i --num_phase %i "
          "--num_way 2";
        sprintf(options1, options_template, metric_type, num_proc_vector, num_proc_repl,
                num_phase);
        sprintf(options2, options_template, metric_type, 1, 1, 1);
        test_2runs(options1, options2);
      }
    }
  }

  //----------
  //---2-way, all2all yes, large, sparse
  //----------

  {
    char options_template_2[] =
      "--sparse yes "
      "--metric_type %s --verbosity %i "
      //"--problem_type analytic "
      "--num_proc_vector %i --num_field %i --num_vector %i "
      "--compute_method %s --all2all %s";

    const int nf = 100;
    const int nv = 48;

    sprintf(options1, options_template_2, metric_type, 1, 1, nf, nv, "REF", "no");
    sprintf(options2, options_template_2, metric_type, 1, 1, nf, nv, "CPU", "yes");
    sprintf(options3, options_template_2, metric_type, 1, 1, nf, nv, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    sprintf(options1, options_template_2, metric_type, 1, 1, nf, nv, "REF", "no");
    sprintf(options2, options_template_2, metric_type, 1, 2, nf, nv, "CPU", "yes");
    sprintf(options3, options_template_2, metric_type, 1, 2, nf, nv, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
  }

  //----------
  //---3-way, all2all yes, large, sparse
  //----------

  if (!is_duo) {
    char options_template_4[] =
        "--sparse yes "
        "--metric_type %s --verbosity %i "
        //"--problem_type analytic "
        "--num_proc_vector %i --num_field %i --num_vector_local %i "
        "--compute_method %s --all2all %s --num_way 3";

    sprintf(options1, options_template_4, metric_type, 1, 1, 100, 48, "REF", "no");
    sprintf(options2, options_template_4, metric_type, 1, 1, 100, 48, "CPU", "yes");
    sprintf(options3, options_template_4, metric_type, 1, 1, 100, 48, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));

    sprintf(options1, options_template_4, metric_type, 1, 1, 100, 48, "REF", "no");
    sprintf(options2, options_template_4, metric_type, 1, 4, 100, 12, "CPU", "yes");
    sprintf(options3, options_template_4, metric_type, 1, 4, 100, 12, "GPU", "yes");
    EXPECT_EQ(true, compare_3runs(options1, options2, options3));
  }

  //----------
  //---file input
  //----------

  for (int num_vector = 13; num_vector <= 13; ++num_vector) {
    for (int num_field = 1; num_field <= 3*300; ++num_field) {
      create_vectors_file("ccc_duo_2way_in.bin", num_field, num_vector,
                          GM_METRIC_TYPE_CCC, 2, problem_type_default(), 1);
      for (int num_proc_vector=3; num_proc_vector<=3; ++num_proc_vector) {
        for (int num_proc_field=1; num_proc_field<=3; ++num_proc_field) {

          const bool is_nearly_multiple_of_32 =
            (num_field/num_proc_field+4) % 32 <= 8;
          if (! is_nearly_multiple_of_32) {
            continue;
          }

          char options1[1024];
          char options2[1024];

          char options_template[] =
                 "--sparse yes "
                 "--metric_type %s "
                 "--num_vector %i --num_field %i "
                 "--num_proc_vector %i --num_proc_field %i "
                 "--compute_method GPU --all2all yes %s --verbosity 1";

          sprintf(options1, options_template, metric_type,
                  num_vector, num_field, num_proc_vector, num_proc_field, "");

          sprintf(options2, options_template, metric_type,
                  num_vector, num_field, num_proc_vector, num_proc_field,
                  "--input_file ccc_duo_2way_in.bin");

          test_2runs(options1, options2);
        }
      }
    }
  }

  //----------
  //---file output, 2-way
  //----------

  {
    char options_template[] =
        "--metric_type %s "
        "--num_proc_vector 1 --num_proc_vector 1 "
        "--num_field 7 --num_vector 100 "
        "--num_way 2 "
        "--all2all yes --sparse %s "
        "--problem_type random %s";
    for (int sparse=0; sparse<2; ++sparse) {
      if (is_duo && sparse == 0) continue;
      sprintf(options1, options_template, metric_type, sparse==1 ? "yes" : "no",
              "--compute_method REF");
      sprintf(options2, options_template, metric_type, sparse==1 ? "yes" : "no",
              "--compute_method GPU "
              "--verbosity 1 "
              "--threshold .5 "
              "--output_file_stub test_ccc_duo_2way");
      test_2runs(options1, options2);
    }
  }

  //----------
  //---file output, 3-way
  //----------

  if (!is_duo) {
    char options_template[] =
        "--metric_type %s "
        "--num_proc_vector 1 --num_proc_vector 1 "
        "--num_field 7 --num_vector 18 "
        "--num_way 3 "
        "--all2all yes --sparse %s "
        "--problem_type random %s";
    for (int sparse=0; sparse<2; ++sparse) {
      if (is_duo && sparse == 0) continue;
      sprintf(options1, options_template, metric_type, sparse==1 ? "yes" : "no",
              "--compute_method REF");
      sprintf(options2, options_template, metric_type, sparse==1 ? "yes" : "no",
              "--compute_method GPU "
              "--verbosity 1 "
              "--threshold .1 "
              "--output_file_stub test_ccc_duo_3way");
      test_2runs(options1, options2);
    }
  }
#endif
} // DriverTest_ccc__duo_

//=============================================================================

void DriverTest_ccc_() {
  DriverTest_ccc_duo_("ccc");
}

//=============================================================================

void DriverTest_duo_() {
  DriverTest_ccc_duo_("duo");
}

//=============================================================================

TEST(DriverTest, duo2_simple_sparse) {
  DriverTest_duo2_simple_sparse_();
}

TEST(DriverTest, duo) {
  DriverTest_duo_();
}

TEST(DriverTest, ccc2_simple) {
  DriverTest_ccc2_simple_();
}

TEST(DriverTest, ccc2_simple_sparse) {
  DriverTest_ccc2_simple_sparse_();
}

TEST(DriverTest, ccc3_simple) {
  DriverTest_ccc3_simple_();
}

TEST(DriverTest, ccc3_simple_sparse) {
  DriverTest_ccc3_simple_sparse_();
}

TEST(DriverTest, ccc) {
  DriverTest_ccc_();
}

TEST(DriverTest, czek) {
  DriverTest_czek_();
}

//=============================================================================

GTEST_API_ int main(int argc, char** argv) {

  ::testing::InitGoogleTest(&argc, argv);

  MPI_Init(&argc, &argv);

  int comm_rank = 0;
  int mpi_code = 0;
  mpi_code = MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  GMInsist(mpi_code == MPI_SUCCESS);

  if (comm_rank != 0) {
    ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release(listeners.default_result_printer());
  }

  int result = RUN_ALL_TESTS();

  int result_g = 11;

  mpi_code = MPI_Allreduce(&result, &result_g,
                           1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  GMInsist(mpi_code == MPI_SUCCESS);

  MPI_Finalize();
  return result_g;
}

//=============================================================================

//-----------------------------------------------------------------------------
