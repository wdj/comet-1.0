//-----------------------------------------------------------------------------
/*!
 * \file   assertions.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Macros and code for assertions.
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
#include "errno.h"

#ifdef TESTING
#include "gtest/gtest.h"
#endif

#include "env.hh"
#include "assertions.hh"

//=============================================================================

namespace CoMet {

//-----------------------------------------------------------------------------
/// \brief For a test build, cause test harness to throw an error.

static void make_test_harness_failure() {
#ifdef TESTING
  ASSERT_TRUE(0);
#endif
}

//-----------------------------------------------------------------------------
/// \brief Function to support the GMAssert, GMInsist macros.

void assert_(const char* condition_string, const char* file, int line) {
  fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Assertion error",
          condition_string, file, line);
  make_test_harness_failure();
#ifdef GM_ASSERTIONS_ON
  //raise(SIGUSR1);
  exit(EXIT_FAILURE);
#else
  exit(EXIT_FAILURE);
#endif
}

//-----------------------------------------------------------------------------
/// \brief Function to support the GMInsistInterface macro.

void insist_interface(const void* const env,
                         const char* condition_string,
                         const char* file,
                         int line) {
  if (GMEnv_proc_num((const GMEnv* const)env) == 0) {
    fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Interface error",
            condition_string, file, line);
  }
  make_test_harness_failure();
  exit(EXIT_FAILURE);
}

//=============================================================================

} // namespace CoMet

//-----------------------------------------------------------------------------
