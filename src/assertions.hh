//-----------------------------------------------------------------------------
/*!
 * \file   assertions.hh
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

#ifndef _gm_assertions_hh_
#define _gm_assertions_hh_

//=============================================================================
// Macros.

//-----------------------------------------------------------------------------
/// \brief Insist macro - assert a condition even for release builds.
///
///        This should be used only for non-performance-sensitive
///        code locations -- e.g., not in a deep loop nest.

#define GMInsist(condition) \
  (void)((condition) || (CoMet::assert_(#condition, __FILE__, __LINE__), 0))

//-----------------------------------------------------------------------------
/// \brief Insist macro specifically for a user-caused error condition.

#define GMInsistInterface(env, condition) \
  (void)((condition) || \
         (CoMet::insist_interface(env, #condition, __FILE__, __LINE__), 0))

//-----------------------------------------------------------------------------
/// \brief Assertion macro (for debug builds only).

#ifndef NDEBUG
#define GM_ASSERTIONS_ON
#define GMAssert(condition) GMInsist(condition)
#else
#define GMAssert(condition)
#endif

//-----------------------------------------------------------------------------
/// \brief Static (i.e., compile time) assertion macro.

// TODO: replace with C++11 equivalent.

#ifndef NDEBUG
// Fail compilation and (hopefully) give a filename/line number
#define GMStaticAssert(condition) \
  {                               \
    int a[(condition) ? 1 : -1];  \
    (void) a;                     \
  }
#else
#define GMStaticAssert(condition)
#endif

//=============================================================================
// Declarations.

namespace CoMet {

//-----------------------------------------------------------------------------
/// \brief Function to support the GMAssert, GMInsist macros.

//        Use trailing underscore to avoid collision with C assert macro.

void assert_(const char* condition_string, const char* file, int line);

//-----------------------------------------------------------------------------
/// \brief Function to support the GMInsistInterface macro.

void insist_interface(void const * const env,
                         const char* condition_string,
                         const char* file,
                         int line);

//=============================================================================

} // namespace CoMet

//-----------------------------------------------------------------------------

#endif // _gm_assertions_hh_

//-----------------------------------------------------------------------------
