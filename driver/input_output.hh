//-----------------------------------------------------------------------------
/*!
 * \file   input_output.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O functions used by driver, header.
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

#ifndef _gm_input_output_hh_
#define _gm_input_output_hh_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"

#include "driver.hh"

//=============================================================================

void set_vectors_from_file(GMVectors* vectors, DriverOptions* do_, GMEnv* env);

void write_vectors_to_file(GMVectors* vectors, const char* vectors_file_path,
                           GMEnv* env);

//=============================================================================
// Class to help output the result metrics values to file

FILE* gm_metrics_file_open(char* metrics_file_path_stub, GMEnv* env);

class MetricWriter {
public:

  MetricWriter(FILE* file, GMMetrics* metrics, GMEnv* env);

  ~MetricWriter() {}

  size_t get_num_written() {return this->num_written_total_;}

  void write(size_t coord0, size_t coord1, GMFloat value);

  void write(size_t coord0, size_t coord1, size_t coord2, GMFloat value);

  void write(size_t coord0, size_t coord1, int i0, int i1, GMFloat value);

  void write(size_t coord0, size_t coord1, size_t coord2,
             int i0, int i1, int i2, GMFloat value);

private:

  FILE* file_;
  const int data_type_;
  int num_way_;
  GMEnv* env_;

  size_t num_written_total_;

  //---Disallowed methods.

  MetricWriter(  const MetricWriter&);
  void operator=(const MetricWriter&);

};

//=============================================================================

class MetricsFile {
public:

  MetricsFile(DriverOptions* do_, GMEnv* env);

  ~MetricsFile();

  void write(GMMetrics* metrics, GMEnv* env);

  size_t get_num_written() {return num_written_;}

private:

  FILE* file_;
  int verbosity_;
  double threshold_;
  size_t num_written_;

  //---Disallowed methods.

  MetricsFile(   const MetricsFile&);
  void operator=(const MetricsFile&);
};

//=============================================================================

#endif // _gm_input_output_hh_

//-----------------------------------------------------------------------------
