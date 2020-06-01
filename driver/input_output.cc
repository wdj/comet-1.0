//-----------------------------------------------------------------------------
/*!
 * \file   input_output.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O functions used by driver.
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

#include "cstdio"
#include "cstdint"
#include "string.h"
//#include "stdlib.h"
//#include "stddef.h"
//#include "float.h"
//#include "errno.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"

#include "driver.hh"
#include "input_output.hh"

//=============================================================================
// Input vectors from files

void set_vectors_from_file_float(GMVectors* vectors, DriverOptions* do_,
                                 GMEnv* env) {
  GMInsist(vectors && do_ && env && do_->input_file_path);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  FILE* input_file = fopen(do_->input_file_path, "r");
  GMInsist(NULL != input_file && "Unable to open file.");

  typedef GMFloat inval_t;

  const size_t proc_num_v = GMEnv_proc_num_vector_i(env);
  const size_t proc_num_f = GMEnv_proc_num_field(env);

  const size_t nva = vectors->dm->num_vector_active;
  const size_t nfa = vectors->dm->num_field_active;

  const size_t bytes_per_vector_file = nfa * sizeof(inval_t);
  const int fl_min = 0;
  const size_t f_min = fl_min + vectors->num_field_local * proc_num_f;

  const size_t invals_to_read = vectors->dm->num_field_active_local;

  // Loop to input vectors

  for (int vl = 0; vl < vectors->num_vector_local; ++vl) {

    const size_t v = vl + vectors->num_vector_local * proc_num_v;
    // Fill the pad vectors with copies of the last vector
    const size_t v_file = gm_min_i8(v, nva-1);
    // Offset into file to first byte to read
    const size_t addr_min_file = f_min * sizeof(inval_t) +
      bytes_per_vector_file * v_file;

    int fseek_success = fseek(input_file, addr_min_file, SEEK_SET);
    GMInsist(0 == fseek_success && "File seek failure.");

    // First location in memory to store into
    inval_t* const addr_min_mem = GMVectors_float_ptr(vectors, fl_min, vl, env);
    // NOTE: the following call is ok since has no side effects
    GMInsist((fl_min+1 >= vectors->num_field_local ||
        GMVectors_float_ptr(vectors, fl_min+1, vl, env) == addr_min_mem + 1)
        && "Vector layout is incompatible with operation.");

      const size_t num_read = fread(addr_min_mem, sizeof(inval_t),
                                    invals_to_read, input_file);
      GMInsist(invals_to_read == num_read && "File read failure.");

  } //---for vl

  // Ensure any end of vector pad set to zero
  GMVectors_initialize_pad(vectors, env);

  fclose(input_file);
}

//-----------------------------------------------------------------------------

void set_vectors_from_file_bits2(GMVectors* vectors, DriverOptions* do_,
                                 GMEnv* env) {
  GMInsist(vectors && do_ && env && do_->input_file_path);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  FILE* input_file = fopen(do_->input_file_path, "r");
  GMInsist(NULL != input_file && "Unable to open file.");

  typedef unsigned char inval_t;

  const size_t proc_num_v = GMEnv_proc_num_vector_i(env);

  const int bits_per_field = vectors->num_bits_per_val; // = 2
  const int bits_per_byte = 8;
  const int fields_per_byte = bits_per_byte / bits_per_field;
  const size_t bytes_per_packedval = vectors->num_bits_per_packedval /
               bits_per_byte;
  // The storage format assumes each vector padded to a whole number of bytes.

  const size_t bytes_per_vector_file
    = gm_ceil_i8(vectors->num_field_active, fields_per_byte);

  const int npvfl = vectors->num_packedval_field_local;
  const size_t nfal = vectors->dm->num_field_active_local;

  const size_t fl_min = 0;
  const size_t f_min = fl_min + vectors->dm->field_base;
  const size_t f_max = f_min + nfal;

  const size_t nva = vectors->dm->num_vector_active;

  const size_t buf_size = nfal * fields_per_byte + 1;
  inval_t* buf = (inval_t*)malloc(buf_size * sizeof(*buf));

  // Loop to input vectors

  for (int vl = 0; vl < vectors->num_vector_local; ++vl) {

    const size_t v = vl + vectors->num_vector_local * proc_num_v;

    // Fill the pad vectors with copies of the last vector
    const size_t v_file = gm_min_i8(v, nva-1);

    // Byte in file where this vector starts
    const size_t v_base = bytes_per_vector_file * v_file;

    // Offset into file to first byte to read
    const size_t addr_min_file = gm_floor_i8(f_min, fields_per_byte) + v_base;

    // Offset into file to (1 plus) last byte to read
    const size_t addr_max_file = gm_ceil_i8(f_max, fields_per_byte) + v_base;

    // total num bytes to read
    const size_t invals_to_read = addr_max_file - addr_min_file;

    const int f_offset = f_min % fields_per_byte;

    // the num bytes read could be at most one additional byte, if straddle
    GMInsist(invals_to_read <= npvfl*bytes_per_packedval + (f_offset ? 1 : 0));

    int fseek_success = fseek(input_file, addr_min_file, SEEK_SET);
    GMInsist(0 == fseek_success && "File seek failure.");

    const size_t num_read = fread(buf, sizeof(inval_t),
                                  invals_to_read, input_file);
    GMInsist(invals_to_read == num_read && "File read failure.");

    // Initialize elements buffer to store into vectors struct
    GMBits2x64 outval = GMBits2x64_null();

    const size_t invals_to_read_1 = gm_min_i8(invals_to_read,
      npvfl*bytes_per_packedval);

    for (size_t i = 0; i < invals_to_read_1; ++i) {
      inval_t* const p = buf + i;
      const int vlo = p[0];
      const int vhi = (i+1 < invals_to_read) ? p[1] : 0;
      const inval_t inval = ( (vhi << (8-2*f_offset)) |
                              (vlo >> (2*f_offset)) ) & (int)255;
      const int wordnum = (i % 16) / 8;
      outval.data[wordnum] += ((uint64_t)inval) << ((i % 8) * 8) ;
      if (i % 16 == 15 || i == invals_to_read_1-1) {
        // Flush buffer
        const size_t pfl = i / 16;
        GMVectors_bits2x64_set(vectors, pfl, vl, outval, env);
        outval = GMBits2x64_null();
      }

    } // for p

  } /*---vl---*/

  // Ensure any end of vector pad set to zero
  GMVectors_initialize_pad(vectors, env);

  free(buf);
  fclose(input_file);
}

//-----------------------------------------------------------------------------

void set_vectors_from_file(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMInsist(vectors && do_ && env && do_->input_file_path);

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/

      set_vectors_from_file_float(vectors, do_, env);

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/

      set_vectors_from_file_bits2(vectors, do_, env);

    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
}

//=============================================================================

void write_vectors_to_file(GMVectors* vectors, const char* vectors_file_path, 
                           GMEnv* env) {
  GMInsist(vectors && vectors_file_path && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  GMInsistInterface(env, GMEnv_num_proc(env) == 1 &&
                    "Only single proc case supported.");

  FILE* vectors_file = fopen(vectors_file_path, "w");
  GMInsist(NULL != vectors_file && "Unable to open file.");

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/

      typedef GMFloat outval_t;

      const size_t nval = vectors->dm->num_vector_active_local;
      const size_t nfal = vectors->dm->num_field_active_local;

      size_t num_written_total = 0;
      const size_t num_written_attempted_total = nval * nfal;

      for (size_t vl = 0 ; vl < nval; ++vl) {
        for (size_t fl = 0 ; fl < nfal; ++fl) {

          const outval_t outv = GMVectors_float_get(vectors, fl, vl, env);
          const size_t num_written = fwrite(&outv, sizeof(outval_t), 1,
                                            vectors_file);
          num_written_total += num_written;
        }
      }

      GMInsist(num_written_attempted_total == num_written_total &&
               "File write failure.");

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/

      typedef unsigned char outval_t;

      const int bits_per_field = vectors->num_bits_per_val; // = 2
      const int bits_per_byte = 8;
      const int fields_per_byte = bits_per_byte / bits_per_field;
      const size_t bytes_per_packedval = vectors->num_bits_per_packedval /
                   bits_per_byte;

      const size_t nval = vectors->dm->num_vector_active_local;
      const size_t npfl = vectors->dm->num_packedfield_local;

      size_t num_written_total = 0;
      size_t num_written_attempted_total = 0;

      for (size_t vl = 0 ; vl < nval; ++vl) {
        for (size_t pfl = 0 ; pfl < npfl; ++pfl) {

          GMBits2x64 val = GMVectors_bits2x64_get(vectors, pfl, vl, env);

          // Calculate begin and end byte numbers of vector to write

          const size_t offset_min = pfl * bytes_per_packedval;

          const size_t byte_max =
            gm_ceil_i8(vectors->dm->num_field_active_local, fields_per_byte);

          const size_t offset_max = gm_min_i8((pfl+1) * bytes_per_packedval,
                                              byte_max);

          // Loop over bytes to output

          for (size_t offset = offset_min, i = 0; offset < offset_max;
               ++offset, ++i) {

            const int index0 = i % 8;
            const int index1 = i / 8;

            const outval_t outval = ((val.data[index1] >> (8*index0)) &
                                    (uint64_t)255);

            const size_t num_to_write = 1;

            const size_t num_written = fwrite(&outval, sizeof(outval_t),
                                              num_to_write, vectors_file);

            num_written_total += num_written;
            num_written_attempted_total += num_to_write;
          }

        } // pfl
      } // vl

      GMInsist(num_written_attempted_total == num_written_total &&
               "File write failure.");

    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/

  fclose(vectors_file);
}

//=============================================================================
// Output results metrics to file: implementation

void output_metrics_tally2x2_bin_(GMMetrics* metrics, FILE* file,
                     double threshold, size_t& num_written, GMEnv* env) {
  GMInsist(metrics && file && env);
  GMInsist(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMInsist(file != stdout);

  MetricWriter writer(file, metrics, env);

  //enum {vals_per_index = 4;}

  // Number of index values values visited for one pass across buffer
  const size_t num_buf_ind = 1000 * 1000;
  // Number of (index, i0, i1) entries (potentially) stored in buffer
  const size_t num_buf = 4 * num_buf_ind;

  // Each buffer entry contains: whether value is to be written,
  // coord0, coord1, i0, i1, and value

  char do_out_buf[num_buf];
  int coord0_buf[num_buf];
  int coord1_buf[num_buf];
  int i01_buf[num_buf];
  GMFloat value_buf[num_buf];

  for (int i=0; i<(int)num_buf; ++i) {
    do_out_buf[i] = 0;
  }

  const GMFloat threshold_eff = threshold<0. ? -1e20 : threshold;

  // Process num_buf_ind index values at a time
  for (size_t ind_base = 0; ind_base < metrics->num_elts_local;
       ind_base += num_buf_ind) {
    const size_t ind_max = gm_min_i8(metrics->num_elts_local,
                                     ind_base + num_buf_ind);

    if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC) {

      // Fill buffer
#pragma omp parallel for schedule(dynamic,1000)
      for (size_t index = ind_base; index < ind_max; ++index) {
        // Do any of the values exceed the threshold
        if (GMMetrics_ccc_get_from_index_2_threshold(
               metrics, index, threshold_eff, env)) {
          for (int i0 = 0; i0 < 2; ++i0) {
            for (int i1 = 0; i1 < 2; ++i1) {
              const GMFloat value =
                GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);
              const bool pass_thresh = value>threshold_eff;
              if (pass_thresh) {
                const size_t coord0 =
                  GMMetrics_coord_global_from_index(metrics, index, 0, env);
                const size_t coord1 =
                  GMMetrics_coord_global_from_index(metrics, index, 1, env);
                const char do_out = coord0 < metrics->num_vector_active &&
                                    coord1 < metrics->num_vector_active;
                const size_t ind_buf = i1 + 2*(i0 + 2*(index-ind_base));
                do_out_buf[ind_buf] = do_out;
                coord0_buf[ind_buf] = coord0;
                coord1_buf[ind_buf] = coord1;
                i01_buf[ind_buf] = i0 + 2*i1;
                value_buf[ind_buf] = value;
              } /*---if---*/
            } /*---i1---*/
          } /*---i0---*/
        }
      } /*---for index---*/

      // Flush buffer

      // Number of buffer entries to visit
      const size_t ind_buf_max = 4 * (ind_max - ind_base);
      // Use 4 byte integer ptr to check 4 chars at a time, for speed
      typedef int multi_t;
      //assert(sizeof(multi_t) == 4 * sizeof(char));
      const multi_t* const do_out_ptr_max = (multi_t*)(do_out_buf + ind_buf_max);
      multi_t* do_out_ptr = (multi_t*)do_out_buf;
      for (; do_out_ptr < do_out_ptr_max;) {
        // Do any of the 4 need to be output
        if (*(do_out_ptr++)) {
          for (int i=0; i<4; ++i) {
            const size_t ind_buf = (do_out_ptr - (multi_t*)do_out_buf - 1)*4 + i;
            if (do_out_buf[ind_buf]) {
              const int i0 = i01_buf[ind_buf] % 2;
              const int i1 = i01_buf[ind_buf] / 2;
              const size_t coord0 = coord0_buf[ind_buf];
              const size_t coord1 = coord1_buf[ind_buf];
              writer.write(coord0, coord1, i0, i1, value_buf[ind_buf]);
              // Reset buffer entry to false
              do_out_buf[ind_buf] = 0;
            }
          }
        }
      } /*---ind_buf---*/

    } else { // DUO

      // Fill buffer
#pragma omp parallel for schedule(dynamic,1000)
      for (size_t index = ind_base; index < ind_max; ++index) {
        // Do any of the values exceed the threshold
        if (GMMetrics_duo_get_from_index_2_threshold(
               metrics, index, threshold_eff, env)) {
          for (int i0 = 0; i0 < 2; ++i0) {
            for (int i1 = 0; i1 < 2; ++i1) {
              const GMFloat value =
                GMMetrics_duo_get_from_index_2(metrics, index, i0, i1, env);
              const bool pass_thresh = value>threshold_eff;
              if (pass_thresh) {
                const size_t coord0 =
                  GMMetrics_coord_global_from_index(metrics, index, 0, env);
                const size_t coord1 =
                  GMMetrics_coord_global_from_index(metrics, index, 1, env);
                const char do_out = coord0 < metrics->num_vector_active &&
                                    coord1 < metrics->num_vector_active;
                const size_t ind_buf = i1 + 2*(i0 + 2*(index-ind_base));
                do_out_buf[ind_buf] = do_out;
                coord0_buf[ind_buf] = coord0;
                coord1_buf[ind_buf] = coord1;
                i01_buf[ind_buf] = i0 + 2*i1;
                value_buf[ind_buf] = value;
              } /*---if---*/
            } /*---i1---*/
          } /*---i0---*/
        }
      } /*---for index---*/

      // Flush buffer

      // Number of buffer entries to visit
      const size_t ind_buf_max = 4 * (ind_max - ind_base);
      // Use 4 byte integer ptr to check 4 chars at a time, for speed
      typedef int multi_t;
      //assert(sizeof(multi_t) == 4 * sizeof(char));
      const multi_t* const do_out_ptr_max = (multi_t*)(do_out_buf + ind_buf_max);
      multi_t* do_out_ptr = (multi_t*)do_out_buf;
      for (; do_out_ptr < do_out_ptr_max;) {
        // Do any of the 4 need to be output
        if (*(do_out_ptr++)) {
          for (int i=0; i<4; ++i) {
            const size_t ind_buf = (do_out_ptr - (multi_t*)do_out_buf - 1)*4 + i;
            if (do_out_buf[ind_buf]) {
              const int i0 = i01_buf[ind_buf] % 2;
              const int i1 = i01_buf[ind_buf] / 2;
              const size_t coord0 = coord0_buf[ind_buf];
              const size_t coord1 = coord1_buf[ind_buf];
              writer.write(coord0, coord1, i0, i1, value_buf[ind_buf]);
              // Reset buffer entry to false
              do_out_buf[ind_buf] = 0;
            }
          }
        }
      } /*---ind_buf---*/

    } // metric_type

  } /*---ind_base---*/

  num_written += writer.get_num_written();
}

//=============================================================================

void output_metrics_tally4x2_bin_(GMMetrics* metrics, FILE* file,
                     double threshold, size_t& num_written, GMEnv* env) {
  GMInsist(metrics && file && env);
  GMInsist(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMInsist(file != stdout);
  GMInsist(GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC);

  MetricWriter writer(file, metrics, env);

  //enum {vals_per_index = 8;}

  // Number of index values values visited for one pass across buffer
  const size_t num_buf_ind = 1000 * 1000;
  // Number of (index, i0, i1) entries (potentially) stored in buffer
  const size_t num_buf = 8 * num_buf_ind;

  // Each buffer entry contains: whether value is to be written,
  // coord0, coord1, coord2, i0, i1, i2, and value

  char do_out_buf[num_buf];
  int coord0_buf[num_buf];
  int coord1_buf[num_buf];
  int coord2_buf[num_buf];
  int i012_buf[num_buf];
  GMFloat value_buf[num_buf];

  for (int i=0; i<(int)num_buf; ++i) {
    do_out_buf[i] = 0;
  }

  const GMFloat threshold_eff = threshold<0. ? -1e20 : threshold;

  // Process num_buf_ind index values at a time
  for (size_t ind_base = 0; ind_base < metrics->num_elts_local;
       ind_base += num_buf_ind) {
    const size_t ind_max = gm_min_i8(metrics->num_elts_local,
                                     ind_base + num_buf_ind);
    // Fill buffer
#pragma omp parallel for schedule(dynamic,1000)
    for (size_t index = ind_base; index < ind_max; ++index) {
      // Do any of the values exceed the threshold
      if (GMMetrics_ccc_get_from_index_3_threshold(
             metrics, index, threshold_eff, env)) {
        for (int i0 = 0; i0 < 2; ++i0) {
          for (int i1 = 0; i1 < 2; ++i1) {
            for (int i2 = 0; i2 < 2; ++i2) {
              const GMFloat value =
                GMMetrics_ccc_get_from_index_3(metrics, index, i0, i1, i2, env);
              const bool pass_thresh = value>threshold_eff;
              if (pass_thresh) {
                const size_t coord0 =
                  GMMetrics_coord_global_from_index(metrics, index, 0, env);
                const size_t coord1 =
                  GMMetrics_coord_global_from_index(metrics, index, 1, env);
                const size_t coord2 =
                  GMMetrics_coord_global_from_index(metrics, index, 2, env);
                const char do_out = coord0 < metrics->num_vector_active &&
                                    coord1 < metrics->num_vector_active &&
                                    coord2 < metrics->num_vector_active;
                const size_t ind_buf = i1 + 2*(i0 + 2*(i2 +2*(index-ind_base)));
                do_out_buf[ind_buf] = do_out;
                coord0_buf[ind_buf] = coord0;
                coord1_buf[ind_buf] = coord1;
                coord2_buf[ind_buf] = coord2;
                i012_buf[ind_buf] = i0 + 2*(i1 + 2*i2);
                value_buf[ind_buf] = value;
              } /*---if---*/
            } /*---i2---*/
          } /*---i1---*/
        } /*---i0---*/
      }
    } /*---for index---*/

    // Flush buffer

    // Number of buffer entries to visit
    const size_t ind_buf_max = 8 * (ind_max - ind_base);
    // Use 4 byte integer ptr to check 4 chars at a time, for speed
    typedef size_t multi_t;
    //assert(sizeof(multi_t) == 8 * sizeof(char));
    const multi_t* const do_out_ptr_max = (multi_t*)(do_out_buf + ind_buf_max);
    multi_t* do_out_ptr = (multi_t*)do_out_buf;
    for (; do_out_ptr < do_out_ptr_max;) {
      // Do any of the 8 need to be output
      if (*(do_out_ptr++)) {
        for (int i=0; i<8; ++i) {
          const size_t ind_buf = (do_out_ptr - (multi_t*)do_out_buf - 1)*8 + i;
          if (do_out_buf[ind_buf]) {
            const int i0 = i012_buf[ind_buf] % 2;
            const int i1 = (i012_buf[ind_buf] / 2) % 2;
            const int i2 = i012_buf[ind_buf] / 4;
            const size_t coord0 = coord0_buf[ind_buf];
            const size_t coord1 = coord1_buf[ind_buf];
            const size_t coord2 = coord2_buf[ind_buf];
            writer.write(coord0, coord1, coord2, i0,i1,i2, value_buf[ind_buf]);
            // Reset buffer entry to false
            do_out_buf[ind_buf] = 0;
          }
        }
      }
    } /*---ind_buf---*/
  } /*---ind_base---*/

  num_written += writer.get_num_written();
}

//=============================================================================

void output_metrics_(GMMetrics* metrics, FILE* file,
                     double threshold, size_t& num_written, GMEnv* env) {
  GMInsist(metrics && file && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  /*---Due to redundancy, only results from some processors are needed---*/
  if (GMEnv_proc_num_field(env) != 0) {
    return;
  }

  switch (GMEnv_data_type_metrics(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/

      /*----------*/
      if (GMEnv_num_way(env) == GM_NUM_WAY_2) {
      /*----------*/

        MetricWriter writer(file, metrics, env);

        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t coord0 =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t coord1 =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          if (coord0 >= metrics->num_vector_active ||
              coord1 >= metrics->num_vector_active) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czek_get_from_index(metrics, index, env);
          if (!(threshold < 0. || value > threshold)) {
            continue;
          }
          /*---Output the value---*/
          if (file == stdout) {

            fprintf(file,
              sizeof(GMFloat) == 8 ?
              "element (%li,%li): value: %.17e\n" :
              "element (%li,%li): value: %.8e\n", coord0, coord1, value);

          } else {

            writer.write(coord0, coord1, value);

          }
        } /*---for index---*/

        num_written += writer.get_num_written();

      } /*---if---*/

      /*----------*/
      if (GMEnv_num_way(env) == GM_NUM_WAY_3) {
      /*----------*/

        MetricWriter writer(file, metrics, env);

        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t coord0 =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t coord1 =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          const size_t coord2 =
            GMMetrics_coord_global_from_index(metrics, index, 2, env);
          if (coord0 >= metrics->num_vector_active ||
              coord1 >= metrics->num_vector_active ||
              coord2 >= metrics->num_vector_active) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czek_get_from_index(metrics, index, env);
          if (!(threshold < 0. || value > threshold)) {
            continue;
          }
          /*---Output the value---*/
          if (file == stdout) {

            fprintf(file,
              sizeof(GMFloat) == 8 ?
              "element (%li,%li,%li): value: %.17e\n" :
              "element (%li,%li,%li): value: %.8e\n",
              coord0, coord1, coord2, value);

          } else {

            writer.write(coord0, coord1, coord2, value);

          }
        } /*---for index---*/

        num_written += writer.get_num_written();

      } /*---if---*/

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
    /*--------------------*/
      GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_2);

      if (file != stdout) {

        output_metrics_tally2x2_bin_(metrics, file, threshold,
                                     num_written, env);

      } else /*---stdout---*/ {
        MetricWriter writer(file, metrics, env);

        size_t index = 0;
        for (index = 0; index < metrics->num_elts_local; ++index) {
          const size_t coord0 =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          if (coord0 >= metrics->num_vector_active) {
            continue;
          }
          const size_t coord1 =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          if (coord1 >= metrics->num_vector_active) {
            continue;
          }
          int num_out_this_line = 0;
          for (int i0 = 0; i0 < 2; ++i0) {
            for (int i1 = 0; i1 < 2; ++i1) {
              const GMFloat value = GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC ?
                GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env) :
                GMMetrics_duo_get_from_index_2(metrics, index, i0, i1, env);
              if (!(threshold < 0. || value > threshold)) {
                continue;
              }

              /*---Output the value---*/
              if (file == stdout) {

                if (num_out_this_line == 0) {
                  fprintf(file,
                    "element (%li,%li): values:", coord0, coord1);
                }

                fprintf(file,
                  sizeof(GMFloat) == 8 ?
                  " %i %i %.17e" :
                  " %i %i %.8e", i0, i1, value);

              } else {

                writer.write(coord0, coord1, i0, i1, value);

              }

              num_out_this_line++;

            } /*---i1---*/
          } /*---i0---*/
          if (file == stdout) {
            if (num_out_this_line > 0) {
              fprintf(file, "\n");
            }
          }
        } /*---index---*/

        num_written += writer.get_num_written();

      } /*---if stdout---*/

    } break;

    /*--------------------*/
    case GM_DATA_TYPE_TALLY4X2: {
    /*--------------------*/

      if (file != stdout) {

        output_metrics_tally4x2_bin_(metrics, file, threshold,
                                     num_written, env);

      } else /*---stdout---*/ {

        GMInsist(GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC);

        MetricWriter writer(file, metrics, env);

        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t coord0 =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t coord1 =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          const size_t coord2 =
            GMMetrics_coord_global_from_index(metrics, index, 2, env);
          if (coord0 >= metrics->num_vector_active ||
              coord1 >= metrics->num_vector_active ||
              coord2 >= metrics->num_vector_active) {
            continue;
          }
          int num_out_this_line = 0;
          for (int i0 = 0; i0 < 2; ++i0) {
            for (int i1 = 0; i1 < 2; ++i1) {
              for (int i2 = 0; i2 < 2; ++i2) {
                const GMFloat value
                  = GMMetrics_ccc_get_from_index_3(metrics, index, i0, i1, i2,
                                                env);
                if (!(threshold < 0. || value > threshold)) {
                  continue;
                }

                /*---Output the value---*/
                if (file == stdout) {

                  if (num_out_this_line == 0) {
                    fprintf(file,
                      "element (%li,%li,%li): values:", coord0, coord1, coord2);
                  }

                  fprintf(file,
                    sizeof(GMFloat) == 8 ?
                    " %i %i %i %.17e" :
                    " %i %i %i %.8e", i0, i1, i2, value);

                } else {

                  writer.write(coord0, coord1, coord2, i0, i1, i2, value);

                }

                num_out_this_line++;

              } /*---i2---*/
            } /*---i1---*/
          } /*---i0---*/
          if (file == stdout) {
            if (num_out_this_line > 0) {
              fprintf(file, "\n");
            }
          }
        } /*---index---*/

        num_written += writer.get_num_written();

      } /*---if stdout---*/

    } break;
    /*--------------------*/
    default:
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
}

//=============================================================================
// MetricWriter member definitions.

MetricWriter::MetricWriter(FILE* file, GMMetrics* metrics, GMEnv* env) :
  file_(file),
  data_type_(GMEnv_data_type_metrics(env)),
  num_way_(GMEnv_num_way(env)),
  num_written_total_(0) {

  if (stdout != file_) {
    GMInsistInterface(env, metrics->num_vector_active ==
                  (uint32_t)metrics->num_vector_active &&
                  "Too many vectors for output format.");
  }
}

//-----------------------------------------------------------------------------
// CZEK 2-way

void
MetricWriter::write(size_t coord0, size_t coord1, GMFloat value) {

  bool success = true;

  const uint32_t outc0 = coord0;
  size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
  success = success && num_written == 1;

  const uint32_t outc1 = coord1;
  num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
  success = success && num_written == 1;

  const GMFp32 outv = value;
  num_written = fwrite(&outv, sizeof(outv), 1, file_);
  success = success && num_written == 1;
//printf("%i %i %f\n", (int)outc0, (int)outc1, (double)outv);

  num_written_total_ += success ? 1 : 0;
  GMInsist(success && "File write failure.");
  }

//-----------------------------------------------------------------------------
// CZEK 3-way

void
MetricWriter::write(size_t coord0, size_t coord1, size_t coord2,
                     GMFloat value) {

  bool success = true;

  const uint32_t outc0 = coord0;
  size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
  success = success && num_written == 1;

  const uint32_t outc1 = coord1;
  num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
  success = success && num_written == 1;

  const uint32_t outc2 = coord2;
  num_written = fwrite(&outc2, sizeof(outc2), 1, file_);
  success = success && num_written == 1;

  const GMFp32 outv = value;
  num_written = fwrite(&outv, sizeof(outv), 1, file_);
  success = success && num_written == 1;

  num_written_total_ += success ? 1 : 0;
  GMInsist(success && "File write failure.");
}

//-----------------------------------------------------------------------------
// CCC, DUO 2-way

void
MetricWriter::write(size_t coord0, size_t coord1, int i0, int i1,
                     GMFloat value) {
  GMAssert(i0 >=0 && i0 < 2);
  GMAssert(i1 >=0 && i1 < 2);

  bool success = true;

  const uint32_t outc0 = i0 + 2 * coord0;
  size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
  success = success && num_written == 1;

  const uint32_t outc1 = i1 + 2 * coord1;
  num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
  success = success && num_written == 1;

  const GMFp32 outv = value;
  num_written = fwrite(&outv, sizeof(outv), 1, file_);
  success = success && num_written == 1;

  num_written_total_ += success ? 1 : 0;
  GMInsist(success && "File write failure.");
}

//-----------------------------------------------------------------------------
// CCC 3-way

void
MetricWriter::write(size_t coord0, size_t coord1, size_t coord2,
              int i0, int i1, int i2, GMFloat value) {
  GMAssert(i0 >=0 && i0 < 2);
  GMAssert(i1 >=0 && i1 < 2);
  GMAssert(i2 >=0 && i2 < 2);

  bool success = true;

  const uint32_t outc0 = i0 + 2 * coord0;
  size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
  success = success && num_written == 1;

  const uint32_t outc1 = i1 + 2 * coord1;
  num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
  success = success && num_written == 1;

  const uint32_t outc2 = i2 + 2 * coord2;
  num_written = fwrite(&outc2, sizeof(outc2), 1, file_);
  success = success && num_written == 1;

  const GMFp32 outv = value;
    num_written = fwrite(&outv, sizeof(outv), 1, file_);
  success = success && num_written == 1;

  num_written_total_ += success ? 1 : 0;
  GMInsist(success && "File write failure.");
}

//=============================================================================
// MetricsFile member definitions.

FILE* gm_metrics_file_open(char* metrics_file_path_stub, GMEnv* env) {

  // Form filename

  size_t len = strlen(metrics_file_path_stub);
  char* path = (char*)malloc((len+50) * sizeof(char));

  int num_digits = 0;
  for (int tmp = 1; ; tmp*=10, ++num_digits) {
    if (tmp > GMEnv_num_proc(env)) {
      break;
    }
  }

  char format[100];
  sprintf(format, "%s0%ii.bin", "%s_%", num_digits);

  sprintf(path, format, metrics_file_path_stub, GMEnv_proc_num(env));

  // Do open

  FILE* const file = fopen(path, "w");
  free(path);

  return file;
}

//-----------------------------------------------------------------------------


MetricsFile::MetricsFile(DriverOptions* do_, GMEnv* env)
  : file_(NULL)
  , num_written_(0) {
  GMInsist(do_ && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  verbosity_ = do_->verbosity;
  threshold_ = do_->threshold;

#if 0
  char* stub = do_->metrics_file_path_stub;
  if (! (NULL != stub && GMEnv_is_proc_active(env) &&
      GMEnv_proc_num_field(env) == 0) ) {
    return;
  }

  // Form filename

  size_t len = strlen(stub);
  char* path = (char*)malloc((len+50) * sizeof(char));

  int num_digits = 0;
  for (int tmp = 1; ; tmp*=10, ++num_digits) {
    if (tmp > GMEnv_num_proc(env)) {
      break;
    }
  }

  char format[100];
  sprintf(format, "%s0%ii.bin", "%s_%", num_digits);

  sprintf(path, format, stub, GMEnv_proc_num(env));

  /*---Do open---*/

  file_ = fopen(path, "w");
  GMInsist(NULL != file_ && "Unable to open file.");
  free(path);
#endif

  if (do_->metrics_file_path_stub) {
    file_ = gm_metrics_file_open(do_->metrics_file_path_stub, env);
    GMInsist(NULL != file_ && "Unable to open file.");
  }
}

//-----------------------------------------------------------------------------

MetricsFile::~MetricsFile() {
  if (file_) {
    fclose(file_);
  }
}

//-----------------------------------------------------------------------------

void MetricsFile::write(GMMetrics* metrics, GMEnv* env) {

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  // Output to file

  if (file_) {
    output_metrics_(metrics, file_, threshold_, num_written_, env);
  }

  // Output to stdout if requested

  if (verbosity_ > 1) {
    double threshold = verbosity_ > 2 ? -1. : threshold_;
    output_metrics_(metrics, stdout, threshold, num_written_, env);
  }
}

//=============================================================================
