/****************************************************************/
/* Czekanowski Similarity Metric                                */
/*                                                              */
/* Code Author: Doug Hyatt                                      */
/* Created: June, 2015                                          */
/****************************************************************/

#include <sys/time.h>

#include "mpi.h"

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "czek.h"

/*===========================================================================*/
/* Timer function */

double get_time() {
  struct timeval tv;
  double result;
  gettimeofday( &tv, NULL );
  result = ( (double) tv.tv_sec +
             (double) tv.tv_usec * 1.e-6 );
  return result;
}

/*===========================================================================*/
/* Czekanowski Similarity Metric between Two Vectors */

/* v1, v2 = double vectors, len = the size of the vectors */
Float_t czekanowski(int len, Float_t * v1, Float_t * v2) {
  int i = 0;
  Float_t numerator = 0.0;
  Float_t denominator = 0.0;

  for(i = 0; i < len; i++) {
    numerator += (v1[i] < v2[i]?v1[i]:v2[i]);
    denominator += (v1[i] + v2[i]);
  }

  return 2.0*numerator/denominator;
}

/*===========================================================================*/
/* Czekanowski Similarity Metric between Two Vectors.
   Variant with restrict keyword */

Float_t czekanowski_alt(int len, Float_t * const __restrict__ v1,
                                 Float_t * const __restrict__ v2) {
  int i = 0;
  Float_t numerator = 0;
  Float_t denominator = 0;

  for(i = 0; i < len; ++i) {
    numerator += ( v1[i] < v2[i] ? v1[i] : v2[i] );
    denominator += ( v1[i] + v2[i] );
  }

  return ((Float_t)2)*numerator/denominator;
}

/*===========================================================================*/
/* Inlinable min function for czek metric. */
///Test comment blah blah blah
Float_t min_op(Float_t a, Float_t b) {
  return a < b ? a : b;
}

/*===========================================================================*/
/* Compute the numerator for the czek metric. */

Float_t czekanowski_numerator(int len, Float_t * const __restrict__ v1, Float_t * const __restrict__ v2) {
  int i = 0;
  Float_t result = (Float_t)0;

  for(i = 0; i < len; ++i) {
    result += min_op(v1[i], v2[i]);
  }

  return result;
}

/*===========================================================================*/
/* Compute simple sum of elements of a vector. */

Float_t vector_sum(int len, Float_t * const __restrict__ v1) {
  int i = 0;
  Float_t result = 0;

  for(i = 0; i < len; ++i) {
    result += ( v1[i] );
  }

  return result;
}

/*===========================================================================*/
/* Compute a checksum of the czek result values. */

Float_t checksum( Float_t* czek_vals, int numvec) {
  int i = 0;
  int j = 0;
  Float_t result = 0;

#if 0
  for (i = 0; i < (numvec/NCOPIES_V)-1; ++i) {
    for (j = i+1; j < (numvec/NCOPIES_V); ++j) {
      result += (1 + j + (numvec/NCOPIES_V)*1.*(i)) *
                czek_vals[j+numvec*i];
    }  
  }
#endif

#if 1
  /*---Compute weighted sum of values---*/
  for (i = 0; i < (numvec)-1; ++i) {
    for (j = i+1; j < (numvec); ++j) {
      result += (1 + j + (numvec)*1.*(i)) *
                czek_vals[j+numvec*i];
    }  
  }
#endif

  /*---Sum the result across nodes if appropriate---*/
  Float_t tmp = result;
  MPI_Allreduce(&tmp, &result, 1, MPI_Float_t, MPI_SUM, MPI_COMM_WORLD);

  return result;
}

/*===========================================================================*/
/*===========================================================================*/
/* Process all pairwise combinations of vectors (We do czek(i,j) and */
/* not czek(j,i) since the result is the same. */

void process_vectors(struct _vector *vectors, int numvec, int numfield,
                     FILE *fp) {
  int i = 0;
  int j = 0;
  Float_t czek = 0.0;   /* Value for the Czekanowski Result */

  for (i = 0; i < numvec-1; i++) {
    for (j = i+1; j < numvec; j++) {
      czek = czekanowski(numfield, vectors[i].data, vectors[j].data);
      fprintf(fp, "%s\t%.4f\t%s\n", vectors[i].id, czek, vectors[j].id);
    }  
  }
}

/*===========================================================================*/
/* Alternate computation of Czekanowski metric for all pairwise combinations.
   For this version, time the computations, without counting I/O. */

void process_vectors_alt(struct _vector *vectors, int numvec, int numfield,
                     FILE *fp) {
  int i = 0;
  int j = 0;
  Float_t* czek_vals = 0;
  double time1 = 0.0;
  double time2 = 0.0;

  czek_vals = malloc(numvec*numvec*sizeof(Float_t));

  time1 = get_time();

  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      const Float_t czek =
                       czekanowski(numfield, vectors[i].data, vectors[j].data);
      czek_vals[j+numvec*i] = czek;
    }  
  }

  time2 = get_time();

#ifndef NO_PRINT
  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      Float_t czek = czek_vals[j+numvec*i];
      fprintf(fp, "%s\t%.4f\t%s\n", vectors[i].id, czek, vectors[j].id);
    }  
  }
#endif

  /*---Output result summary on proc 0---*/
  Float_t cksum = checksum(czek_vals, numvec);
  int comm_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  if (comm_rank == 0) {
    printf("alt  numvec %i numfield %i "
           "time: %.6f "
           "checksum %.15e\n",
           numvec, numfield, time2-time1, cksum);
  }

  free(czek_vals);
}

/*===========================================================================*/
/* Alternate computation of Czekanowski metric for all pairwise combinations.
   For this version, change the ordering of computations, memory layout to
   attempt to improve speed. (wasn't very effective) */

void process_vectors_alt2(struct _vector *vectors, int numvec, int numfield,
                     FILE *fp) {
  int i = 0;
  int j = 0;
  int k = 0;
  Float_t* czek_vals = 0;
  Float_t* czek_numerators = 0;
  Float_t* czek_denominators = 0;
  Float_t* vectors_reordered = 0;
  double time1 = 0.0;
  double time2 = 0.0;

  czek_vals = malloc(numvec*numvec*sizeof(Float_t));
  czek_numerators = malloc(numvec*numvec*sizeof(Float_t));
  czek_denominators = malloc(numvec*numvec*sizeof(Float_t));
  vectors_reordered = malloc(numvec*numfield*sizeof(Float_t));

  /* Copy vectors into reordered array */

  for (i = 0; i < numvec; ++i) {
    for (k = 0; k < numfield; ++k) {
      vectors_reordered[i+numvec*k] = vectors[i].data[k];
    }  
  }

  time1 = get_time();

  /* Initialize */

  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      czek_numerators[j+numvec*i] = 0.;
      czek_denominators[j+numvec*i] = 0.;
    }  
  }

  /* Do Czek computation */

  for (k = 0; k < numfield; ++k) {
    const Float_t* const __restrict__ vector_field_k
                                            = &(vectors_reordered[0+numvec*k]);
    for (i = 0; i < numvec-1; ++i) {
      Float_t* const __restrict__ czek_num_i = &(czek_numerators[0+numvec*i]);
      Float_t* const __restrict__ czek_denom_i = &(czek_denominators[0+numvec*i]);
      const Float_t vik = vector_field_k[i];
      for (j = i+1; j < numvec; ++j) {
        const Float_t vjk = vector_field_k[j];
        czek_num_i[j] += vik < vjk ? vik : vjk;
        czek_denom_i[j] += (vik + vjk);
      }  
    }
  }

  /* Get final result */

  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      czek_vals[j+numvec*i] = ((Float_t)2.0) *
          czek_numerators[j+numvec*i] /
          czek_denominators[j+numvec*i];
    }  
  }

  time2 = get_time();

  /* Output */

#ifndef NO_PRINT
  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      const Float_t czek = czek_vals[j+numvec*i];
      fprintf(fp, "%s\t%.4f\t%s\n", vectors[i].id, czek, vectors[j].id);
    }  
  }
#endif

  /*---Output result summary on proc 0---*/
  Float_t cksum = checksum(czek_vals, numvec);
  int comm_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  if (comm_rank == 0) {
    printf("alt  numvec %i numfield %i "
           "time: %.6f "
           "checksum %.15e\n",
           numvec, numfield, time2-time1, cksum);
  }

  free(czek_vals);
  free(czek_numerators);
  free(czek_denominators);
  free(vectors_reordered);
}

/*===========================================================================*/
/* Alternate computation of Czekanowski metric for all pairwise combinations.
   For this version, apply blocking of vectors, to try to get better
   reuse of elements read into cache. (helped some) */

void process_vectors_alt3(struct _vector *vectors, int numvec, int numfield,
                     FILE *fp) {
  int i = 0;
  int ibase = 0;
  int j = 0;
  Float_t czek = (Float_t)0;
  Float_t* __restrict__ czek_vals = 0;
  double time1 = 0.0;
  double time2 = 0.0;
  const int blocksize = 8;

  czek_vals = malloc(numvec*numvec*sizeof(Float_t));

  time1 = get_time();

  for (ibase = 0; ibase < numvec-1; ibase+=blocksize) {
    for (j = ibase+1; j < numvec; ++j) {
      for (i=ibase; i<ibase+blocksize && i<j; ++i) {
        czek = czekanowski_alt(numfield, vectors[i].data, vectors[j].data);
        czek_vals[j+numvec*i] = czek;
      }  
    }  
  }

  time2 = get_time();

#ifndef NO_PRINT
  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      const Float_t czek = czek_vals[j+numvec*i];
      fprintf(fp, "%s\t%.4f\t%s\n", vectors[i].id, czek, vectors[j].id);
    }  
  }
#endif

  /*---Output result summary on proc 0---*/
  Float_t cksum = checksum(czek_vals, numvec);
  int comm_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  if (comm_rank == 0) {
    printf("alt  numvec %i numfield %i "
           "time: %.6f "
           "checksum %.15e\n",
           numvec, numfield, time2-time1, cksum);
  }

  free(czek_vals);
}

/*===========================================================================*/
/* Alternate computation of Czekanowski metric for all pairwise combinations.
   For this version, separate out the denominator computation, compute
   using individual vector sums. (helped some) */

void process_vectors_alt4(struct _vector *vectors, int numvec, int numfield,
                     FILE *fp) {
  int i = 0;
  int j = 0;
  Float_t czek_numerator = (Float_t)0;
  Float_t* czek_vals = 0;
  Float_t* vector_sums = 0;
  double time1 = 0.0;
  double time2 = 0.0;

  czek_vals = malloc(numvec*numvec*sizeof(Float_t));

  vector_sums = malloc(numvec*sizeof(Float_t));

  time1 = get_time();

  for (i = 0; i < numvec; ++i) {
      vector_sums[i] = vector_sum(numfield, vectors[i].data);
  }

  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      czek_numerator =
             czekanowski_numerator(numfield, vectors[i].data, vectors[j].data);
      czek_vals[j+numvec*i] = ((Float_t)2) * czek_numerator /
                                          ( vector_sums[i] + vector_sums[j] );
    }  
  }

  time2 = get_time();

#ifndef NO_PRINT
  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      const Float_t czek = czek_vals[j+numvec*i];
      fprintf(fp, "%s\t%.4f\t%s\n", vectors[i].id, czek, vectors[j].id);
    }  
  }
#endif

  /*---Output result summary on proc 0---*/
  Float_t cksum = checksum(czek_vals, numvec);
  int comm_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  if (comm_rank == 0) {
    printf("alt  numvec %i numfield %i "
           "time: %.6f "
           "checksum %.15e\n",
           numvec, numfield, time2-time1, cksum);
  }

  free(czek_vals);
  free(vector_sums);
}

/*===========================================================================*/
/* Compute the strict triangular part of A^T * A, with scalar multiplication
   replaced with the "min" operation. */

void matrix_matrix_min_product(int numvec, int numfield,
  Float_t* const __restrict__ vector_matrix,
  Float_t* const __restrict__ czek_vals) {

  int i = 0;
  int j = 0;
  int k = 0;
  int ibase = 0;
  const int blocksize = 6;

  for (ibase = 0; ibase < numvec-1; ibase+=blocksize) {
    for (j = ibase+1; j < numvec; ++j) {
      for (i=ibase; i<ibase+blocksize && i<j; ++i) {
        Float_t czek_numerator = (Float_t)0;
        for(k = 0; k < numfield; ++k) {
          czek_numerator += min_op(vector_matrix[k+numfield*i],
                                   vector_matrix[k+numfield*j]);
        }
        czek_vals[j+numvec*i] = czek_numerator;
      }
    }
  }
}

/*===========================================================================*/
/* Alternate computation of Czekanowski metric for all pairwise combinations.
   For this version, separate the numerator calculation and use a
   matrix-matrix-product-like computation. */

void process_vectors_alt5(struct _vector *vectors, int numvec, int numfield,
                     FILE *fp) {
  int i = 0;
  int j = 0;
  int k = 0;
  Float_t* __restrict__ czek_vals = 0;
  Float_t* __restrict__ vector_sums = 0;
  Float_t* __restrict__ vector_matrix = 0;
  double time1 = 0.0;
  double time2 = 0.0;

  czek_vals = malloc(numvec*numvec*sizeof(Float_t));
  vector_sums = malloc(numvec*sizeof(Float_t));
  vector_matrix = malloc(numvec*numfield*sizeof(Float_t));

  for (i = 0; i < numvec; ++i) {
    for (k = 0; k < numfield; ++k) {
      vector_matrix[k+numfield*i] = vectors[i].data[k];
    }  
  }

  time1 = get_time();

  for (i = 0; i < numvec; ++i) {
      vector_sums[i] = vector_sum(numfield, vectors[i].data);
  }

  matrix_matrix_min_product(numvec, numfield, vector_matrix, czek_vals);

  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      czek_vals[j+numvec*i] = ((Float_t)2) * czek_vals[j+numvec*i] /
                                          ( vector_sums[i] + vector_sums[j] );
    }  
  }

  time2 = get_time();

#ifndef NO_PRINT
  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      const Float_t czek = czek_vals[j+numvec*i];
      fprintf(fp, "%s\t%.4f\t%s\n", vectors[i].id, czek, vectors[j].id);
    }  
  }
#endif

  /*---Output result summary on proc 0---*/
  Float_t cksum = checksum(czek_vals, numvec);
  int comm_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  if (comm_rank == 0) {
    printf("alt  numvec %i numfield %i "
           "time: %.6f "
           "checksum %.15e\n",
           numvec, numfield, time2-time1, cksum);
  }

  free(czek_vals);
  free(vector_sums);
  free(vector_matrix);
}

/*===========================================================================*/
/* Alternate computation of Czekanowski metric for all pairwise combinations.
   For this version, use modified magma library.

   The magma library was modified as follows.  The "daxpy" function near
   the beginning of the file magma-1.6.2/magmablas/dgemm_T_N.cu
   was modified to replace the multiplies with a "min" operation.

   This version uses the code for the older NVIDIA Tesla GPUS.  It is likely
   that other code in magma e.g. for Fermi would give higher performance.
 */

void process_vectors_alt6(struct _vector *vectors, int numvec, int numfield,
                     FILE *fp) {
  int i = 0;
  int j = 0;
  int k = 0;
  Float_t* __restrict__ vector_sums = 0;
  Float_t* czek_vals = 0;
  Float_t* vector_matrix = 0;
  double time1 = 0.0;
  double time2 = 0.0;

  vector_sums = malloc(numvec*sizeof(Float_t));

  magma_minproduct_init();

  /* Allocate magma CPU memory for vectors and for czek_vals result */
  magma_minproduct_dmalloc_pinned(&czek_vals,numvec*numvec);
  magma_minproduct_dmalloc_pinned(&vector_matrix,numvec*numfield);

  /* Copy in vectors */

  for (i = 0; i < numvec; ++i) {
    for (k = 0; k < numfield; ++k) {
      vector_matrix[k+numfield*i] = vectors[i].data[k];
    }  
  }

  /* Allocate GPU mirrors for CPU arrays */

  Float_t* d_vector_matrix = 0;
  Float_t* d_czek_vals = 0;

  magma_minproduct_dmalloc(&d_vector_matrix, numvec*numfield);
  magma_minproduct_dmalloc(&d_czek_vals, numvec*numvec);

  /* Initialize result to zero (apparently is required) */
  for (i = 0; i < numvec; ++i) {
    for (j = 0; j < numvec; ++j) {
      czek_vals[j+numvec*i] = 0;
    }
  }

  /* Send matrix to GPU */

  magma_minproduct_dsetmatrix(numvec, numvec, czek_vals, numvec,
                                 d_czek_vals, numvec);

/*
  magma_minproduct_dgetmatrix(numvec, numvec, d_czek_vals, numvec, 
                                     czek_vals, numvec);
*/

  time1 = get_time();

  /* Calculate vector sums for denom */

  for (i = 0; i < numvec; ++i) {
      vector_sums[i] = vector_sum(numfield, vectors[i].data);
  }

  /* Send matrix to GPU */

  magma_minproduct_dsetmatrix(numfield, numvec, vector_matrix, numfield,
                                   d_vector_matrix, numfield);

  /* Perform pseudo matrix-matrix product */

  magma_minproductblas_dgemm(
    Magma_minproductTrans,
    Magma_minproductNoTrans,
    numvec,
    numvec,
    numfield,
    1.0,
    d_vector_matrix,
    numfield,
    d_vector_matrix,
    numfield,
    0.0,
    d_czek_vals,
    numvec );

  /* Copy result from GPU */

  magma_minproduct_dgetmatrix(numvec, numvec, d_czek_vals, numvec, 
                                     czek_vals, numvec);

/*
  matrix_matrix_min_product(numvec, numfield, vector_matrix, czek_vals);
*/

  /* Compute final result */

  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      czek_vals[j+numvec*i] = ((Float_t)2) * czek_vals[j+numvec*i] /
                                          ( vector_sums[i] + vector_sums[j] );
    }  
  }

  time2 = get_time();

#ifndef NO_PRINT
  for (i = 0; i < numvec-1; ++i) {
    for (j = i+1; j < numvec; ++j) {
      const Float_t czek = czek_vals[j+numvec*i];
      fprintf(fp, "%s\t%.4f\t%s\n", vectors[i].id, czek, vectors[j].id);
    }  
  }
#endif

  /*---Output result summary on proc 0---*/
  Float_t cksum = checksum(czek_vals, numvec);
  int comm_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  if (comm_rank == 0) {
    printf("alt  numvec %i numfield %i "
           "time: %.6f "
           "checksum %.15e\n",
           numvec, numfield, time2-time1, cksum);
  }

  /* Free memory */

  magma_minproduct_free(d_czek_vals);
  magma_minproduct_free(d_vector_matrix);

  magma_minproduct_free_pinned(czek_vals);
  magma_minproduct_free_pinned(vector_matrix);

  magma_minproduct_finalize();

  free(vector_sums);
}

/*===========================================================================*/
/*===========================================================================*/
/* Read in the vectors (ids and data) and dynamically allocate */
/* memory as needed. */

struct _vector *read_vectors(FILE *fp, int *numvec, int *numfield) {
  int i = 0;
  int c = 0;               /* Character we read in */
  int field_ctr = 0;       /* Count the number of tab-delimited fields */
  int vec_ctr = 0;         /* Count the number of vectors (lines) */
  int ch_ctr = 0;          /* Counter for individual chars */
  char buf[MAXFIELDSIZE] = "";    /* Buffer for reading fields */
  char *conv_ptr = NULL;   /* Pointer to check conversion */
  struct _vector *vectors = NULL; /* Pointer to the vectors */

  int numvec_orig = 0;
  int numfield_orig = 0;

  /* Read in the first line to determine the number of fields */
  do {
    c = getc(fp);
    if (c == '\t') field_ctr++;
  } while(c != EOF && c != '\n');
  numfield_orig = field_ctr;
  *numfield = numfield_orig * NCOPIES_F;

  /* Now handle various error situations */
  if (feof(fp) != 0) {
    fprintf(stderr, "Read error: Saw end of file before expected.\n");
    return NULL;
  }
  if (ferror(fp) != 0 || c != '\n') {
    fprintf(stderr, "Read error: Error reading the file.\n");
    return NULL;
  }

  /* Allocate initial memory for the vectors and ids */
  vectors = (struct _vector *)malloc(VECCHUNK * sizeof(struct _vector) * NCOPIES_V);
  if (vectors == NULL) {
    fprintf(stderr, "Failed to allocate memory for the vectors.\n");
    return NULL;
  }
  for (i = 0; i < VECCHUNK; i++)
    memset(&vectors[i], 0, sizeof(struct _vector));

  /* Read in the remaining lines of the file and add data to the vectors */
  while(c != EOF) {

    /* Allocate memory for a new vector */
    vectors[vec_ctr].data = (Float_t *)malloc((*numfield)*sizeof(Float_t));
    if (vectors[vec_ctr].data == NULL) {
      fprintf(stderr, "Failed to allocate memory for the vectors.\n");
      free_vectors(vectors, vec_ctr);
      return NULL;
    }

    /* Reset field/character counters */
    field_ctr = 0;
    ch_ctr = 0;

    /* Read characters until we hit a tab or newline, then process */
    /* the buffer, either adding it to ID (field 1) or vector data. */
    do {
      c = getc(fp);
      if (c == '\t' || c == '\n') {
        if (field_ctr == 0) vectors[vec_ctr].id[ch_ctr] = '\0';
        else if (field_ctr <= numfield_orig) {
          buf[ch_ctr] = '\0';
          vectors[vec_ctr].data[field_ctr-1] = (Float_t)strtod(buf, &conv_ptr);
          if (conv_ptr != buf+ch_ctr) {
            fprintf(stderr, "Error converting '%s' to double, line %d.\n",
                    buf, vec_ctr+2);
            free_vectors(vectors, vec_ctr);
            return NULL;
          }
        }
        if (c == '\t') field_ctr++; 
        ch_ctr = 0;
      }
      else if (field_ctr == 0)
        vectors[vec_ctr].id[ch_ctr++] = c;
      else buf[ch_ctr++] = c;
    } while(c != EOF && c != '\n');
    if (c == EOF) break;

    /* If we see wrong number of fields... */
    if (field_ctr != numfield_orig) {
      fprintf(stderr, "Did not see correct number of fields.\n");
      fprintf(stderr, "Expected %d, saw %d, line %d.\n", numfield_orig, field_ctr, vec_ctr+2);
      free_vectors(vectors, vec_ctr);
      return NULL;
    }

    /* Increment vector counter.  If hit our limit, realloc for more. */
    vec_ctr++;
    if (vec_ctr%VECCHUNK == 0) {
      vectors = (struct _vector *)realloc(vectors,
                (vec_ctr+VECCHUNK) * sizeof(struct _vector) * NCOPIES_V);
      if (vectors == NULL) {
        free_vectors(vectors, vec_ctr);
        fprintf(stderr, "Failed to realloc memory for vectors.\n");
        return NULL;
      }
      for (i = vec_ctr; i < vec_ctr + VECCHUNK; i++)
        memset(&vectors[i], 0, sizeof(struct _vector));
    }
  }
  if (ferror(fp) != 0) {
    fprintf(stderr, "Read error: Error reading the file.\n");
    free_vectors(vectors, vec_ctr);
    return NULL;
  }

  numvec_orig = vec_ctr;
  *numvec = numvec_orig * NCOPIES_V;

  /* Replicate fields */
  for (i=0; i<numvec_orig; ++i) {
    int k = 0;
    for (k=0; k<numfield_orig; ++k) {
      int j = 0;
      for (j=1; j<NCOPIES_F; ++j) {
        vectors[i].data[k+j*numfield_orig] = vectors[i].data[k];
      }
    }
  }

  /* Replicate vectors */
  for (i=0; i<numvec_orig; ++i) {
    int j = 0;
    for (j=1; j<NCOPIES_V; ++j) {
      int k = 0;
      strcpy(vectors[i+j*numvec_orig].id, vectors[i].id);
      vectors[i+j*numvec_orig].data = (Float_t *)malloc((*numfield)*sizeof(Float_t));
      for (k=0; k<*numfield; ++k) {
        vectors[i+j*numvec_orig].data[k] = vectors[i].data[k];
      }
    }
  }

  return vectors;
}

/* Free routine for vectors */
void free_vectors(struct _vector *vectors, int numvec) {
  int i = 0;

  if (vectors != NULL) {
    for (i = 0; i < numvec; i++)
      if (vectors[i].data != NULL) free(vectors[i].data);
    free(vectors);
  } 
}

/*===========================================================================*/
