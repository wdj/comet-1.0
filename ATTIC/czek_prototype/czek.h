/****************************************************************/
/* Czekanowski Similarity Metric                                */
/*                                                              */
/* Code Author: Doug Hyatt                                      */
/* Created: June, 2015                                          */
/****************************************************************/

#ifndef _CZEK_H_
#define _CZEK_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VECCHUNK 100
#define MAXFIELDSIZE 10000

#ifndef NCOPIES_V
#define NCOPIES_V 1
#endif

#ifndef NCOPIES_F
#define NCOPIES_F 1
#endif

#ifndef PROCESS_VECTORS_FN
#define PROCESS_VECTORS_FN process_vectors
#endif

#define NO_PRINT

typedef double Float_t;
enum{MPI_Float_t = MPI_DOUBLE};

struct _vector {
  char id[MAXFIELDSIZE];
  Float_t *data;
};

Float_t czekanowski(int, Float_t *, Float_t *);
void process_vectors(struct _vector *, int, int, FILE *);
void process_vectors_alt(struct _vector *, int, int, FILE *);
void process_vectors_alt2(struct _vector *, int, int, FILE *);
void process_vectors_alt3(struct _vector *, int, int, FILE *);
void process_vectors_alt4(struct _vector *, int, int, FILE *);
void process_vectors_alt5(struct _vector *, int, int, FILE *);
void process_vectors_alt6(struct _vector *, int, int, FILE *);
struct _vector *read_vectors(FILE *, int *, int *);
void free_vectors(struct _vector *, int);

#endif
