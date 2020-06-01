//-----------------------------------------------------------------------------
// "Manually" calculate a CCC metric value directly from the original
// SNP file.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

enum{MAX_LABEL_LEN = 11};

//-----------------------------------------------------------------------------

int process_line(int argc, char** argv, FILE* snptxtfile, FILE* lifile) {

  if (argc != 1+4 && argc != 1+6) {
    // Input format expectd from each line of stdin:
    printf("Line format: lineno0 bitno0 lineno1 bitno1 [lineno2 bitno2]\n");
    printf("Line and bit numbers are 0-based\n");
    return 1;
  }

  const int num_way = argc == 1+4 ? 2 : 3;
  enum {NUM_WAY_MAX = 3};

  int lineno[NUM_WAY_MAX];
  lineno[0] = atoi(argv[1]);
  lineno[1] = atoi(argv[3]);
  lineno[2] = num_way == 3 ? atoi(argv[5]) : 0;

  int bitno[NUM_WAY_MAX];
  bitno[0] = atoi(argv[2]);
  bitno[1] = atoi(argv[4]);
  bitno[2] = num_way == 3 ? atoi(argv[6]) : 0;

  const int no_al = 'X';

  const int label_len = MAX_LABEL_LEN;
  unsigned char label[NUM_WAY_MAX][label_len+1];
  for (int i=0; i<num_way; ++i) {
    for (int j=0; j<label_len+1; ++j) {
      label[i][j] = 0;
    }
  }

  int num_read = 0;
  int fseek_success = 0;

  // Initial loop to get upper bound on fields per line.

  fseek_success = fseek(snptxtfile, 0, SEEK_SET);
  if (0 != fseek_success) {
    printf("Error: error reading SNP data file (0).\n");
    return 1;
  }

  int c = 0;

  while ((c = fgetc(snptxtfile)) != '\n') {
    if (c == EOF) {
      printf("Error snp file has no newlines.");
      return 1;
    }
    num_read++;
  }

  const int num_field_max = num_read / 2 + 1;

  // Storage for allelel labels and alleles for each needed SNP.

  int al[NUM_WAY_MAX][2];
  int* elts[NUM_WAY_MAX];
  for (int i=0; i<num_way; ++i) {
    al[i][0] = no_al;
    al[i][1] = no_al;
    elts[i] = malloc(num_field_max * 2 * sizeof(int));
  }

  int num_field = 0;

  // Loop over num_way to input the required vectors.

  for (int i=0; i<num_way; ++i) {

    fseek_success = fseek(lifile, lineno[i]*sizeof(size_t), SEEK_SET);
    if (0 != fseek_success) {
      printf("Error: error reading label_indices file (1).\n");
      return 1;
    }

    // Get the index to this line in the file

    size_t line_index = 0;
    num_read = fread(&line_index, sizeof(size_t), 1, lifile);
    if (1 != num_read) {
      printf("Error: error reading label_indices file (2).\n");
      return 1;
    }

    fseek_success = fseek(snptxtfile, line_index, SEEK_SET);
    if (0 != fseek_success) {
      printf("Error: error reading SNP data file (1).\n");
      return 1;
    }

    int num_tabs = 0;
    int j = 0;

    // Loop to read in line

    while ((c = fgetc(snptxtfile)) != '\n') {
      if (c == '\t') {
        num_tabs++;
        continue;
      }

      // Get line label

      if (num_tabs == 1) {
        label[i][j++] = c; //check
        continue;
      }

      if (num_tabs == 3) {
        j = 0;
      }

      if (num_tabs < 4) {
        continue;
      }

      // Store allele

      elts[i][j++] = c;
      if (i == 0) {
        num_field++;
      }

      // Record allele label

      if (c != '0') {
        if (al[i][0] == no_al) {
          al[i][0] = c;
        } else if (al[i][1] == no_al) {
          if (al[i][0] < c) {
            al[i][1] = c;
          } else if (al[i][0] > c) {
            al[i][1] = al[i][0];
            al[i][0] = c;
          }
        }
      } // if
    } // while

  } // for i

  num_field /= 2;

  // We account for sparsity of the data

 // First get sum_i's

  int count1[NUM_WAY_MAX];
  int sum1[NUM_WAY_MAX];

  for (int i=0; i<num_way; ++i) {
    count1[i] = 0;
    sum1[i] = 0;
    for (int f=0; f<num_field; ++f) {
      const int e0 = elts[i][2*f];
      const int e1 = elts[i][2*f+1];

      if (e0 == '0') {
        // assert(e1 == '0');
        continue;
      }

      count1[i] += 1;

      // Calculate row and add to sum - see paper

      const int rho = (e0 == al[i][bitno[i]]) + (e1 == al[i][bitno[i]]);
      sum1[i] += rho;
    } // for f
  } // for i

  // Now get sum_{ij}'s (or sum_{ijk}'s if 3-way)

  int countijk = 0;
  int sumijk = 0;

  for (int f=0; f<num_field; ++f) {
    const int e00 = elts[0][2*f];
    const int e01 = elts[0][2*f+1];
    if (e00 == '0') {
      continue;
    }

    const int e10 = elts[1][2*f];
    const int e11 = elts[1][2*f+1];
    if (e10 == '0') {
      continue;
    }

    const int e20 = num_way == 2 ? 1 :  elts[2][2*f];
    const int e21 = num_way == 2 ? 1 :  elts[2][2*f+1];
    if (num_way == 3 && e20 == '0') {
      continue;
    }

    countijk += 1;

    const int rho0 = (e00 == al[0][bitno[0]]) + (e01 == al[0][bitno[0]]);
    const int rho1 = (e10 == al[1][bitno[1]]) + (e11 == al[1][bitno[1]]);
    const int rho2 = num_way == 2 ? 1 :
                     (e20 == al[2][bitno[2]]) + (e21 == al[2][bitno[2]]);

    sumijk += rho0 * rho1 * rho2;
  } // for f

  // substitute into formula

  double f1[NUM_WAY_MAX];
  for (int i=0; i<num_way; ++i) {
    f1[i] = 0 == count1[i] ? 0 : sum1[i] * 1. / (double)( 2 * count1[i] );
  }

  const double fijk = 0 == countijk ? 0 :
                     sumijk * 1. / (double)( (1 << num_way) * countijk );

  // NOTE hard-wired constant here

  const double ccc_multiplier = 9. / (double)2.;
  const double ccc_param = 2. / (double)3.;

  const double value = 2 == num_way ?
    ccc_multiplier * fijk * (1 - ccc_param*f1[0]) *  (1 - ccc_param*f1[1]) :
    ccc_multiplier * fijk * (1 - ccc_param*f1[0]) *  (1 - ccc_param*f1[1])
                          * (1 - ccc_param*f1[2]);


  // Permute labels to output each result in a uniform order

  int perm[3];

  if (num_way == 2) {
    if (lineno[0] > lineno[1]) {
      perm[0] = 0;
      perm[1] = 1;
    } else {
      perm[0] = 1;
      perm[1] = 0;
    }
  } else {
    if (lineno[0] > lineno[1] && lineno[1] > lineno[2]) {
      perm[0] = 0;
      perm[1] = 1;
      perm[2] = 2;
    } else if (lineno[0] > lineno[2] && lineno[2] > lineno[1]) {
      perm[0] = 0;
      perm[1] = 2;
      perm[2] = 1;
    } else if (lineno[1] > lineno[0] && lineno[0] > lineno[2]) {
      perm[0] = 1;
      perm[1] = 0;
      perm[2] = 2;
    } else if (lineno[1] > lineno[2] && lineno[2] > lineno[0]) {
      perm[0] = 1;
      perm[1] = 2;
      perm[2] = 0;
    } else if (lineno[2] > lineno[0] && lineno[0] > lineno[1]) {
      perm[0] = 2;
      perm[1] = 0;
      perm[2] = 1;
    } else if (lineno[2] > lineno[1] && lineno[1] > lineno[0]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 0;
    }
  } // if num_way

  // Do output to stdout

  for (int i=0; i<num_way; ++i) {
    int iperm = perm[i];
    printf(0 != i ? " " : "");
    printf("%i %i", lineno[iperm], bitno[iperm]);
  } // for i

  for (int i=0; i<num_way; ++i) {
    int iperm = perm[i];
    printf(" %s_%c", label[iperm], al[iperm][bitno[iperm]]);
  } // for i

  printf(" %f\n", value);

  for (int i=0; i<num_way; ++i) {
    free(elts[i]);
  }
  return 0;
}

//-----------------------------------------------------------------------------

int main(int argc, char** argv) {

  if (argc < 4) {
    printf("ccc_validate: create validation data for CCC calculations\n");
    printf("Usage: ccc_validate <num_way> <snptxtfile> <line_index_file>\n");
    return 0;
  }

  if (!( ('2' == argv[1][0] || '3' == argv[1][0]) && 0 == argv[1][1] )) {
    printf("Error: invalid num_way\n");
    return 1;
  }

  const int num_way = atoi(argv[1]);

  FILE* snptxtfile = fopen(argv[2], "r");
  if (!snptxtfile) {
    printf("Error: unable to open file. %s\n", argv[2]);
    return 1;
  }

  FILE* lifile = fopen(argv[3], "rb");
  if (!lifile) {
    printf("Error: unable to open file. %s\n", argv[2]);
    return 1;
  }

  enum{LINE_LEN_MAX = 4096};
  unsigned char line[LINE_LEN_MAX];

  int lineptr = 0;
  int argc_ = 1;
  char* argv_[LINE_LEN_MAX];
  int c = 0;

  // Loop over chars from stdin

  while ((c = fgetc(stdin)) != EOF) {

    if (c != '\n') {
      line[lineptr++] = c;
      if (lineptr >= LINE_LEN_MAX) {
        printf("Error: line too long");
        return 1;
      }
      continue;
    }

    // Process full line
    // Create an arg list consisting of the tokens

    line[lineptr] = 0;
    argv_[0] = (char*)&line[lineptr];
    for (int i=0; i<lineptr; ++i) {
      if (line[i] == ' ') {
        line[i] = 0;
      }
    }
    for (int i=0; i<lineptr; ++i) {
      if (line[i] != 0 && (i == 0 || line[i-1] == 0)) {
          argv_[argc_++] = (char*)&line[i];
      }
    }

    // Only use the first several tokens, specifying line and bit numbers
    argc_ = argc_ < 2*num_way+1 ? argc_ : 2*num_way+1;
    int result = process_line(argc_, argv_, snptxtfile, lifile);
    if (result != 0) {
      return 1;
    }

    // Prepare for next line

    lineptr = 0;
    argc_ = 1;
  }

  fclose(snptxtfile);
  fclose(lifile);

  return 0;
}

//-----------------------------------------------------------------------------
