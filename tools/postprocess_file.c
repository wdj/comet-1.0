//-----------------------------------------------------------------------------
// Postprocess a single output file from CoMet - convert binary to text.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

enum{MAX_LABEL_LEN = 11};

//-----------------------------------------------------------------------------

int main(int argc, char** argv) {

  if (argc < 6) {
    printf("postprocess_file: convert metrics file from binary to text\n");
    printf("Usage: postprocess_file <num_way> <allele_label_file> <label_file> <metrics_bin_file> <metrics_txt_file>\n");
    return 0;
  }

  if (!( ('2' == argv[1][0] || '3' == argv[1][0]) && 0 == argv[1][1] )) {
    printf("Error: invalid num_way\n");
    return 1;
  }

  const int num_way = atoi(argv[1]);

  FILE* alfile = fopen(argv[2], "r");
  if (!alfile) {
    printf("Error: unable to open file. %s\n", argv[2]);
    return 1;
  }

  FILE* labelfile = fopen(argv[3], "r");
  if (!labelfile) {
    printf("Error: unable to open file. %s\n", argv[3]);
    return 1;
  }

  char* metricsbinfilename = argv[4];

  FILE* metricsbinfile = fopen(argv[4], "rb");
  if (!metricsbinfile) {
    printf("Error: unable to open file. %s\n", argv[4]);
    return 1;
  }

  FILE* metricstxtfile = fopen(argv[5], "w");
  if (!metricstxtfile) {
    printf("Error: unable to open file. %s\n", argv[5]);
    return 1;
  }

  // Initializations

  int num_read = 0;
  int fseek_success = 0;
  float value = 0;
  int coord[3];
  const int no_al = 'X'; // Marker for case of only one allele label in a line

  // Loop over CCC result values in specified file

  while ((num_read = fread(&coord[0], sizeof(int), 1, metricsbinfile)) == 1) {

    // Get (rest of) coordinates and metric value

    num_read = fread(&coord[1], sizeof(int), 1, metricsbinfile);
    if (1 != num_read) {
      printf("Error: error reading file. 1\n");
      return 1;
    }

    if (3 == num_way) {
      num_read = fread(&coord[2], sizeof(int), 1, metricsbinfile);
      if (num_read != 1) {
        printf("Error: error reading file. 2\n");
        return 1;
      }
    }

    num_read = fread(&value, sizeof(float), 1, metricsbinfile);
    if (1 != num_read) {
      printf("Error: error reading file. 3\n");
      return 1;
    }

    const int al_len = 2;
    const int label_len = MAX_LABEL_LEN;
    const int num_way_max = 3;

    int lineno[num_way_max];
    int bitno[num_way_max];
    unsigned char al[num_way_max][al_len];
    unsigned char label[num_way_max][label_len+1];

    // Loop over coordinates

    for (int i=0; i<num_way; ++i) {

      // Decode vector number (line number in file) and whether lo or hi bit

      lineno[i] = coord[i] / 2;
      bitno[i] = coord[i] % 2;

      // Get allele labels

      fseek_success = fseek(alfile, (al_len+1)*lineno[i], SEEK_SET);
      if (0 != fseek_success) {
        printf("Error: error reading file. 4\n");
        return 1;
      }

      num_read = fread(&al[i], sizeof(unsigned char), al_len, alfile);
      if (al_len != num_read) {
        printf("Error: error reading file. 5 %s\n", metricsbinfilename);
        return 1;
      }

      // If repeated then disregard upper value, replace with "X"
      if (al[i][0] == al[i][1]) {
        al[i][1] = no_al;
      }

      // Get (vector) label

      fseek_success = fseek(labelfile, (label_len+1)*lineno[i], SEEK_SET);
      if (0 != fseek_success) {
        printf("Error: error reading file. 6\n");
        return 1;
      }

      num_read = fread(&label[i], sizeof(unsigned char), label_len, labelfile);
      if (label_len != num_read) {
        printf("Error: error reading file. 7\n");
        return 1;
      }

      // Remove end padding
      for (int j=0; j<label_len; ++j) {
        if (' ' == (int)label[i][j]) {
          label[i][j] = 0;
        }
      }
      label[i][label_len] == 0;

    } // for i

    // Permute labels to output each result with a uniform order of labels
    // By convention let line numbers bein increasing, e.g. "0 1" not "1 0"

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

    // Output line and bit numbers

    for (int i=0; i<num_way; ++i) {
      int iperm = perm[i];
      fprintf(metricstxtfile, 0 != i ? " " : "");
      fprintf(metricstxtfile, "%i %i", lineno[iperm], bitno[iperm]);
    } // for i

    // Output label

    for (int i=0; i<num_way; ++i) {
      int iperm = perm[i];
      fprintf(metricstxtfile, " %s_%c", label[iperm], al[iperm][bitno[iperm]]);
    } // for i

    // Output value

    fprintf(metricstxtfile, " %f\n", value);

  } // while

  fclose(metricsbinfile);
  fclose(alfile);
  fclose(labelfile);
  fclose(metricstxtfile);

  return 0;
}

//-----------------------------------------------------------------------------
