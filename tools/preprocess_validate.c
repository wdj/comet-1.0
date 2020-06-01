//-----------------------------------------------------------------------------
// Validate the result of running the preprocess command.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

  if (argc != 4) {
    printf("preprocess_validate: check binary SNP file against text SNP file.\n");
    printf("Usage: preprocess_validate <snp_bin_file> <snp_text_file> <allele_label_file>\n");
    return 0;
  }

  FILE* snpbinfile = fopen(argv[1], "rb");
  FILE* snptxtfile = fopen(argv[2], "r");
  FILE* alfile = fopen(argv[3], "r");
  //FILE* snpbinfile = fopen("28M.bin", "rb");
  //FILE* snptxtfile = fopen("28M.txt", "r");
  //FILE* alfile = fopen("28M_allele_labels.txt", "r");

  if (!snpbinfile) {
    printf("Error: unable to open file. %s\n", argv[1]);
    return 1;
  }

  if (!snptxtfile) {
    printf("Error: unable to open file. %s\n", argv[2]);
    return 1;
  }

  if (!alfile) {
    printf("Error: unable to open file. %s\n", argv[3]);
    return 1;
  }

  //const int num_vector = 28342758;
  //const int num_vector = 10;
  //const int num_field = 882;

  size_t num_checked = 0;
  size_t num_validated = 0;

  // Loop over lines

  for (int line_num=0; ; ++line_num) {

    // Read from allele label file to get info for this line

    int c = 0;
    if ((c = fgetc(alfile)) == EOF) {
      break; // Check complete
    }

    const int clo = c;
    const int chi = fgetc(alfile);
    int discard = fgetc(alfile); // skip newline

    // Skip first four tokens of snptxtfile

    int num_tabs = 0;
    while ((c = fgetc(snptxtfile)) != EOF) {
      if (c == '\t') {
        num_tabs++;
      }
      if (num_tabs == 4) {
        break;
      }
    } // while

    // Loop over elements of this line

    int inbyte = 0;
    for (int elt_num=0; ; ++elt_num) {

      const int c0true = fgetc(snptxtfile);
      discard = fgetc(snptxtfile); // skip tab
      const int c1true = fgetc(snptxtfile);
      discard = fgetc(snptxtfile); // skip tab
      const int is_end_of_line = discard == 10;

      // Extract 2-bit seminibble from bin file

      const int sn_num = elt_num % 4;

      if (sn_num == 0) {
        inbyte = fgetc(snpbinfile); // Get new byte
      }

      const int sn = (inbyte >> (6 - 2*sn_num) ) & 3;

      // Map

      int c0 = 0;
      int c1 = 0;

      if (sn == 0) {
        c0 = clo;
        c1 = clo;
      } else if (sn == 1) {
        c0 = clo;
        c1 = chi;
      } else if (sn == 2) {
        c0 = '0';
        c1 = '0';
      } else if (sn == 3) {
        c0 = chi;
        c1 = chi;
      }

      // Check

      //printf("%c %c %c %c\n", c0, c1, c0true, c1true);
      if ((c0==c0true && c1==c1true) || (c0==c1true && c1==c0true)) {
        num_validated++;
      } else {
        printf("Error: invalid value detected, line %i\n", line_num);
        return 1;
      }
      num_checked++;

      if (is_end_of_line) {
        break;
      }

    } // elt_num

  } // line_num

  printf("Number of elements validated: %ul of %ul\n",
         num_validated, num_checked);

  // Finish

  fclose(alfile);
  fclose(snpbinfile);
  fclose(snptxtfile);

  return 0;
}

//-----------------------------------------------------------------------------
