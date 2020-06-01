//-----------------------------------------------------------------------------
// Convert a SNP file in tped format to a binary file.
//-----------------------------------------------------------------------------

#include <stdio.h>

int main(int argc, char** argv) {

  if (argc != 3) {
    printf("preprocess: convert text SNP file to a packed binary SNP file.\n");
    printf("Usage: preprocess <snp_text_file> <snp_bin_file>\n");
    return 0;
  }

  FILE* snptxtfile = fopen(argv[1], "r");
  FILE* snpbinfile = fopen(argv[2], "wb");
  //FILE* snptxtfile = fopen("28M_permuted/28M.txt", "r");
  //FILE* snpbinfile = fopen("28M_permuted/28M.bin", "wb");

  if (!snptxtfile) {
    printf("Error: unable to open file. %s\n", argv[1]);
    return 1;
  }

  if (!snpbinfile) {
    printf("Error: unable to open file. %s\n", argv[2]);
    return 1;
  }

  const int max_line_len = 4096;
  int line[max_line_len];

  int c = 0;
  int cprev = 0;
  int line_len = 0;
  int clo = 0;
  int chi = 0;
  int col = 0;

  // Loop to input chars from input file

  while ((c = fgetc(snptxtfile)) != EOF) {

    line[line_len] = c;
    line_len++;
    if (line_len > max_line_len) {
      printf("Error: input line too long\n");
      return 1;
    }

    // Reached end of a line, so process it

    if (c != '\n') {
      continue;
    }

    line_len--;

    // First pass: loop over tokens to get the allele labels

    int i, num_tabs;
    for (i=0, num_tabs=0; i<line_len; ++i) {

      c = line[i];
      if (c == '\t') {
        // Found tab, go to next token
        num_tabs++;
        continue;
      }
      if (num_tabs < 4) {
        // Skip first four tokens - pick these up with another command
        continue;
      }

      // Get token number
      col = num_tabs - 4;

      if (col % 2 == 1) {
        if ((c=='0' && cprev!='0') ||
            (c!='0' && cprev=='0')) {
          printf("Error: token pair must be both zero or both nonzero\n");
          return 1;
        }
      }

      // Record values of the tokens encountered

      if (c == '0') {
      } else if (clo == 0) {
        clo = c;
      } else if (clo < c && chi == 0) {
        chi = c;
      } else if (clo > c && chi == 0) {
        chi = clo;
        clo = c;
      }

      cprev = c;
    } // for i - first pass loop

    //----------

    if (chi == 0) {
      chi = clo;
    }

    if (col%2 == 0) {
      printf("Error: line has invalid number of tokens\n");
      return 1;
    }

    int out_buf = 0;
    int num_buf = 0;

    // Second pass: loop to output results

    for (i=0, num_tabs=0; i<line_len; ++i) {

      c = line[i];
      if (c == '\t') {
        num_tabs++;
        continue;
      }
      if (num_tabs < 4) {
        continue;
      }

      col = num_tabs - 4;

      if (col % 2 == 0) {
        cprev = c;
        continue;
      }

      const int c0 = cprev;
      const int c1 = c;

      // Next map the token pair to a pair of bits needed by CoMet.

      // Note the alphabetically first allele is mapped to the zero bit

      // NOTE: order of lines matters below.
      int sn =
          c0 == '0' && c1 == '0' ? 2 * (1) + 1 * (0)  //  00  ->  1 0
        : c0 == clo && c1 == clo ? 2 * (0) + 1 * (0)  //  AA  ->  0 0
        : c0 == clo && c1 == chi ? 2 * (0) + 1 * (1)  //  AB  ->  0 1
        : c0 == chi && c1 == clo ? 2 * (0) + 1 * (1)  //  BA  ->  0 1
        : c0 == chi && c1 == chi ? 2 * (1) + 1 * (1)  //  BB  ->  1 1
        : -1;

      if (sn == -1) {
        printf("Error: unknown error encountered mapping tokens\n");
        return 1;
      }

      // Insert the 2 bits into the ouput buffer

      out_buf = (out_buf << 2) | sn;
      num_buf++;

      // Output buffer has 8 bits, so flush

      if (num_buf == 4) {
        fputc(out_buf, snpbinfile);
        out_buf = 0;
        num_buf = 0;
      }

    } // for i - second pass loop

    //----------

    // Final flush of buffer

    if (num_buf != 0) {
      out_buf = out_buf << (2 * (4 - num_buf));
      fputc(out_buf, snpbinfile);
      out_buf = 0;
      num_buf = 0;
    }

    // Process next line

    line_len = 0;
    clo = 0;
    chi = 0;

  } // while c

  // Finish

  fclose(snptxtfile);
  fclose(snpbinfile);
  return 0;
}

//-----------------------------------------------------------------------------
