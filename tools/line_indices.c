//-----------------------------------------------------------------------------
// Calculate the byte offset for each line of the SNP text file; store as
// 8-byte unsigned int.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

  if (argc != 3) {
    printf("line_indices: get byte offsets of lines of the SNP text file.\n");
    printf("Usage: line_indices <snp_text_file> <line_indices_file>\n");
    return 0;
  }

  FILE* snptxtfile = fopen(argv[1], "r");
  FILE* lifile = fopen(argv[2], "wb");
  //FILE* snptxtfile = fopen("28M_permuted/28M.txt", "r");
  //FILE* lifile = fopen("28M_permuted/28M_line_indices.txt", "wb");

  if (!snptxtfile) {
    printf("Error: unable to open file. %s\n", argv[1]);
    return 1;
  }

  if (!lifile) {
    printf("Error: unable to open file. %s\n", argv[1]);
    return 1;
  }

  int c = 0;
  int cprev = '\n';
  size_t line_index = 0;

  while ((c = fgetc(snptxtfile)) != EOF) {
    if (cprev == '\n') {
      size_t num_written = fwrite(&line_index, sizeof(size_t), 1, lifile);
      if (1 != num_written) {
        printf("Error: failure to write result");
        return 1;
      }
    }
    line_index++;
    cprev = c;
  } // while

  fclose(snptxtfile);
  fclose(lifile);
  return 0;
}

//-----------------------------------------------------------------------------
