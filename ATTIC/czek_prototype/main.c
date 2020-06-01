/****************************************************************/
/* Czekanowski Similarity Metric                                */
/*                                                              */
/* Code Author: Doug Hyatt                                      */
/* Created: June, 2015                                          */
/****************************************************************/

#include <string.h>

#include "mpi.h"

#include "czek.h"

#define VERSION "0.0.2"
#define DATE "June 7, 2015"
#define TEXTSIZE 10000

struct _option
{
  char letter;                 /* Single letter identifier for option */
  char optarg[TEXTSIZE];       /* Parameter (if option takes one) */
  int index;                   /* Index of next argument to parse */
};

void version();
void usage(char *);
void help();
void parse_arguments(int, char **, char *, char *);
void get_option(int, char **, struct _option *);

/*===========================================================================*/
/* Main Routine */

int main(int argc, char **argv) {
  int num_vectors = 0;       /* Number of vectors read */
  int num_fields = 0;        /* Number of fields per line */
  char input_file[TEXTSIZE] = "";
  char output_file[TEXTSIZE] = "";
  char err_msg[TEXTSIZE] = "";
  FILE *input_ptr = stdin;   /* Input file pointer */
  FILE *output_ptr = stdout; /* Output file pointer */
  struct _vector *vectors = 0;   /* Vectors with ids/data */
  int comm_rank = 0;
  int comm_size = 1;
  int vec_ctr = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  /* parse arguments */
  parse_arguments(argc, argv, input_file, output_file);

  /*---Open files, read vectors on proc 0---*/

  if (comm_rank == 0) {
    /* open files */
    if (input_file[0] != '\0')
    {
      input_ptr = fopen(input_file, "r");
      if (input_ptr == NULL)
      {
        sprintf(err_msg, "\nError: can't open input file %s", input_file);
        perror(err_msg);
        exit(2);
      }
    }
    if (output_file[0] != '\0')
    {
      output_ptr = fopen(output_file, "w");
      if (output_ptr == NULL)
      {
        sprintf(err_msg, "\nError: can't open output file %s", output_file);
        perror(err_msg);
        exit(2);
      }
    }

    /* read in the vectors */
    vectors = read_vectors(input_ptr, &num_vectors, &num_fields);
    if (vectors == NULL) {
      fprintf(stderr, "Error reading in the vectors.\n");
      exit(3);
    }
  }

  /*---Broadcast sizes to other procs---*/

  MPI_Bcast(&num_vectors, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&num_fields, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /*---Allocate vecs on other procs---*/

  if (comm_rank != 0) {
    vectors = (struct _vector*)malloc(num_vectors * sizeof(struct _vector));
    for (vec_ctr=0; vec_ctr<num_vectors; ++vec_ctr) {
      vectors[vec_ctr].id[0] = '\0';
      vectors[vec_ctr].data = (Float_t*)malloc(num_fields * sizeof(Float_t));
    }
  }

  /*---Broadcast vecs to other procs---*/

  for (vec_ctr=0; vec_ctr<num_vectors; ++vec_ctr) {
    MPI_Bcast(vectors[vec_ctr].data, num_fields,
              MPI_Float_t, 0, MPI_COMM_WORLD);
  }

  if (0) {
    int i = 0;
    int j = 0;
    for (i = 0; i < num_vectors; i++) {
      printf("%s", vectors[i].id);
      for (j = 0; j < num_fields; j++) printf("\t%.2f", vectors[i].data[j]);
      printf("\n");
    }
  }

  /* do all the czekanowki calculations and output them */
#if 0
  process_vectors(vectors, num_vectors, num_fields, output_ptr);
  process_vectors_alt(vectors, num_vectors, num_fields, output_ptr);
  process_vectors_alt2(vectors, num_vectors, num_fields, output_ptr);
  process_vectors_alt3(vectors, num_vectors, num_fields, output_ptr);
  process_vectors_alt4(vectors, num_vectors, num_fields, output_ptr);
  process_vectors_alt5(vectors, num_vectors, num_fields, output_ptr);
  process_vectors_alt6(vectors, num_vectors, num_fields, output_ptr);
#endif
  PROCESS_VECTORS_FN (vectors, num_vectors, num_fields, output_ptr);

  /*---Close files---*/

  if (comm_rank == 0) {
    /* close filehandles */
    if (input_ptr != stdin) fclose(input_ptr); 
    if (output_ptr != stdout) fclose(output_ptr); 
  }

  /* free vectors */
  free_vectors(vectors, num_vectors);

  MPI_Finalize();
  exit(0);
}

/*===========================================================================*/

void parse_arguments(int argc, char **argv, char *input, char *output) {
/*  char err_msg[TEXTSIZE] = "";
  char digits[10] = "0123456789"; */
  struct _option option = {0};

  option.index = 1;

  /* Parse command line options */
  while(1)
  {
    if (argc == 1)
    {
      break;
    }
    get_option(argc, argv, &option);
    switch (option.letter)
    {
      case 'h':
        help();
      case 'i':
        strcpy(input, option.optarg);
        break;
      case 'o':
        strcpy(output, option.optarg);
        break;
      case 'v':
        version();
      default:
        break;
    }
    if (option.index == -1) break;
  }
}

/*===========================================================================*/
/* Get next option and parse out argument, update index */

void get_option(int argc, char **argv, struct _option *opt)
{
  int opt_len = 0;
  int requires_arg = 0;           /* Set to 1 if option requires an arg */
  int arg_included = 0;           /* Set to 1 if user included arg in option */
  char *parse = NULL;
  char *tmp_ptr = NULL;
  char tmp_string[TEXTSIZE] = "";

  opt->optarg[0] = '\0';
  parse = argv[opt->index];
  /* Not an option */
  if (parse == NULL || parse[0] != '-' || strlen(parse) <= 1)
  {
    sprintf(tmp_string, "Unrecognized option '%s'.", parse);
    usage(tmp_string);
  }
  /* Long Options */
  if (strlen(parse) > 2 && strncmp(parse, "--", 2) == 0)
  {
    tmp_ptr = strchr(parse, '=');
    opt_len = strlen(parse);
    if (tmp_ptr != NULL)
    {
      opt_len -= strlen(tmp_ptr);
    }
    strncpy(tmp_string, parse, opt_len);
    tmp_string[opt_len] = '\0';
    if (strcmp(tmp_string, "--help") == 0)
    {
      opt->letter = 'h';
      requires_arg = 0;
    }
    else if (strcmp(tmp_string, "--input") == 0)
    {
      opt->letter = 'i';
      requires_arg = 1;
    }
    else if (strcmp(tmp_string, "--output") == 0)
    {
      opt->letter = 'o';
      requires_arg = 1;
    }
    else if (strcmp(tmp_string, "--version") == 0)
    {
      opt->letter = 'v';
      requires_arg = 0;
    }
    else
    {
      sprintf(tmp_string, "Unrecognized option '%s'.", parse);
      usage(tmp_string);
    }
    if (tmp_ptr != NULL)
    {
      arg_included = 1;
      strcpy(opt->optarg, tmp_ptr+1);
    }
  }
  /* Short option */
  if (parse[0] == '-' && parse[1] != '-')
  {
    switch (parse[1])
    {
      case 'i':
      case 'o':
        requires_arg = 1;
        opt->letter = parse[1];
        if (strlen(parse) > 2)
        {
          arg_included = 1;
          strcpy(opt->optarg, parse+2);
        }
        break;
      case 'h':
      case 'v':
        requires_arg = 0;
        opt->letter = parse[1];
        if (strlen(parse) > 2)
        {
          arg_included = 1;
          strcpy(opt->optarg, parse+2);
        }
        break;
      default:
        sprintf(tmp_string, "Unrecognized option '%s'.", parse);
        usage(tmp_string);
        break;
    }
  }
  /* Option takes no argument but one was given */
  if (requires_arg == 0 && (arg_included == 1 ||
      (opt->index < argc-1 && argv[opt->index+1][0] != '-')))
  {
    sprintf(tmp_string, "Option -%c should not have an argument.",
            opt->letter);
    usage(tmp_string);
  }
  /* Option requires an argument but none was given */
  if (requires_arg == 1 && (arg_included == 0 && (opt->index >= argc-1 ||
      argv[opt->index+1][0] == '-')))
  {
    sprintf(tmp_string, "Option -%c requires an argument.",
            opt->letter);
    usage(tmp_string);
  }
  /* Option requires an argument and needs opt->optarg to be copied */
  if (requires_arg == 1 && arg_included == 0)
  {
    strcpy(opt->optarg, argv[opt->index+1]);
  }
  else if (requires_arg == 1 && arg_included == 1)
  {
    strcpy(opt->optarg, opt->optarg);
  }
  /* Update index, -1 if no more arguments to parse */
  opt->index += (1 + requires_arg - arg_included);
  if (opt->index >= argc)
  {
    opt->index = -1;
  }
}

/*===========================================================================*/

void usage(char *msg) {
  fprintf(stderr, "\nError: %s\n", msg);
  fprintf(stderr, "\nUsage:  czek [-i input_file] [-h]");
  fprintf(stderr, " [-o output_file] [-v]\n");
  exit(1);
}

/*===========================================================================*/

void help() {
  printf("\nUsage:  czek [-i input_file] [-h]");
  printf(" [-o output_file] [-v]\n\n");
  printf("     -h, --help:     Print help information and exit.\n");
  printf("     -i, --input:    Specify input file (default stdin).\n");
  printf("     -o, --output:   Specify output file (default stdout).\n");
  printf("     -v, --version:  Print version information and exit.\n\n");
  exit(0);
}

/*===========================================================================*/

/* Print version number and exit */
void version()
{
  printf("\nczekanowski v%s: %s\n\n", VERSION, DATE);
  exit(0);
}

/*===========================================================================*/
