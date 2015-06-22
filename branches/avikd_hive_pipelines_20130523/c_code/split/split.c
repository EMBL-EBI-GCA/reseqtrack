#include <stdio.h>
#include <zlib.h>
#include <unistd.h> /* for getopt */
#include <stdlib.h>
#include <string.h>

#define BUFFSIZE 131072 /* 128 kb */
//#define BUFFSIZE 8192 /* 8 kb */
//#define BUFFSIZE 32768 /* 32 kb */
//#define BUFFSIZE 1048576 /* 1 Mb */

void usage();
gzFile open_file(const char* const prefix, const char* const suffix, int filenum, int label_length);
void safe_close(const gzFile file);
void safe_write(const gzFile outfh, const char * const out_start, int out_length);

static int buffsize = BUFFSIZE;

int main(int argc, char *argv[])
{
  char *infile;
  int max_lines;
  int label_length = 1;
  char *outfile_prefix = "output";
  char *outfile_suffix = "";

  int c;
  while((c = getopt(argc, argv, "n:a:p:s:b:")) != -1)
    switch (c) {
      case 'n':
        max_lines = atoi(optarg);
        break;
      case 'a':
        label_length = atoi(optarg);
        break;
      case 'p':
        outfile_prefix = optarg;
        break;
      case 's':
        outfile_suffix = optarg;
        break;
      case 'b':
        buffsize = atoi(optarg);
        if (buffsize <= 0) {
          printf("Invalid value for buffer size: %s\n", optarg);
          usage();
        }
        break;
      case '?':
        break;
    }

  if (optind < argc)
    infile = argv[optind++];
  else
    usage();

  if (max_lines <= 0) {
    printf("Need to know maximum number lines in output files\n");
    usage();
  }

  int filenum = 1;

  gzFile outfh;
  gzFile infh = gzopen(infile, "r");
  if (! infh) {
    printf("input file %s could not be opened\n", infile);
    exit(1);
  }
  gzbuffer(infh, buffsize);
  int file_open = 0;

  char buffer[buffsize];
  int line_count = 0;
  while (1) {
    int buffer_length = gzread(infh, buffer, buffsize);
    if (buffer_length < 0) {
      printf("error reading input file %s", infile);
      exit(1);
    }
    if (buffer_length == 0) {
      break; // end of input file
    }

    if (! file_open) {
      outfh = open_file(outfile_prefix, outfile_suffix, filenum, label_length);
      file_open = 1;
    }

    char * out_start = buffer;
    char * out_end = out_start;
    for (; buffer_length > 0; buffer_length--, out_end++) {
      if (*out_end == '\n')
        line_count ++;
      if (line_count == max_lines) {
        int out_length = out_end - out_start + 1;
        safe_write(outfh, out_start, out_length);
        line_count = 0;
        safe_close(outfh);
        filenum++;
        out_start = out_end + 1;
        if (buffer_length >1)
          outfh = open_file(outfile_prefix, outfile_suffix, filenum, label_length);
        else
          file_open = 0;
      }
    }
    if (out_start < out_end) {
      int out_length = out_end - out_start;
      safe_write(outfh, out_start, out_length);
    }
  }
  
  safe_close(infh);
  if (file_open)
    safe_close(outfh);

  return 0;
}


void usage()
{
  printf("usage: split -n line_count [-p outfile_prefix] [-s outfile_suffix] [-a label_length] [-b buffer_size] file\n");
  exit(1);
}

gzFile open_file(const char* const prefix, const char* const suffix, int filenum, int label_length)
{
  int min_length = 1;
  int max_filenum = 10;
  while (max_filenum <= filenum) {
    min_length ++;
    max_filenum *= 10;
  }
  if (min_length > label_length)
    label_length = min_length;

  int num_zeros = label_length - min_length;

  int name_length = strlen(prefix) + label_length + strlen(suffix) + 4;
  char outfile[name_length];

  strcpy(outfile, prefix);
  char * fpnt = outfile + strlen(outfile);

  int i;
  for (i=0; i<num_zeros; i++) {
    fpnt[0] = '0';
    fpnt ++;
  }
  sprintf(fpnt, "%d%s.gz", filenum, suffix);

  gzFile outfh = gzopen(outfile, "w");
  if (! outfh) {
    printf("output file %s could not be opened\n", outfile);
    exit(1);
  }
  gzbuffer(outfh, buffsize);

  printf("now writing output to %s\n", outfile);

  return outfh;
}


void safe_close(const gzFile file) {
  int success = gzclose(file);
  if (success != Z_OK) {
    printf("Error closing file\n");
    exit(1);
  }
  return;
}

void safe_write(const gzFile outfh, const char * const out_start, int out_length) {
  int written_length = gzwrite(outfh, out_start, out_length);
  if (written_length != out_length) {
    printf("error writing output\n");
    exit(1);
  }
  return;
}
