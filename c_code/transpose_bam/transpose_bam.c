#include <stdio.h>
#include <unistd.h> /* for getopt */
#include <stdlib.h>
#include <string.h>
#include <sam.h>
#include <sam_header.h>


void usage();

void append_line(bam_header_t *out_header, char *line, size_t l_length) {
  out_header->text = (char*) realloc(out_header->text, sizeof(char*) * (out_header->l_text + l_length));
  memcpy((char*)out_header->text + out_header->l_text, line, l_length);
  out_header->l_text += l_length;
}

void merge_bam_headers(int num_infiles, samfile_t **in_bams, bamFile out_bam) {
  int i;
  char* ptr;
  char* header_line = NULL;
  size_t l_line;

  bam_header_t *out_header = bam_header_init();
  for (i=0; i<num_infiles; i++) {
    char* ptr = in_bams[i]->header->text;
    int read = 0;
    while (read < in_bams[i]->header->l_text) {
      size_t length;
      char *end;
      if (ptr[0] != '@') {
        ptr++;
        read++;
        continue;
      }
      end = strchr(ptr, '\n');
      length = end - ptr + 1;
      if (length > 4) {
        if (ptr[1]=='R' && ptr[2] == 'G')
          append_line(out_header, ptr, length);
        else if (i==0)
          if ((ptr[1] == 'H' && ptr[2] == 'D') || (ptr[1] == 'S' && ptr[2] == 'Q'))
            append_line(out_header, ptr, length);
      }
      read += length;
      ptr = end + 1;
    }
  }

  out_header->n_targets = in_bams[0]->header->n_targets;
  out_header->target_name = in_bams[0]->header->target_name;
  out_header->target_len = in_bams[0]->header->target_len;
  bam_header_write(out_bam, out_header);
  out_header->n_targets = 0;
  out_header->target_name = NULL;
  out_header->target_len = NULL;

  bam_header_destroy(out_header);

  return;
}


int merge_region(int num_infiles, samfile_t **in_bams, bam_index_t **bam_indices, bamFile out_bam, char *region) {
  int i;
  bam_iter_t *bam_iters = (bam_iter_t*) malloc(sizeof(bam_iter_t) * num_infiles);
  bam1_t **b = (bam1_t**) malloc(sizeof(bam1_t*) * num_infiles);
  int *rets = (int*) malloc(sizeof(int) * num_infiles);
  for (i=0; i<num_infiles; i++) {
    int tid, beg, end;
    bam_parse_region(in_bams[0]->header, region, &tid, &beg, &end);
    if (tid < 0) {
      fprintf(stderr, "could not recognise region %s\n", region);
      return -1;
    }
    bam_iters[i] = bam_iter_query(bam_indices[i], tid, beg, end);
    b[i] = bam_init1();
    rets[i] = bam_iter_read(in_bams[i]->x.bam, bam_iters[i], b[i]);
  }

  while (1) {
    int i_write = -1;
    for (i=0; i<num_infiles; i++) {
      if (rets[i] <= 0) continue;
      if (i_write < 0) {
        i_write = i;
        continue;
      }
      if (b[i]->core.pos < b[i_write]->core.pos)
        i_write = i;
    }
    if (i_write < 0) break;
    bam_write1(out_bam, b[i_write]);
    rets[i_write] = bam_iter_read(in_bams[i_write]->x.bam, bam_iters[i_write], b[i_write]);

  }
  

  for (i=0; i<num_infiles; i++) {
    bam_iter_destroy(bam_iters[i]);
    bam_destroy1(b[i]);
  }
  free(bam_iters);
  free(b);
  return 0;
}


int main(int argc, char *argv[])
{
  char *out_fname = NULL;
  char *region = NULL;
  int c, i;
  int num_infiles;
  int ret = 0;
  samfile_t **in_bams = NULL;
  bamFile out_bam = NULL;
  bam_index_t **bam_indices = NULL;

  while((c = getopt(argc, argv, "o:r:")) != -1)
    switch (c) {
      case 'o':
        out_fname = optarg;
        break;
      case 'r':
        region = optarg;
        break;
      case '?':
        usage();
        break;
    }

  if ( !out_fname || !region ) {
    usage();
  }

  if (optind < argc)
    num_infiles = argc-optind;
  else
    usage();

  
  in_bams = (samfile_t**) malloc(sizeof(samfile_t*) * num_infiles);
  bam_indices = (bam_index_t**) malloc(sizeof(bam_index_t*) * num_infiles);
  for (i=0; i<num_infiles; i++) {
    if ((in_bams[i] = samopen(argv[optind], "rb", NULL)) == 0) {
      printf("Failed to open file %s\n", argv[optind]);
      exit(-1);
    }
    if ((bam_indices[i] = bam_index_load(argv[optind])) == 0) {
      printf("Failed to load index for file %s\n", argv[optind]);
      exit(-1);
    }
    optind++;
  }

  if ((out_bam = bam_open(out_fname, "wb")) == 0) {
    printf("Failed to open file %s\n", out_fname);
    exit(-1);
  }

  merge_bam_headers(num_infiles, in_bams, out_bam);

  ret = merge_region(num_infiles, in_bams, bam_indices, out_bam, region);


  bam_close(out_bam);
  for (i=0; i<num_infiles; i++) {
    samclose(in_bams[i]);
    bam_index_destroy(bam_indices[i]);
  }

  free(in_bams);
  free(bam_indices);

  return (ret == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


void usage()
{
  printf("usage: transpose_bam -o output_file -r region input_file1 [input_file2 ...]\n");
  exit(-1);
}
