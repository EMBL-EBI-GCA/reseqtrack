#include <stdio.h>
#include <unistd.h> /* for getopt */
#include <stdlib.h>
#include <string.h>
#include <sam.h>
#include <sam_header.h>
#include <ksort.h>


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

typedef struct {
  int i;
  bam1_t *b;
} heap1_t;

static inline int pos_cmp(const heap1_t a, const heap1_t b) {
  if (a.i < 0) {
    if (b.i < 0) return 0;
    return 1;
  }
  return (a.b->core.pos > b.b->core.pos || (a.b->core.pos == b.b->core.pos && bam1_strand(a.b) > bam1_strand(b.b)));
}

KSORT_INIT(heap, heap1_t, pos_cmp);

int merge_region(int num_infiles, samfile_t **in_bams, char **in_fnames, bamFile out_bam, char *region) {
  int i;
  bam_iter_t *bam_iters = (bam_iter_t*) malloc(sizeof(bam_iter_t) * num_infiles);
  heap1_t *heap = (heap1_t*) malloc(sizeof(heap1_t) * num_infiles);

  for (i=0; i<num_infiles; i++) {
    int tid, beg, end;
    bam_index_t *idx;
    heap1_t *h = heap + i;
    bam_parse_region(in_bams[0]->header, region, &tid, &beg, &end);
    if (tid < 0) {
      fprintf(stderr, "could not recognise region %s\n", region);
      return -1;
    }
    if ((idx = bam_index_load(in_fnames[i])) == 0) {
      printf("Failed to load index for file %s\n", in_fnames[i]);
      return -1;
    }
    bam_iters[i] = bam_iter_query(idx, tid, beg, end);
    bam_index_destroy(idx);
    h->b = bam_init1();
    h->i = i;
    if (bam_iter_read(in_bams[i]->x.bam, bam_iters[i], h->b) <= 0)
      h->i = -1;
  }

  ks_heapmake(heap, num_infiles, heap);

  while (heap->i >= 0) {
    bam_write1(out_bam, heap->b);
    if (bam_iter_read(in_bams[heap->i]->x.bam, bam_iters[heap->i], heap->b) <= 0) {
      heap->i = -1;
    }
    ks_heapadjust(heap, 0, num_infiles, heap);
  }

  for (i=0; i<num_infiles; i++) {
    heap1_t *h = heap + i;
    bam_iter_destroy(bam_iters[i]);
    bam_destroy1(h->b);
  }
  free(bam_iters);
  free(heap);
  return 0;
}


int main(int argc, char *argv[])
{
  char *out_fname = NULL;
  char **in_fnames = NULL;
  char *region = NULL;
  int c, i;
  int num_infiles;
  int ret = 0;
  samfile_t **in_bams = NULL;
  bamFile out_bam = NULL;

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

  
  in_fnames = argv + optind;
  in_bams = (samfile_t**) malloc(sizeof(samfile_t*) * num_infiles);
  for (i=0; i<num_infiles; i++)
    if ((in_bams[i] = samopen(in_fnames[i], "rb", NULL)) == 0) {
      printf("Failed to open file %s\n", in_fnames[i]);
      exit(-1);
    }

  if ((out_bam = bam_open(out_fname, "wb")) == 0) {
    printf("Failed to open file %s\n", out_fname);
    exit(-1);
  }

  merge_bam_headers(num_infiles, in_bams, out_bam);

  ret = merge_region(num_infiles, in_bams, in_fnames, out_bam, region);


  bam_close(out_bam);
  for (i=0; i<num_infiles; i++) {
    samclose(in_bams[i]);
  }
  free(in_bams);

  return (ret == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


void usage()
{
  printf("usage: transpose_bam -o output_file -r region input_file1 [input_file2 ...]\n");
  printf("bam files must be indexed\n");
  printf("all bam files must be mapped to the same reference sequences\n");
  exit(-1);
}
