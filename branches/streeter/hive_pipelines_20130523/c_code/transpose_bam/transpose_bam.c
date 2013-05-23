#include <stdio.h>
#include <unistd.h> /* for getopt */
#include <stdlib.h>
#include <string.h>
#include <sam.h>
#include <sam_header.h>
#include <ksort.h>


void usage();
int bam_index_build(const char *fn);

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
  int j;
  bam1_t *b;
} heap1_t;

static inline int pos_cmp(const heap1_t a, const heap1_t b) {
  if (a.j < 0) {
    if (b.j < 0) return 0;
    return 1;
  }
  return (a.b->core.pos > b.b->core.pos || (a.b->core.pos == b.b->core.pos && bam1_strand(a.b) > bam1_strand(b.b)));
}

KSORT_INIT(heap, heap1_t, pos_cmp);

int merge_regions(int num_infiles, samfile_t **in_bams, bam_index_t **in_indexes, bamFile out_bam, int num_regions, char **regions) {
  int *tid, *beg, *end;
  int i,j,num_regions_collapsed;
  bam_iter_t *bam_iters = (bam_iter_t*) malloc(sizeof(bam_iter_t) * num_infiles);
  heap1_t *heap = (heap1_t*) malloc(sizeof(heap1_t) * num_infiles);

  tid = (int*) malloc(sizeof(int) * num_regions);
  beg = (int*) malloc(sizeof(int) * num_regions);
  end = (int*) malloc(sizeof(int) * num_regions);

  for (i=0; i<num_regions; i++) {
    int insert_i, insert_tid, insert_beg, insert_end;
    bam_parse_region(in_bams[0]->header, regions[i], &insert_tid, &insert_beg, &insert_end);
    if (insert_tid < 0) {
      fprintf(stderr, "could not recognise region %s\n", regions[i]);
      return -1;
    }
    for (insert_i = i; insert_i > 0; insert_i--) {
      if (insert_tid > tid[insert_i-1])
        break;
      if (insert_tid == tid[insert_i-1] && insert_beg > beg[insert_i-1])
        break;
      tid[insert_i] = tid[insert_i-1];
      beg[insert_i] = beg[insert_i-1];
      end[insert_i] = end[insert_i-1];
    }
    tid[insert_i] = insert_tid;
    beg[insert_i] = insert_beg;
    end[insert_i] = insert_end;
  }
  while (1) {
    int old_num_regions = num_regions;
    for(i=num_regions-1; i>0; i--) {
      if (tid[i-1] == tid[i] && end[i-1] >= beg[i]-1) {
        int i_shift;
        if (end[i] > end[i-1])
          end[i-1] = end[i];
        for (i_shift = i+1; i_shift < num_regions; i_shift++) {
          tid[i_shift-1] = tid[i_shift];
          beg[i_shift-1] = beg[i_shift];
          end[i_shift-1] = end[i_shift];
        }
        num_regions --;
      }
    }
    if (old_num_regions == num_regions)
      break;
  }


  for (j=0; j<num_infiles; j++) {
    heap1_t *h = heap + j;
    h->b = bam_init1();
  }

  for (i=0; i<num_regions; i++) {
    for (j=0; j<num_infiles; j++) {
      heap1_t *h = heap + j;
      bam_iters[j] = bam_iter_query(in_indexes[j], tid[i], beg[i], end[i]);
      h->j = j;
      while (1) {
        if (bam_iter_read(in_bams[j]->x.bam, bam_iters[j], h->b) <= 0) {
          h->j = -1;
          break;
        }
        if (h->b->core.tid == tid[i] && h->b->core.pos >= beg[i])
          break;
      }
    }
    ks_heapmake(heap, num_infiles, heap);
    while (heap->j >= 0) {
      bam_write1(out_bam, heap->b);
      if (bam_iter_read(in_bams[heap->j]->x.bam, bam_iters[heap->j], heap->b) <= 0)
        heap->j = -1;
      ks_heapadjust(heap, 0, num_infiles, heap);
    }

    for (j=0; j<num_infiles; j++)
      bam_iter_destroy(bam_iters[j]);
  }

  for (j=0; j<num_infiles; j++) {
    heap1_t *h = heap + j;
    bam_destroy1(h->b);
  }

  free(tid);
  free(beg);
  free(end);
  free(bam_iters);
  free(heap);
  return 0;
}


int main(int argc, char *argv[])
{
  char *out_fname = NULL;
  char **in_fnames = NULL;
  char **regions = NULL;
  int num_regions = 0;
  int build_index = 0;
  int c, i;
  int num_infiles;
  int ret = 0;
  samfile_t **in_bams = NULL;
  bam_index_t **in_indexes = NULL;
  bamFile out_bam = NULL;

  regions = (char**) malloc(sizeof(char*) * argc);

  while((c = getopt(argc, argv, "o:r:i")) != -1)
    switch (c) {
      case 'o':
        out_fname = optarg;
        break;
      case 'r':
        regions[num_regions] = optarg;
        num_regions ++;
        break;
      case 'i':
        build_index = 1;
        break;
      case '?':
        usage();
        break;
    }

  if ( !out_fname || num_regions == 0 ) {
    usage();
  }

  if (optind < argc)
    num_infiles = argc-optind;
  else
    usage();

  
  in_fnames = argv + optind;
  in_bams = (samfile_t**) malloc(sizeof(samfile_t*) * num_infiles);
  in_indexes = (bam_index_t**) malloc(sizeof(bam_index_t*) * num_infiles);
  for (i=0; i<num_infiles; i++) {
    if ((in_bams[i] = samopen(in_fnames[i], "rb", NULL)) == 0) {
      printf("Failed to open file %s\n", in_fnames[i]);
      exit(-1);
    }
    if ((in_indexes[i] = bam_index_load(in_fnames[i])) == 0) {
      printf("Failed to load index for file %s\n", in_fnames[i]);
      exit(-1);
    }
  }

  if ((out_bam = bam_open(out_fname, "wb")) == 0) {
    printf("Failed to open file %s\n", out_fname);
    exit(-1);
  }

  merge_bam_headers(num_infiles, in_bams, out_bam);

  ret = merge_regions(num_infiles, in_bams, in_indexes, out_bam, num_regions, regions);


  bam_close(out_bam);
  for (i=0; i<num_infiles; i++) {
    samclose(in_bams[i]);
    bam_index_destroy(in_indexes[i]);
  }
  free(in_bams);
  free(in_indexes);
  free(regions);

  if (build_index && ret == 0)
    ret = bam_index_build(out_fname);

  return (ret == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


void usage()
{
  printf("usage: transpose_bam -o output_file -r region [-i] input_file1 [input_file2 ...]\n");
  printf("input bam files must be indexed\n");
  printf("all input bam files must be mapped to the same reference sequences\n");
  printf("flag -i is to index the output");
  exit(-1);
}
