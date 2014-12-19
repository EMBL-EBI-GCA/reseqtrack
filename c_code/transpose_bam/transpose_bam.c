#include <stdio.h>
#include <unistd.h> /* for getopt */
#include <stdlib.h>
#include <string.h>
#include <sam_header.h>
#include <htslib/sam.h>
#include <sam.h>
#include <htslib/ksort.h>


void usage();
//int bam_index_build(const char *fn);

void append_line(bam_hdr_t *out_header, char *line, size_t l_length) {
  out_header->text = (char*) realloc(out_header->text, sizeof(char*) * (out_header->l_text + l_length));
  memcpy((char*)out_header->text + out_header->l_text, line, l_length);
  out_header->l_text += l_length;
}

void uniquify_bam1_rg(bam1_t *b, int j) {
      uint8_t *rg_name;
      uint8_t *ptr;
      int ori_length = b->l_data;
      int extra_length = 2;
      int j_test = j;
      while (j_test /= 10) {
        extra_length++;
      }
      b->l_data += extra_length;
      if (b->m_data < b->l_data) {
              b->m_data = b->l_data;
              kroundup32(b->m_data);
              b->data = (uint8_t*)realloc(b->data, b->m_data);
      }
      rg_name = bam_aux_get(b, "RG") +1;
      for (ptr=b->data + ori_length -1; ptr >= rg_name; ptr--) {
        *(ptr+extra_length) = *ptr;
      }
      sprintf(rg_name, "%u", j);
      *(rg_name + extra_length -1) = '.';
}

void merge_bam_headers(int num_infiles, samfile_t **in_bams, bamFile out_bam, int uniquify_rg) {
  int i;
  char* ptr;
  char* header_line = NULL;
  size_t l_line;

  bam_hdr_t *out_header = bam_header_init();
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
        if (ptr[1]=='R' && ptr[2] == 'G') {
          if (uniquify_rg) {
            int extra_length = 2;
            int i_test = i;
            char *new_line, *new_line_ptr;
            while (i_test /= 10) {
              extra_length++;
            }
            new_line = (char*) malloc(sizeof(char) * (length + extra_length) );
            new_line_ptr = new_line;
            while (ptr <= end) {
              *new_line_ptr = *ptr;
              new_line_ptr ++;
              if (ptr[0] == ':' && ptr[-1] == 'D' && ptr[-2] == 'I') {
                sprintf(new_line_ptr, "%d.", i);
                new_line_ptr += extra_length;
              }
              ptr ++;
            }

            //sprintf(new_line, "@RG\tID:%d.", i);
            //memcpy(new_line+7+extra_length, ptr+7, length-7);
            //new_line[length + extra_length] = '\0';
            append_line(out_header, new_line, length+extra_length);
          }
          else
            append_line(out_header, ptr, length);
        }
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

int merge_regions(int num_infiles, samfile_t **in_bams, char **in_fnames, bamFile out_bam, int num_regions, char **regions, int uniquify_rg) {
  int *tid, *beg, *end;
  int i,j;
  bam_iter_t **bam_iters = (bam_iter_t**) malloc(sizeof(bam_iter_t*) * num_infiles);
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
    bam_index_t* in_index = bam_index_load(in_fnames[j]);
    if (in_index == 0) {
      printf("Failed to load index for file %s\n", in_fnames[j]);
      return -1;
    }
    h->b = bam_init1();

    bam_iters[j] = (bam_iter_t*) malloc(sizeof(bam_iter_t) * num_regions);
    for (i=0; i<num_regions; i++) {
      bam_iters[j][i] = bam_iter_query(in_index, tid[i], beg[i], end[i]);
    }
    bam_index_destroy(in_index);
  }

  for (i=0; i<num_regions; i++) {
    for (j=0; j<num_infiles; j++) {
      heap1_t *h = heap + j;
      h->j = j;
      while (1) {
        if (bam_iter_read(in_bams[j]->x.bam, bam_iters[j][i], h->b) <= 0) {
          h->j = -1;
          break;
        }
        if (h->b->core.tid == tid[i] && h->b->core.pos >= beg[i])
          break;
      }
    }
    ks_heapmake(heap, num_infiles, heap);
    if (uniquify_rg)
      while (heap->j >= 0) {
        uniquify_bam1_rg(heap->b, heap->j);
        bam_write1(out_bam, heap->b);
        if (bam_iter_read(in_bams[heap->j]->x.bam, bam_iters[heap->j][i], heap->b) <= 0)
          heap->j = -1;
        ks_heapadjust(heap, 0, num_infiles, heap);
      }
    else
      while (heap->j >= 0) {
        bam_write1(out_bam, heap->b);
        if (bam_iter_read(in_bams[heap->j]->x.bam, bam_iters[heap->j][i], heap->b) <= 0)
          heap->j = -1;
        ks_heapadjust(heap, 0, num_infiles, heap);
      }

    for (j=0; j<num_infiles; j++)
      bam_iter_destroy(bam_iters[j][i]);
  }

  for (j=0; j<num_infiles; j++) {
    heap1_t *h = heap + j;
    bam_destroy1(h->b);
    free(bam_iters[j]);
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
  int uniquify_rg = 0;
  int c, i;
  int num_infiles;
  int ret = 0;
  samfile_t **in_bams = NULL;
  bamFile out_bam = NULL;

  regions = (char**) malloc(sizeof(char*) * argc);

  while((c = getopt(argc, argv, "o:r:iu")) != -1)
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
      case 'u':
        uniquify_rg = 1;
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
  for (i=0; i<num_infiles; i++) {
    if ((in_bams[i] = samopen(in_fnames[i], "rb", NULL)) == 0) {
      printf("Failed to open file %s\n", in_fnames[i]);
      exit(-1);
    }
  }

  if ((out_bam = bam_open(out_fname, "wb")) == 0) {
    printf("Failed to open file %s\n", out_fname);
    exit(-1);
  }

  merge_bam_headers(num_infiles, in_bams, out_bam, uniquify_rg);

  ret = merge_regions(num_infiles, in_bams, in_fnames, out_bam, num_regions, regions, uniquify_rg);


  bam_close(out_bam);
  for (i=0; i<num_infiles; i++) {
    samclose(in_bams[i]);
  }
  free(in_bams);
  free(regions);

  if (build_index && ret == 0)
    ret = bam_index_build(out_fname);

  return (ret == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


void usage()
{
  printf("usage: transpose_bam -o output_file -r region [-i] [-u] input_file1 [input_file2 ...]\n");
  printf("input bam files must be indexed\n");
  printf("all input bam files must be mapped to the same reference sequences\n");
  printf("flag -i is to index the output\n");
  printf("flag -u changes the RG tag to ensure uniqueness\n");
  exit(-1);
}
