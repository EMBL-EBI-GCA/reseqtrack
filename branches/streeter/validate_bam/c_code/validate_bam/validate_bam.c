#include <stdio.h>
#include <unistd.h> /* for getopt */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sam.h>
#include <khash.h>
#include <openssl/md5.h> /* for md5 */
#include <sys/mman.h> /* for md5 */
#include <fcntl.h> /* for md5 */


void usage();

KHASH_MAP_INIT_INT(int32, uint32_t); // used to store insert sizes
KHASH_MAP_INIT_STR(tags, char*); // used to store bam header tags

// struct containing statistics for a read group
typedef struct
{
  char *library; // points to a bam_header string
  char *sample; // points to a bam_header string
  char *platform; // points to a bam_header string
  char *description; // points to a bam_header string
  uint32_t num_total_bases;
  uint32_t num_mapped_bases;
  uint32_t num_total_reads;
  uint32_t num_mapped_reads;
  uint32_t num_mapped_reads_paired_in_sequencing;
  uint32_t num_mapped_reads_properly_paired;
  uint32_t num_NM_bases;
  uint32_t num_NM_mismatched_bases;
  uint32_t sum_quality_mapped_bases;
  uint32_t sum_insert_sizes;
  uint32_t num_insert_sizes;
  uint32_t num_duplicate_reads;
  uint32_t num_duplicate_bases;

  khash_t(int32) *insert_sizes; // key is the insert size, value is the number of occurrences

  uint32_t median_insert_size;
}
rg_stats_t;

KHASH_MAP_INIT_STR(rg_stats, rg_stats_t*);

// Initialise a khash containing: key is read group ID, value is an initialised rg_stats_t
void *rg_stats_init(bam_header_t *bam_header) {
  khash_t(rg_stats) *rg_stats_hash = kh_init(rg_stats);

  bam_header->dict = sam_header_parse2(bam_header->text);
  int n_rg;
  char **ID_array = sam_header2list(bam_header->dict, "RG", "ID", &n_rg);
  khash_t(tags) *LB_hash = sam_header2tbl(bam_header->dict, "RG", "ID", "LB");
  khash_t(tags) *SM_hash = sam_header2tbl(bam_header->dict, "RG", "ID", "SM");
  khash_t(tags) *PL_hash = sam_header2tbl(bam_header->dict, "RG", "ID", "PL");
  khash_t(tags) *DS_hash = sam_header2tbl(bam_header->dict, "RG", "ID", "DS");

  int i;
  for(i=0; i<n_rg; i++) {
    rg_stats_t *new_rg_stats = malloc(sizeof(rg_stats_t));

    int ret;
    khiter_t k = kh_put(rg_stats, rg_stats_hash, ID_array[i], &ret);
    kh_value(rg_stats_hash, k) = new_rg_stats;

    k = kh_get(tags, LB_hash, ID_array[i]);
    new_rg_stats->library = (k == kh_end(LB_hash)) ? NULL : kh_value(LB_hash, k);

    k = kh_get(tags, SM_hash, ID_array[i]);
    new_rg_stats->sample = (k == kh_end(SM_hash)) ? NULL : kh_value(SM_hash, k);

    k = kh_get(tags, PL_hash, ID_array[i]);
    new_rg_stats->platform = (k == kh_end(PL_hash)) ? NULL : kh_value(PL_hash, k);

    k = kh_get(tags, DS_hash, ID_array[i]);
    new_rg_stats->description = (k == kh_end(DS_hash)) ? NULL : kh_value(DS_hash, k);

    new_rg_stats->insert_sizes = kh_init(int32);
    new_rg_stats->num_total_bases = 0;
    new_rg_stats->num_mapped_bases = 0;
    new_rg_stats->num_total_reads = 0;
    new_rg_stats->num_mapped_reads = 0;
    new_rg_stats->num_mapped_reads_paired_in_sequencing = 0;
    new_rg_stats->num_mapped_reads_properly_paired = 0;
    new_rg_stats->num_NM_bases = 0;
    new_rg_stats->num_NM_mismatched_bases = 0;
    new_rg_stats->sum_quality_mapped_bases = 0;
    new_rg_stats->sum_insert_sizes = 0;
    new_rg_stats->num_insert_sizes = 0;
    new_rg_stats->num_duplicate_reads = 0;
    new_rg_stats->num_duplicate_bases = 0;
    new_rg_stats->median_insert_size = 0;
  }

  free(ID_array);
  kh_destroy(tags, LB_hash);
  kh_destroy(tags, SM_hash);
  kh_destroy(tags, PL_hash);
  kh_destroy(tags, DS_hash);
  return rg_stats_hash;
}

// destroy the khash of rg_stats.  Also destroy each rg_stats.
void rg_stats_destroy(void *_rg_stats_hash) {
  khash_t(rg_stats) * rg_stats_hash = _rg_stats_hash;
  khiter_t k;
  for (k = kh_begin(rg_stats_hash); k != kh_end(rg_stats_hash); k++)
    if (kh_exist(rg_stats_hash, k)) {
      rg_stats_t *rg_stats = kh_value(rg_stats_hash, k);
      kh_destroy(int32, rg_stats->insert_sizes);
    }
  kh_destroy(rg_stats, rg_stats_hash);
  return;
}

// update the stats for a single alignment
void update_stats(rg_stats_t *rg_stats, bam1_t *bam_line) {
  rg_stats->num_total_reads ++;
  rg_stats->num_total_bases += bam_line->core.l_qseq;


  if (! (bam_line->core.flag & BAM_FUNMAP)) {
    rg_stats->num_mapped_reads ++;

    uint32_t *cigar = bam1_cigar(bam_line);
    uint8_t *qual = bam1_qual(bam_line);
    uint8_t i,j;
    for (i=0; i < bam_line->core.n_cigar; i++){
      int op = bam_cigar_op(cigar[i]);
      int oplen = bam_cigar_oplen(cigar[i]);
      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        rg_stats->num_mapped_bases += oplen;
        for (j=0; j<oplen; j++)
          rg_stats->sum_quality_mapped_bases += qual[j];
      }
      if (op != BAM_CDEL && op != BAM_CREF_SKIP)
        qual += oplen;
    }

    if (bam_line->core.flag & BAM_FPAIRED)
      rg_stats->num_mapped_reads_paired_in_sequencing ++;

    if (bam_line->core.flag & BAM_FPROPER_PAIR) {
      rg_stats->num_mapped_reads_properly_paired ++;
      if (bam_line->core.qual > 0 && bam_line->core.isize > 0) {
        rg_stats->num_insert_sizes ++;
        rg_stats->sum_insert_sizes += bam_line->core.isize;

        int ret;
        khiter_t k = kh_put(int32, rg_stats->insert_sizes, bam_line->core.isize, &ret);
        if (!ret)
          kh_value(rg_stats->insert_sizes, k) ++;
        else
          kh_value(rg_stats->insert_sizes, k) = 1;
      }
    }

/*
    uint8_t *qual = bam1_qual(bam_line);
    for (i=0; i < bam_line->core.l_qseq; i++)
      rg_stats->sum_quality_mapped_bases += qual[i];
      */
  }

  uint8_t *nm_with_type = bam_aux_get(bam_line, "NM");
  if (nm_with_type != 0) {
    rg_stats->num_NM_mismatched_bases += bam_aux2i(nm_with_type);
    rg_stats->num_NM_bases += bam_line->core.l_qseq;
  }


  if (bam_line->core.flag & BAM_FDUP) {
    rg_stats->num_duplicate_reads ++;
    rg_stats->num_duplicate_bases += bam_line->core.l_qseq;
  }
}

float calc_percentage_mismatched_bases(rg_stats_t *rg_stats) {
  return 100 * (float) rg_stats->num_NM_mismatched_bases / rg_stats->num_NM_bases;
}

float calc_average_quality_mapped_bases(rg_stats_t *rg_stats) {
  return (float) rg_stats->sum_quality_mapped_bases / rg_stats->num_total_bases;
}

double calc_mean_insert_size(rg_stats_t *rg_stats) {
  return (double) (rg_stats->sum_insert_sizes) / rg_stats->num_insert_sizes;
}

double calc_insert_size_sd(rg_stats_t *rg_stats) {
  double mean = calc_mean_insert_size(rg_stats);
  double sd;
  khiter_t k;
  for (k = kh_begin(rg_stats->insert_sizes); k != kh_end(rg_stats->insert_sizes); k++)
    if (kh_exist(rg_stats->insert_sizes, k)) {
      double diff = kh_key(rg_stats->insert_sizes, k) - mean;
      sd += diff * diff * kh_value(rg_stats->insert_sizes, k);
    }
  sd /= rg_stats->num_insert_sizes;
  sd = sqrt(sd);
  return sd;
}

// needed for the qsort function
int compare_integers (const void * a, const void *b) {
  return ( *(uint32_t*)a - *(uint32_t*)b );
}
uint32_t calc_median_insert_size(rg_stats_t *rg_stats) {
  uint32_t *insert_sizes = malloc(sizeof(uint32_t) * kh_size(rg_stats->insert_sizes));

  khiter_t k;
  uint8_t arr_size=0;
  for (k = kh_begin(rg_stats->insert_sizes); k != kh_end(rg_stats->insert_sizes); k++)
    if (kh_exist(rg_stats->insert_sizes, k)) {
      insert_sizes[arr_size] = kh_key(rg_stats->insert_sizes, k);
      arr_size++;
    }

  qsort(insert_sizes, arr_size, sizeof(uint32_t), compare_integers);

  uint32_t median;
  uint32_t num_inserts_counted = 0;
  uint32_t stop_counting = (uint32_t) (0.5 * rg_stats->num_insert_sizes);
  uint8_t i;
  for (i=0; num_inserts_counted < stop_counting ; i++) {
    k = kh_get(int32, rg_stats->insert_sizes, insert_sizes[i]);
    num_inserts_counted += kh_value(rg_stats->insert_sizes, k);
    median = insert_sizes[i];
  }

  free(insert_sizes);
  rg_stats->median_insert_size = median;
  return median;
}

uint32_t calc_insert_size_abs_dev(rg_stats_t *rg_stats) {
  khash_t(int32) *dev_hash = kh_init(int32);
  khiter_t k;
  for (k = kh_begin(rg_stats->insert_sizes); k != kh_end(rg_stats->insert_sizes); k++)
    if (kh_exist(rg_stats->insert_sizes, k)) {
      uint32_t deviation = abs(rg_stats->median_insert_size - kh_key(rg_stats->insert_sizes, k));
      int ret;
      khiter_t k_dev = kh_put(int32, dev_hash, deviation, &ret);
      if (!ret)
        kh_value(dev_hash, k_dev) += kh_value(rg_stats->insert_sizes, k);
      else
        kh_value(dev_hash, k_dev) = kh_value(rg_stats->insert_sizes, k);
    }

  uint32_t *dev_arr = malloc(sizeof(uint32_t) * kh_size(dev_hash));
  uint8_t arr_size=0;
  for (k = kh_begin(dev_hash); k != kh_end(dev_hash); k++)
    if (kh_exist(dev_hash, k)) {
      dev_arr[arr_size] = kh_key(dev_hash, k);
      arr_size++;
    }

  qsort(dev_arr, arr_size, sizeof(uint32_t), compare_integers);

  uint32_t median_dev;
  uint32_t num_devs_counted = 0;
  uint32_t stop_counting = (uint32_t) (0.5 * rg_stats->num_insert_sizes);
  uint8_t i;
  for (i=0; num_devs_counted < stop_counting ; i++) {
    k = kh_get(int32, dev_hash, dev_arr[i]);
    num_devs_counted += kh_value(dev_hash, k);
    median_dev = dev_arr[i];
  }

  kh_destroy(int32, dev_hash);
  free(dev_arr);
  return median_dev;
}




void output_rg_stats(void *_rg_stats_hash, char* out_fname, char* bam_basename, unsigned char* md5) {
  khash_t(rg_stats) *rg_stats_hash = (khash_t(rg_stats)*) _rg_stats_hash;
  FILE *fp = fopen(out_fname, "w");
  uint8_t i;
  fprintf(fp, "bam_filename\tmd5\tstudy\tsample\tplatform\tlibrary\treadgroup");
  fprintf(fp, "\t#_total_bases\t#_mapped_bases\t#_total_reads\t#_mapped_reads");
  fprintf(fp, "\t#_mapped_reads_paired_in_sequencing\t#_mapped_reads_properly_paired");
  fprintf(fp, "\t\%_of_mismatched_bases\taverage_quality_of_mapped_bases\tmean_insert_size");
  fprintf(fp, "\tinsert_size_sd\tmedian_insert_size\tinsert_size_median_absolute_deviation");
  fprintf(fp, "\t#_duplicate_reads\t#_duplicate_bases\n");

  khiter_t k;
  for (k = kh_begin(rg_stats_hash); k != kh_end(rg_stats_hash); k++)
    if (kh_exist(rg_stats_hash, k)) {
      rg_stats_t *rg_stats = kh_value(rg_stats_hash, k);
      fprintf(fp, "%s\t", bam_basename);
      for (i=0; i < MD5_DIGEST_LENGTH; i++)
        fprintf(fp, "%02x", md5[i]);
      fprintf(fp, "\t%s", (rg_stats->description ? rg_stats->description : "-"));
      fprintf(fp, "\t%s", (rg_stats->sample ? rg_stats->sample : "-"));
      fprintf(fp, "\t%s", (rg_stats->platform ? rg_stats->platform : "-"));
      fprintf(fp, "\t%s", (rg_stats->library ? rg_stats->library : "-"));
      fprintf(fp, "\t%s", kh_key(rg_stats_hash, k));
      fprintf(fp, "\t%d", rg_stats->num_total_bases);
      fprintf(fp, "\t%d", rg_stats->num_mapped_bases);
      fprintf(fp, "\t%d", rg_stats->num_total_reads);
      fprintf(fp, "\t%d", rg_stats->num_mapped_reads);
      fprintf(fp, "\t%d", rg_stats->num_mapped_reads_paired_in_sequencing);
      fprintf(fp, "\t%d", rg_stats->num_mapped_reads_properly_paired);
      fprintf(fp, "\t%.2f", calc_percentage_mismatched_bases(rg_stats));
      fprintf(fp, "\t%.2f", calc_average_quality_mapped_bases(rg_stats));
      fprintf(fp, "\t%d", (int) calc_mean_insert_size(rg_stats));
      fprintf(fp, "\t%.2f", calc_insert_size_sd(rg_stats));
      fprintf(fp, "\t%d", calc_median_insert_size(rg_stats));
      fprintf(fp, "\t%d", calc_insert_size_abs_dev(rg_stats));
      fprintf(fp, "\t%d", rg_stats->num_duplicate_reads);
      fprintf(fp, "\t%d", rg_stats->num_duplicate_bases);
      fprintf(fp, "\n");
    }
  fclose(fp);
  return;
}

// Takes the input filename and creates a string containing the basename
// The calling function will need to free/destroy this string
char *get_basename(char *full_filename) {
  char *filename = strrchr(full_filename, '/');
  if (! filename)
    filename = full_filename;
  else
    filename ++;

  size_t filename_length = strlen(filename);
  char* extension = filename + filename_length -4;

  char* basename;
  if (filename_length >4 && (strcmp(extension, ".bam") ==0 || strcmp(extension, ".sam") ==0)) {
    basename = malloc(filename_length - 3);
    strncpy(basename, filename, filename_length -4);
    basename[filename_length -4] = '\0';
  }
  else {
    basename = malloc(filename_length +1);
    strcpy(basename, filename);
  }
  return basename;
}

// Reads a file and creates a string containing md5
// The calling function will need to free/destroy this string
unsigned char * calc_md5(char* filename) {
  int inFile = open(filename, O_RDONLY);
  if (inFile < 0) {
    printf("Could not open file\n");
  }

  struct stat statbuf;
  fstat(inFile, &statbuf);

  unsigned char *md5 = malloc(MD5_DIGEST_LENGTH);

  unsigned char* file_buffer = (unsigned char*) mmap(0, statbuf.st_size, PROT_READ, MAP_SHARED, inFile, 0);
  MD5(file_buffer, statbuf.st_size, md5);
  close(inFile);
  return md5;
}

int main(int argc, char *argv[])
{
  char *bam_fname = NULL;
  char *out_fname = NULL;
  char in_mode[5];
  strcpy(in_mode, "rb");

  int c;
  while((c = getopt(argc, argv, "so:")) != -1)
    switch (c) {
      case 'o':
        out_fname = optarg;
        break;
      case 's':
        strcpy(in_mode, "r");
        break;
      case '?':
        usage();
        break;
    }

  if (optind < argc)
    bam_fname = argv[optind++];

  if ( !bam_fname ) {
    if ( isatty(fileno((FILE *)stdin)) )
        usage();
    bam_fname = "-";
  }

  if ( !out_fname ) {
    usage();
  }

  char* bam_basename = get_basename(bam_fname);
  unsigned char* md5 = calc_md5(bam_fname);

  bamFile bam = NULL;
  if ((bam = bam_open(bam_fname, in_mode)) == 0)
    printf("Failed to open file %s\n", bam_fname);
  bam_header_t *bam_header = bam_header_read(bam);

  khash_t(rg_stats) *rg_stats_hash = rg_stats_init(bam_header);
  khiter_t k;
  int ret;


  bam1_t *bam_line = bam_init1();
  while (bam_read1(bam, bam_line) > 0) {
    uint8_t *rg_with_type = bam_aux_get(bam_line, "RG");
    k = kh_put(rg_stats, rg_stats_hash, bam_aux2Z(rg_with_type), &ret);
    rg_stats_t *rg_stats = kh_value(rg_stats_hash, k);

    update_stats(rg_stats, bam_line);

  }

  bam_destroy1(bam_line);
  bam_close(bam);

  output_rg_stats(rg_stats_hash, out_fname, bam_basename, md5);
  rg_stats_destroy(rg_stats_hash);
  bam_header_destroy(bam_header);
  free(bam_basename);
  free(md5);

  return 0;
}


void usage()
{
  printf("usage: validate_bam -o output_file [-s] input_file\n");
  printf("-s flag means input is in sam format\n");
  exit(1);
}
