#include <stdio.h>
#include <unistd.h> /* for getopt */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sam.h>
#include <sam_header.h>
#include <khash.h>
#include <openssl/md5.h> /* for md5 */


void usage();

KHASH_MAP_INIT_INT(int32, uint32_t) /* used to store insert sizes */
KHASH_MAP_INIT_STR(tags, char*) /* used to store bam header tags */

/* struct containing statistics for a read group */
typedef struct
{
  char *platform_unit; /* points to a bam_header string */
  char *library; /* points to a bam_header string */
  char *sample; /* points to a bam_header string */
  char *platform; /* points to a bam_header string */
  char *description; /* points to a bam_header string */

  uint64_t num_total_bases; /* The sum of the length of all reads in this readgroup */
  uint64_t num_mapped_bases; /* The total number of mapped bases for all reads in this readgroup that did not have flag 4 (== unmapped) */
  uint64_t num_bases_in_mapped_reads; /* The total number of bases for all reads in this readgroup that did not have flag 4 (== unmapped) */
  uint64_t num_total_reads; /* The total number of reads in this readgroup */
  uint64_t num_mapped_reads; /* The total number of reads in this readgroup that did not have flag 4 (== unmapped) */
  uint64_t num_mapped_reads_paired_in_sequencing; /* Same as num_mapped_reads, also requiring flag 1 (== reads paired in sequecing) */
  uint64_t num_mapped_reads_properly_paired; /* Same as num_mapped_reads, also requiring flag 2 (== mapped in a proper pair, inferred during alignment). */
  uint64_t num_NM_bases; /* The sum of the length of all reads in this readgroup that have a NM tag */
  uint64_t num_NM_mismatched_bases; /* The sum of the value of the NM tags for all reads in this readgroup */
  uint64_t sum_quality_mapped_bases; /* Sum of the base qualities for all bases for all reads in this readgroup that did not have the flag 4 (== unmapped) */
  uint64_t num_duplicate_reads; /* Number of reads which were marked as duplicates */
  uint64_t num_duplicate_bases; /* The sum of the length of all reads which were marked as duplicates */

  /* Insert sizes, counts all insert sizes greater than 0 for properly paired reads (requires flag 2) and with a mapping quality greater than 0 */
  khash_t(int32) *insert_sizes; /* key is the insert size, value is the number of occurrences in this read group */
  uint32_t num_insert_sizes; /* The number of reads in this read group meeting the criteria */
  uint64_t sum_insert_sizes; /* Sum of all valid insert sizes for this read group */

  /* The following are calculated after the entire bam file has been read  */
  uint16_t median_insert_size;
  double percent_mismatched_bases; /* = num_NUM_mismatched_bases / num_NM_bases */
  double avg_quality_mapped_bases; /* = sum_quality_mapped_bases / num_total_bases */
  double mean_insert_size; /* = sum_insert_sizes / num_insert_sizes */
  double insert_size_sd; /* the standard deviation from the mean of insert sizes */
  uint16_t insert_size_median_abs_dev; /* The median absolute deviation of insert sizes */
}
rg_stats_t;

KHASH_MAP_INIT_STR(rg_stats, rg_stats_t*)

/* Initialise a khash containing: key is read group ID, value is an initialised rg_stats_t */
void *rg_stats_init(bam_header_t *bam_header) {
  khash_t(rg_stats) *rg_stats_hash = kh_init(rg_stats);
  int n_rg, i, ret;
  char **ID_array;
  khash_t(tags) *LB_hash, *SM_hash, *PL_hash, *DS_hash, *PU_hash;
  khiter_t k;

  bam_header->dict = sam_header_parse2(bam_header->text);
  ID_array = sam_header2list(bam_header->dict, "RG", "ID", &n_rg);
  LB_hash = sam_header2tbl(bam_header->dict, "RG", "ID", "LB");
  SM_hash = sam_header2tbl(bam_header->dict, "RG", "ID", "SM");
  PL_hash = sam_header2tbl(bam_header->dict, "RG", "ID", "PL");
  DS_hash = sam_header2tbl(bam_header->dict, "RG", "ID", "DS");
  PU_hash = sam_header2tbl(bam_header->dict, "RG", "ID", "PU");

  for(i=0; i<n_rg; i++) {
    rg_stats_t *new_rg_stats = (rg_stats_t*) malloc(sizeof(rg_stats_t));
    if (new_rg_stats == NULL) {
      fprintf(stderr, "Out of memory\n");
      exit(-1);
    }

    k = kh_put(rg_stats, rg_stats_hash, ID_array[i], &ret);
    kh_value(rg_stats_hash, k) = new_rg_stats;

    k = kh_get(tags, LB_hash, ID_array[i]);
    new_rg_stats->library = (k == kh_end(LB_hash)) ? NULL : kh_value(LB_hash, k);

    k = kh_get(tags, SM_hash, ID_array[i]);
    new_rg_stats->sample = (k == kh_end(SM_hash)) ? NULL : kh_value(SM_hash, k);

    k = kh_get(tags, PL_hash, ID_array[i]);
    new_rg_stats->platform = (k == kh_end(PL_hash)) ? NULL : kh_value(PL_hash, k);

    k = kh_get(tags, DS_hash, ID_array[i]);
    new_rg_stats->description = (k == kh_end(DS_hash)) ? NULL : kh_value(DS_hash, k);

    k = kh_get(tags, PU_hash, ID_array[i]);
    new_rg_stats->platform_unit = (k == kh_end(PU_hash)) ? NULL : kh_value(PU_hash, k);

    new_rg_stats->insert_sizes = kh_init(int32);
    new_rg_stats->num_total_bases = 0;
    new_rg_stats->num_mapped_bases = 0;
    new_rg_stats->num_bases_in_mapped_reads = 0;
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
  kh_destroy(tags, PU_hash);
  return rg_stats_hash;
}

/* destroy the khash of rg_stats.  Also destroy each rg_stats. */
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

/* update the stats for a single alignment */
void update_stats(rg_stats_t *rg_stats, bam1_t *bam_line) {
  rg_stats->num_total_reads ++;
  rg_stats->num_total_bases += bam_line->core.l_qseq;


  if (! (bam_line->core.flag & BAM_FUNMAP)) {
    uint32_t *cigar = bam1_cigar(bam_line);
    uint8_t *qual = bam1_qual(bam_line);
    uint8_t i,j;

    rg_stats->num_mapped_reads ++;
    rg_stats->num_bases_in_mapped_reads += bam_line->core.l_qseq;
    for (i=0; i < bam_line->core.n_cigar; i++){
      int op = bam_cigar_op(cigar[i]);
      int oplen = bam_cigar_oplen(cigar[i]);
      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
        rg_stats->num_mapped_bases += oplen;
      if (op != BAM_CDEL && op != BAM_CREF_SKIP) {
        for (j=0; j<oplen; j++)
          rg_stats->sum_quality_mapped_bases += qual[j];
        qual += oplen;
      }

    }

    if (bam_line->core.flag & BAM_FPAIRED)
      rg_stats->num_mapped_reads_paired_in_sequencing ++;

    if (bam_line->core.flag & BAM_FPROPER_PAIR) {
      rg_stats->num_mapped_reads_properly_paired ++;
      if (bam_line->core.qual > 0 && bam_line->core.isize > 0) {
        int ret;
        khiter_t k;

        rg_stats->num_insert_sizes ++;
        rg_stats->sum_insert_sizes += bam_line->core.isize;
        k = kh_put(int32, rg_stats->insert_sizes, bam_line->core.isize, &ret);
        if (!ret)
          kh_value(rg_stats->insert_sizes, k) ++;
        else
          kh_value(rg_stats->insert_sizes, k) = 1;
      }
    }

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


/* calculate the mean_insert_size (after the entire bam file has been read) */
void calc_insert_size_sd(rg_stats_t *rg_stats) {
  double mean = rg_stats->mean_insert_size;
  double variance = 0;
  khiter_t k;

  if (! isfinite(mean)) {
    rg_stats->insert_size_sd = 0.0;
    return;
  }

  for (k = kh_begin(rg_stats->insert_sizes); k != kh_end(rg_stats->insert_sizes); k++)
    if (kh_exist(rg_stats->insert_sizes, k)) {
      double diff = kh_key(rg_stats->insert_sizes, k) - mean;
      variance += diff * diff * kh_value(rg_stats->insert_sizes, k);
    }
  variance /= rg_stats->num_insert_sizes;
  rg_stats->insert_size_sd = sqrt(variance);
  return;
}

/* needed for the qsort function */
int compare_integers (const void * a, const void *b) {
  return ( *(uint16_t*)a - *(uint16_t*)b );
}

/* calculate the median_insert_size (after the entire bam file has been read) */
void calc_median_insert_size(rg_stats_t *rg_stats) {
  khiter_t k;
  uint16_t arr_size, i;
  uint16_t median = 0;
  uint32_t num_inserts_counted = 0;
  uint32_t stop_counting = (uint32_t) (0.5 * rg_stats->num_insert_sizes);

  uint16_t *insert_sizes = (uint16_t*) malloc(sizeof(uint16_t) * kh_size(rg_stats->insert_sizes));
  if (insert_sizes == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(-1);
  }

  for (arr_size=0, k = kh_begin(rg_stats->insert_sizes); k != kh_end(rg_stats->insert_sizes); k++)
    if (kh_exist(rg_stats->insert_sizes, k)) {
      insert_sizes[arr_size] = kh_key(rg_stats->insert_sizes, k);
      arr_size++;
    }

  qsort(insert_sizes, arr_size, sizeof(uint16_t), compare_integers);

  for (i=0; num_inserts_counted < stop_counting ; i++) {
    k = kh_get(int32, rg_stats->insert_sizes, insert_sizes[i]);
    num_inserts_counted += kh_value(rg_stats->insert_sizes, k);
    median = insert_sizes[i];
  }

  free(insert_sizes);
  rg_stats->median_insert_size = median;
  return;
}

/* calculate the insert size median absolute deviation (after the entire bam file has been read) */
void calc_insert_size_abs_dev(rg_stats_t *rg_stats) {
  khash_t(int32) *dev_hash = kh_init(int32);
  khiter_t k;
  uint8_t arr_size, i;
  uint16_t median_dev = 0;
  uint32_t num_devs_counted = 0;
  uint32_t stop_counting = (uint32_t) (0.5 * rg_stats->num_insert_sizes);
  uint16_t *dev_arr;

  for (k = kh_begin(rg_stats->insert_sizes); k != kh_end(rg_stats->insert_sizes); k++)
    if (kh_exist(rg_stats->insert_sizes, k)) {
      uint16_t deviation = abs(rg_stats->median_insert_size - kh_key(rg_stats->insert_sizes, k));
      int ret;
      khiter_t k_dev = kh_put(int32, dev_hash, deviation, &ret);
      if (!ret)
        kh_value(dev_hash, k_dev) += kh_value(rg_stats->insert_sizes, k);
      else
        kh_value(dev_hash, k_dev) = kh_value(rg_stats->insert_sizes, k);
    }

  dev_arr = (uint16_t*) malloc(sizeof(uint16_t) * kh_size(dev_hash));
  if (dev_arr == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(-1);
  }
  for (arr_size=0, k = kh_begin(dev_hash); k != kh_end(dev_hash); k++)
    if (kh_exist(dev_hash, k)) {
      dev_arr[arr_size] = kh_key(dev_hash, k);
      arr_size++;
    }

  qsort(dev_arr, arr_size, sizeof(uint16_t), compare_integers);

  for (i=0; num_devs_counted < stop_counting ; i++) {
    k = kh_get(int32, dev_hash, dev_arr[i]);
    num_devs_counted += kh_value(dev_hash, k);
    median_dev = dev_arr[i];
  }

  kh_destroy(int32, dev_hash);
  free(dev_arr);
  rg_stats->insert_size_median_abs_dev = median_dev;
  return;
}

/* calculate statistics that can only be calculated after the entire bam file has been read */
void calc_final_rg_stats(void *_rg_stats_hash) {
  khash_t(rg_stats) *rg_stats_hash = (khash_t(rg_stats)*) _rg_stats_hash;
  khiter_t k;
  for (k = kh_begin(rg_stats_hash); k != kh_end(rg_stats_hash); k++)
    if (kh_exist(rg_stats_hash, k)) {
      rg_stats_t *rg_stats = kh_value(rg_stats_hash, k);

      rg_stats->percent_mismatched_bases = 100 * (double) rg_stats->num_NM_mismatched_bases / rg_stats->num_NM_bases;
      rg_stats->avg_quality_mapped_bases =  (double) rg_stats->sum_quality_mapped_bases / rg_stats->num_bases_in_mapped_reads;
      rg_stats->mean_insert_size =  (double) rg_stats->sum_insert_sizes / rg_stats->num_insert_sizes;

      calc_insert_size_sd(rg_stats);
      calc_median_insert_size(rg_stats);
      calc_insert_size_abs_dev(rg_stats);
    }
}




/* opens a file, writes statistics for all read groups */
void output_rg_stats(void *_rg_stats_hash, char* out_fname, char* bam_basename, char* md5) {
  khash_t(rg_stats) *rg_stats_hash = (khash_t(rg_stats)*) _rg_stats_hash;
  khiter_t k;
  FILE *fp = fopen(out_fname, "w");
  if (fp == NULL) {
    fprintf(stderr, "Could not open file for writing: %s\n", out_fname);
    exit(-1);
  }

  fprintf(fp, "bam_filename\tmd5\tstudy\tsample\tplatform\tlibrary\treadgroup");
  fprintf(fp, "\t#_total_bases\t#_mapped_bases\t#_total_reads\t#_mapped_reads");
  fprintf(fp, "\t#_mapped_reads_paired_in_sequencing\t#_mapped_reads_properly_paired");
  fprintf(fp, "\t%c_of_mismatched_bases\taverage_quality_of_mapped_bases\tmean_insert_size", '%');
  fprintf(fp, "\tinsert_size_sd\tmedian_insert_size\tinsert_size_median_absolute_deviation");
  fprintf(fp, "\t#_duplicate_reads\t#_duplicate_bases\n");

  for (k = kh_begin(rg_stats_hash); k != kh_end(rg_stats_hash); k++)
    if (kh_exist(rg_stats_hash, k)) {
      uint8_t i;
      rg_stats_t *rg_stats = kh_value(rg_stats_hash, k);
      fprintf(fp, "%s", bam_basename);
      fprintf(fp, "\t%s", md5);
      fprintf(fp, "\t%s", (rg_stats->description ? rg_stats->description : "unknown_study"));
      fprintf(fp, "\t%s", (rg_stats->sample ? rg_stats->sample : "-"));
      fprintf(fp, "\t%s", (rg_stats->platform ? rg_stats->platform : "-"));
      fprintf(fp, "\t%s", (rg_stats->library ? rg_stats->library : "-"));
      fprintf(fp, "\t%s", (rg_stats->platform_unit ? rg_stats->platform_unit : kh_key(rg_stats_hash, k)));
      fprintf(fp, "\t%llu", rg_stats->num_total_bases);
      fprintf(fp, "\t%llu", rg_stats->num_mapped_bases);
      fprintf(fp, "\t%llu", rg_stats->num_total_reads);
      fprintf(fp, "\t%llu", rg_stats->num_mapped_reads);
      fprintf(fp, "\t%llu", rg_stats->num_mapped_reads_paired_in_sequencing);
      fprintf(fp, "\t%llu", rg_stats->num_mapped_reads_properly_paired);
      fprintf(fp, "\t%.2f", rg_stats->percent_mismatched_bases);
      fprintf(fp, "\t%.2f", rg_stats->avg_quality_mapped_bases);
      fprintf(fp, "\t%d", isfinite(rg_stats->mean_insert_size) ? (int) rg_stats->mean_insert_size : 0);
      fprintf(fp, "\t%.2f", rg_stats->insert_size_sd);
      fprintf(fp, "\t%d", rg_stats->median_insert_size);
      fprintf(fp, "\t%d", rg_stats->insert_size_median_abs_dev);
      fprintf(fp, "\t%llu", rg_stats->num_duplicate_reads);
      fprintf(fp, "\t%llu", rg_stats->num_duplicate_bases);
      fprintf(fp, "\n");
    }
  fclose(fp);
  return;
}

/* Takes the input filename and creates a string containing the basename
The calling function will need to free/destroy this string */
char *get_basename(char *full_filename) {
  size_t filename_length, basename_length;
  char *filename, *basename, *extension;

  filename = strrchr(full_filename, '/');
  if (! filename)
    filename = full_filename;
  else
    filename ++;

  filename_length = strlen(filename);
  extension = filename + filename_length -4;

  if (filename_length >4 && (strcmp(extension, ".bam") ==0 || strcmp(extension, ".sam") ==0))
    basename_length = filename_length -4;
  else
    basename_length = filename_length;

  basename = (char*) malloc(basename_length + 1);
  if (basename == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(-1);
  }
  strncpy(basename, filename, basename_length);
  basename[basename_length] = '\0';

  return basename;
}

void md5_to_str(unsigned char* md5, char* md5_string) {
  char* current_position = md5_string;
  int i;
  for(i=0; i< MD5_DIGEST_LENGTH; i++) {
    sprintf(current_position, "%02x", md5[i]);
    current_position += 2;
  }
}

/* Reads a file and creates md5
Writes md5 as a string to md5_string */
void calc_md5(char* filename, char* md5_string) {
  unsigned char md5[MD5_DIGEST_LENGTH], file_buffer[1024];
  MD5_CTX md_context;
  int bytes;

  FILE *inFile = fopen(filename, "rb");
  if (inFile == NULL) {
    fprintf(stderr, "Could not open file for reading: %s\n", filename);
    exit(-1);
  }

  MD5_Init(&md_context);
  while ((bytes = fread (file_buffer, 1, 1024, inFile)) != 0)
    MD5_Update (&md_context, file_buffer, bytes);
  MD5_Final (md5, &md_context);
  if (fclose(inFile) != 0) {
    fprintf(stderr,"Could not close file\n");
    exit(-1);
  }

  md5_to_str(md5, md5_string);

  return;
}

int main(int argc, char *argv[])
{
  char *bam_fname = NULL;
  char *out_fname = NULL;
  char in_mode[5];
  int c;
  char* bam_basename;
  char md5[33];
  bamFile bam = NULL;
  bam_header_t *bam_header = NULL;
  khash_t(rg_stats) *rg_stats_hash;
  bam1_t *bam_line = bam_init1();

  strcpy(in_mode, "rb");
  md5[0]='\0';

  while((c = getopt(argc, argv, "so:m:")) != -1)
    switch (c) {
      case 'o':
        out_fname = optarg;
        break;
      case 's':
        strcpy(in_mode, "r");
        break;
      case 'm':
        if (strlen(optarg) != 32) {
          fprintf(stderr, "invalid md5: %s\n", optarg);
          exit(-1);
        }
        strcpy(md5, optarg);
        break;
      case '?':
        usage();
        break;
    }

  if (optind < argc)
    bam_fname = argv[optind++];
  else
    usage();

  if ( !out_fname ) {
    usage();
  }

  bam_basename = get_basename(bam_fname);
  if (! strlen(md5))
    calc_md5(bam_fname, md5);

  if ((bam = bam_open(bam_fname, in_mode)) == 0)
  {
    printf("Failed to open file %s\n", bam_fname);
    exit(-1);
  }
  bam_header = bam_header_read(bam);

  rg_stats_hash = rg_stats_init(bam_header);



  while (bam_read1(bam, bam_line) > 0) {
    uint8_t *rg_with_type = bam_aux_get(bam_line, "RG");
    khiter_t k = kh_get(rg_stats, rg_stats_hash, bam_aux2Z(rg_with_type));
    if (k == kh_end(rg_stats_hash)) {
      fprintf(stderr, "Unknown read group in bam: %s\n", bam_aux2Z(rg_with_type));
      exit(-1);
    }
    rg_stats_t *rg_stats = kh_value(rg_stats_hash, k);

    update_stats(rg_stats, bam_line);

  }

  bam_destroy1(bam_line);
  bam_close(bam);

  calc_final_rg_stats(rg_stats_hash);

  output_rg_stats(rg_stats_hash, out_fname, bam_basename, md5);
  rg_stats_destroy(rg_stats_hash);
  bam_header_destroy(bam_header);
  free(bam_basename);

  return 0;
}


void usage()
{
  printf("usage: validate_bam -o output_file [-s] -m md5 input_file\n");
  printf("-s flag means input is in sam format\n");
  printf("md5 will be calculated if the -m flag is not given\n");
  exit(-1);
}
