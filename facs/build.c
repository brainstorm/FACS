#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>

#include "build.h"
#include "bloom.h"
#include "file_dir.h"
#include "tool.h"

//@lh3 FASTA/FASTQ reading magic
#include <zlib.h>
#include "readfq/kseq.h"
KSEQ_INIT(gzFile, gzread);

static int
build_usage (void)
{
  fprintf (stderr, "\nUsage: ./facs build [options]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "\t-r reference FASTA/FASTQ file\n");
  fprintf (stderr, "\t-o output bloom filter file\n");
  fprintf (stderr,
	   "\t-l text file containing all reference files, will build individual bloom filters for each\n");
  fprintf (stderr,
	   "\t-k k-mer size, default is automatically estimated from the reference file\n");
  fprintf (stderr,
	   "\t-e allowed false positive frequency, default is 0.005\n");
  fprintf (stderr, "\n");
  return 1;
}

int
build_main (int argc, char **argv)
{
  if (argc < 2)
    build_usage ();

  BIGNUM capacity;

  int opt;
  int k_mer = 0;
  float error_rate = 0.005;

  char *bloom_file = NULL;
  char *reference = NULL;

  while ((opt = getopt (argc, argv, "e:k:o:r:h")) != -1) {
      switch (opt) {
        case 'e':
          (optarg) && (error_rate = atof (optarg), 1);
          break;
        case 'k':
          (optarg) && (k_mer = atoi (optarg), 1);
          break;
        case 'o':
          (optarg) && (bloom_file = optarg, 1);
          break;
        case 'r':
          (optarg) && (reference = optarg, 1);
          break;
        case 'h':
          return build_usage ();
        default:
          printf ("Unknown option: -%c\n", (char) optopt);
          return build_usage ();
	  }
  }

  // Nothing to do here if no bloom file is specified
  // XXX create same input filename but with ".bloom" extension?
  if (!bloom_file) {
      fprintf(stderr, "error: No bloom file supplied (use -o)\n");
      return build_usage();
  }

  build (reference, bloom_file, k_mer, error_rate, NULL);
  return 0;
}



int
build (char *ref_name, char *bloom_file, int k_mer, double error_rate, char *prefix)
{
  gzFile fp;
  kseq_t *seq = NULL;

  bloom *bl = NEW (bloom);

  if (k_mer != 0)
    bl->k_mer = k_mer;
  else
    bl->k_mer = kmer_suggestion (get_filesize(ref_name));

  bl->stat.e = error_rate;
  bl->dx = dx_add(bl->k_mer);
  bl->stat.capacity = get_filesize(ref_name);
  get_rec(&bl->stat);

  bloom_init(bl, bl->stat.elements, bl->stat.capacity,
             bl->stat.e, bl->stat.ideal_hashes, NULL, 3);

  // Files or stdin input
  if(!ref_name)
    fp = gzdopen(fileno(stdin), "r");
  else
  	fp = gzdopen(open(ref_name, O_RDONLY), "r");

  // Init uncompressed reads stream
  seq = kseq_init(fp);

  int cur;

  // Iterates through reads (FASTA/FASTQ).
  while (kseq_read(seq) >= 0) {
    // adding k_mers to a bloom filter via sliding window.
    cur = 0;
    while (cur <= (seq->seq.l - bl->k_mer)) {
        //printf("%s\n", strndup(seq->seq.s + cur, bl->k_mer));
        bloom_add(bl, strndup(seq->seq.s + cur, bl->k_mer), bl->k_mer);
        cur++;
    }
  }

  if(ref_name)
    save_bloom(ref_name, bl, NULL, bloom_file);

  return 0;
}
