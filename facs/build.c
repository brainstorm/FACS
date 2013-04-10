#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
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

  char *position;
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

  build (reference, bloom_file, k_mer, error_rate, NULL);
  return 0;
}



int
build (char *ref_name, char *bloom_file, int k_mer, double error_rate, char *prefix)
{
  char* read_chunk;
  char* position = mmaping(ref_name);
  float parts;
  gzFile fp;
  kseq_t *seq = NULL;
  
  bloom *bl = NEW (bloom);
 
  if (k_mer != 0)
    bl->k_mer = k_mer;
  else
    bl->k_mer = kmer_suggestion (get_size (ref_name));

  bl->stat.e = error_rate;
  bl->dx = dx_add(bl->k_mer);
  bl->stat.capacity = strlen(position);
  get_rec(&bl->stat);

  bloom_init(bl, bl->stat.elements, bl->stat.capacity,
             bl->stat.e, bl->stat.ideal_hashes, NULL, 3);

  // Files or stdin input
  if(!ref_name || !bloom_file)
    fp = gzdopen(fileno(stdin), "r");
  else {
  	fp = gzdopen(open(ref_name, O_RDONLY), "r");
  }

  seq = kseq_init(fp);

  //k_mer sliding window counter
  int k_mer_window = 0;

  // Reads sequences, partitions them in k_mer size and
  // adds them to the bloom filter
  while (kseq_read(seq) >= 0) {
    parts = seq->seq.l/bl->k_mer;

    while (parts <= seq->seq.l) {
        read_chunk = substr(seq->seq.s, k_mer_window,
                            sizeof(bl->k_mer));
        printf("%f,%d,%s\n", parts, bl->k_mer, read_chunk);
        bloom_add(bl, read_chunk, seq->seq.l);
        k_mer_window++;
        parts++;
    }
  }

  if(ref_name)
    save_bloom(ref_name, bl, NULL, bloom_file);

  return 0;
}
