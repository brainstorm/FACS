#include <zlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>
#include "tool.h"
#include "bloom.h"
#include "file_dir.h"
#include "query.h"

//@lh3 FASTA/FASTQ reading magic
#include "readfq/kseq.h"
KSEQ_INIT(gzFile, gzread);

#ifndef __clang__
  // openMP not yet ported to clang: http://www.phoronix.com/scan.php?page=news_item&px=MTI2MjU
  #include <omp.h>
#endif

static int
query_usage (void)
{
  fprintf (stderr, "\nUsage: ./facs query [options]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "\t-r reference Bloom filter to query against\n");
  fprintf (stderr, "\t-q FASTA/FASTQ file containing the query\n");
  fprintf (stderr, "\t-l input list containing all Bloom filters, one per line\n");
  fprintf (stderr, "\t-t threshold value\n");
  fprintf (stderr, "\t-f report output format, valid values are: 'json' and 'tsv'\n");
  fprintf (stderr, "\t-s sampling rate, default is 1 so it reads the whole query file\n");
  fprintf (stderr, "\n");
  return 1;
}

int
bq_main (int argc, char **argv)
{
  if (argc < 3)
    return query_usage ();

/*-------defaults for bloom filter building-------*/
  int opt;
  double tole_rate = 0;
  double sampling_rate = 1;

  char *ref = NULL;
  char *list = NULL;
  char *target_path = NULL;
  char *source = NULL;
  char *report_fmt = NULL;

  // XXX: make r and l mutually exclusive
  while ((opt = getopt (argc, argv, "s:t:r:o:q:l:f:h")) != -1) {
      switch (opt) {
	case 't':
	  (optarg) && ((tole_rate = atof (optarg)), 1);
	  break;
	case 's':
	  (optarg) && ((sampling_rate = atof (optarg)), 1);
	  break;
	case 'o':
	  (optarg) && ((target_path = optarg), 1);
	  break;
	case 'q':
	  (optarg) && (source = optarg, 1);
	  break;
	case 'r':
	  (optarg) && (ref = optarg, 1);
	  break;
	case 'l':
	  (optarg) && (list = optarg, 1);
	  break;
	case 'f': // "json", "tsv" or none
	  (optarg) && (report_fmt = optarg, 1);
	  break;
	case 'h':
	  return query_usage ();
	case '?':
	  printf ("Unknown option: -%c\n", (char) optopt);
	  return query_usage ();
	}
  }

  if (!target_path && !source) {
      fprintf (stderr, "\nPlease, at least specify a bloom filter (-r) and a query file (-q)\n");
      exit (-1);
  }

  if (target_path == NULL) {
      target_path = argv[0];
  }

  return query(source, ref, tole_rate, sampling_rate, list,
               target_path, report_fmt);
}

int
query (char *query, char *bloom_filter, double tole_rate, double sampling_rate,
       char *list, char *target_path, char *report_fmt)
{
  /* TODO: Implement sampling
   *       OpenMP it
   *       MPI it
   */
  gzFile fp;
  kseq_t *seq;

  bloom *bl = NEW (bloom);
  F_set *File_head = make_list (bloom_filter, list);
  int read = 0;

  File_head->reads_num = 0;
  File_head->reads_contam = 0;
  File_head->hits = 0;
  File_head->filename = bloom_filter;
 
  load_bloom(File_head->filename, bl);

  if(!query) {
  	fp = gzdopen(fileno(stdin), "r");
  } else {
  	fp = gzdopen(open(query, O_RDONLY), "r");
  }

  seq = kseq_init(fp);
  
  if (tole_rate == 0)
    tole_rate = mco_suggestion (bl->k_mer);

  while (kseq_read(seq) >= 0) {
    File_head->reads_num++;

    read = query_read(seq->seq.s, seq->seq.l, 'n', bl, tole_rate, File_head);
    if(read){
    	File_head->reads_contam++;
    }
  }

  report(File_head, query, report_fmt, target_path);

  kseq_destroy(seq);
  gzclose(fp);
  bloom_destroy(bl);
  return 0;
}

char *
strrstr (char *s, char *str)
{
/* Matches string in reverse */
  char *p;
  int len = strlen (s);
  for (p = s + len - 1; p >= s; p--)
    {
      if ((*p == *str) && (memcmp (p, str, strlen (str)) == 0))
	return p;
    }
  return NULL;
}

void
clean_list (Queue * head, Queue * tail)
{
  Queue *element;
  while (head != tail)
    {
      element = head->next;
      memset (head, 0, sizeof (Queue));
      free (head);
      head = element;
    }
  free (tail);
}

char *
bac_2_n (char *filename)
{
  while (*filename != '\n')
    filename--;
  filename--;			//move from \n
  while (*filename != '\n')
    filename--;
  filename++;
  return filename;
}
