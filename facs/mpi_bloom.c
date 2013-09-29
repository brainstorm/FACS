#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <math.h>
#include <time.h>
#include <errno.h>
#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <string.h>

#include<fcntl.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/time.h>
#include<sys/mman.h>
#include<sys/types.h>

#include "tool.h"
#include "prob.h"
#include "bloom.h"
#include "big_query.h"
#include "check.h"
#include "hashes.h"
#include "mpi_bloom.h"
#include<omp.h>
#include<mpi.h>
/*-------------------------------------*/
static int mpicheck_usage (void)
{
  fprintf (stderr, "\nUsage: mpirun -n Nr_of_nodes ./facs [options]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "\t-r reference Bloom filter to query against\n");
  fprintf (stderr, "\t-q FASTA/FASTQ file containing the query\n");
  fprintf (stderr, "\t-l input list containing all Bloom filters,\
           one per line\n");
  fprintf (stderr, "\t-t threshold value\n");
  fprintf (stderr, "\t-f report output format, valid values are:\
           'json' and 'tsv'\n");
  fprintf (stderr, "\t-s sampling rate, default is 1 so it reads the whole\
           query file\n");
  fprintf (stderr, "\n");
  exit(1);
}
/*-------------------------------------*/
main (int argc, char **argv)
{
/*------------variables----------------*/
  double tole_rate = 0, sampling_rate = 1;
  char *bloom_filter = NULL, *list = NULL, *target_path = NULL, *position = NULL, *query = NULL, *report_fmt = "json";
  int opt=0, ntask = 0, mytask = 0, exit_sign = 0;
  BIGCAST share=0, offset=0;
  char type = '@';
  gzFile zip = NULL;
  bloom *bl_2 = NEW (bloom);
  Queue *head = NEW (Queue), *tail = NEW (Queue);
  head->location = NULL;
  head->next = tail;
  /*----------MPI initialize------------*/
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &ntask);
  MPI_Comm_size (MPI_COMM_WORLD, &mytask);
/*------------get opt------------------*/
  
  if (argc<3)
  {
	if (ntask == 0)
	{
        	return mpicheck_usage();
	}
  }
  while ((opt = getopt (argc, argv, "s:t:r:o:q:l:f:h")) != -1)
  {
      switch (opt)
      {
        case 't':
          tole_rate = atof(optarg);
          break;
        case 's':
          sampling_rate = atof(optarg);
          break;
        case 'o':
          target_path = optarg;
          break;
        case 'q':
          query = optarg;
          break;
        case 'r':
          bloom_filter = optarg;
          break;
        case 'l':
          list = optarg;
          break;
        case 'f': // "json", "tsv" or none
          (optarg) && (report_fmt = optarg, 1);
          break;
        case 'h':
	  if (ntask==0)
          	return mpicheck_usage();
          break;
        case '?':
          printf ("Unknown option: -%c\n", (char) optopt);
          if (ntask==0)
	  	return mpicheck_usage();
          break;
      }
  }
  if (!bloom_filter && !query)
  {
	if (ntask==0)
  		fprintf (stderr, "\nPlease, at least specify a bloom filter (-r) and a query file (-q)\n");
	MPI_Finalize ();
	exit (EXIT_FAILURE);
  }
  if (target_path == NULL)
  {
        target_path = argv[0];
  }
  if ((zip = gzopen (query, "rb")) < 0)
  {
	if (ntask==0)
	{
        	fprintf(stderr, "%s\n", strerror(errno));
	}
		MPI_Finalize ();
        	exit(EXIT_FAILURE);
  }
  if (strstr (query, ".fastq") != NULL || strstr (query, ".fq") != NULL)
        type = '@';
  else
        type = '>';
  /*initialize emtpy string for query*/
  position = (char *) calloc ((2*ONEG + 1), sizeof (char));
  share = struc_init (query,ntask,mytask);
  F_set *File_head = make_list (bloom_filter, list);
  File_head->reads_num = 0;
  File_head->reads_contam = 0;
  File_head->hits = 0;
  File_head->all_k = 0;
  File_head->filename = bloom_filter;
  load_bloom (File_head->filename, bl_2);	//load a bloom filter
  if (tole_rate == 0)
  {
  	tole_rate = mco_suggestion (bl_2->k_mer);
  }
  while (offset<share)
  {
	if ((share-offset)<2*ONEG)
		exit_sign = 1;
	//printf ("offset->%lld left->%lld\n",offset+share*ntask,share-offset);
        offset+= gz_mpi (zip, offset+share*ntask, share-offset, position, type);
	// put offset += ntask* share inside
	get_parainfo (position, head, type);
        //head = head->next;
#pragma omp parallel
  	{
#pragma omp single nowait
	{
	  	while (head != tail)
		{
#pragma omp task firstprivate(head)
		{
	      		if (head->location != NULL)
                	{
                        	read_process (bl_2, head, tail, File_head, sampling_rate, tole_rate, 'c', type);
			}
		}
	      	head = head->next;
		}
      	}	
	}
     		memset (position, 0, strlen(position));
      		if (exit_sign==1)
			break;
  }
  printf ("finish processing...\n");
  MPI_Barrier (MPI_COMM_WORLD);	//wait until all nodes finish
  gather (File_head,mytask,ntask);			//gather info from all nodes
  if (mytask == 0)		
  {
  	char *result =  report(File_head, query, report_fmt, target_path, prob_suggestion(bl_2->k_mer));
  	printf("%s\n",result);
  }
  MPI_Finalize ();
  return 0;
}
/*-------------------------------------*/
BIGCAST struc_init (char *filename, int ntask, int mytask)
{
  
  BIGCAST total_size = get_size(filename);
  BIGCAST share = 0;
  share = total_size / ntask;	//every task gets an euqal piece
  if (total_size%mytask!=0 && ntask==(mytask-1))
  {
  	share += (total_size % mytask);	//last node takes extra job
  }
  printf("share->%lld\n",share);
  return share;
}
/*-------------------------------------*/
/*current sacrifice file mapping, use gzip instead*/
/*
char *ammaping (char *source)
{
  int src;
  char *sm;

  if ((src = open (source, O_RDONLY | O_LARGEFILE)) < 0)
    {
      perror (" open source ");
      exit (EXIT_FAILURE);
    }

  if (fstat (src, &statbuf) < 0)
    {
      perror (" fstat source ");
      exit (EXIT_FAILURE);
    }

  printf ("share->%d PAGES per node\n", share);

  if (share >= CHUNK)
    buffer = CHUNK;
  else
    buffer = share;
  printf ("total pieces->%d\n", total_piece);
  printf ("PAGE->%d\n", PAGE);
  printf ("node %d chunk size %d buffer size %d offset %d\n", mytask, CHUNK,
	  buffer, offset);

  sm = mmap (0, buffer * PAGE, PROT_READ, MAP_SHARED | MAP_NORESERVE, src, offset * PAGE);	//everytime we process a chunk of data

  //sm = mmap (0,share*PAGE, PROT_READ, MAP_SHARED | MAP_NORESERVE,src, offsetmytask*share*PAGE); //last time we process the rest

  if (MAP_FAILED == sm)
    {
      perror (" mmap source ");
      exit (EXIT_FAILURE);
    }

  return sm;
}
*/
/*-------------------------------------*/
int gather (F_set *File_head, int mytask, int ntask)
{
  printf ("gathering...\n");
  if (mytask == 0)
  {
        // The master thread will need to receive all computations from all other threads.
  	MPI_Status status;
        // MPI_Recv(void *buf, int count, MPI_DAtatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
        // We need to go and receive the data from all other threads.
        // The arbitrary tag we choose is 1, for now.
  	int i = 0;
     	for (i = 1; i < ntask; i++)
      	{
		BIGCAST temp, temp2, temp3, temp4;
		MPI_Recv (&temp, 1, MPI_LONG_LONG_INT, i, 1, MPI_COMM_WORLD, &status);
	  	MPI_Recv (&temp2, 2, MPI_LONG_LONG_INT, i, 1, MPI_COMM_WORLD, &status);
	  	MPI_Recv (&temp3, 3, MPI_LONG_LONG_INT, i, 1, MPI_COMM_WORLD, &status);
		MPI_Recv (&temp4, 4, MPI_LONG_LONG_INT, i, 1, MPI_COMM_WORLD, &status);
	  	printf ("RECEIVED %lld from thread %d\n", temp, i);
	  	File_head->reads_num += temp;
	  	File_head->reads_contam += temp2;
	  	File_head->hits += temp3;
		File_head->all_k+=temp4;
		
      	}
  }
  else
  {
      	// We are finished with the results in this thread, and need to send the data to thread 1.
      	// MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
      	// The destination is thread 0, and the arbitrary tag we choose for now is 1.
      	MPI_Send (File_head->reads_num, 1, MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD);
      	MPI_Send (File_head->reads_contam, 2, MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD);
      	MPI_Send (File_head->hits, 3, MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send (File_head->all_k, 4, MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD);
  }
  return 1;
}

/*-------------------------------------*/
BIGCAST gz_mpi (gzFile zip, BIGCAST offset, BIGCAST left, char *data, char type)
{
  char *start=NULL, *end = NULL;
  BIGCAST complete = 0;
  gzseek (zip, offset, SEEK_SET);
  if (left>2*ONEG)
  {
  	gzread (zip, data, 2*ONEG);
  	complete = 2*ONEG;
  }
  else
  {
	gzread (zip, data, left);
  	complete = left;
  }
  printf("here\n");
  if (type == '@')
  {
	if (offset!=0)
	{
		start = strstr (data,"\n+");
		start = strchr (strchr(start+1,'\n')+1,'\n')+1;
        }
	printf("dick\n"); 
	end = strrstr (data, "\n+");
        end = bac_2_n (end - 1);
  }
  else
  {
	if (offset!=0)
	{
		start = strchr (data, '>');
	}
        end = strrchr (data, '>') - 1;
  }
  if (start)
  {
	complete -= (start - data);
  }
  if (end)
  {
        complete += (end - data);
        memset (end, 0, strlen (end));
  }
  return complete;
}
