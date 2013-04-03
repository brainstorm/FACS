#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libgen.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>

#include "tool.h"
#include "check.h"
#include "bloom.h"
#include "file_dir.h"

#ifndef __clang__ 
#include <omp.h>
#endif

void
report(char *detail, char *filename, F_set * File_head, char* query,
       char* fmt, char* prefix)
{
  char buffer[200] = { 0 };
  float contamination_rate = (float) (File_head->reads_contam) /
                             (float) (File_head->reads_num);


  if(!fmt){
      return;
  // JSON output format (via stdout)
  } else if(!strcmp(fmt, "json")) {
      isodate(buffer);

      printf("{\n");
      printf("\t\"timestamp\": \"%s\"\n", buffer);
      printf("\t\"sample\": \"%s\"\n", basename(query)); //sample (query)
      printf("\t\"bloom_filter\": \"%s\"\n", basename(filename)); //reference
      printf("\t\"total_read_count\": %lld,\n", File_head->reads_num);
      printf("\t\"contaminated_reads\": %lld,\n", File_head->reads_contam);
      printf("\t\"total_hits\": %lld,\n", File_head->hits);
      printf("\t\"contamination_rate\": %f,\n", contamination_rate);
      printf("}\n");

  // TSV output format (via file in CWD)
  } else if (!strcmp(fmt, "tsv")) {
      strcat (detail, "sample\tbloom_filter\ttotal_read_count\t\
contaminated_reads\tcontamination_rate\n");

      sprintf(buffer, "%s\t%s\t%lld\t%lld\t%f\n", basename(query),
              basename(filename), File_head->reads_num,
              File_head->reads_contam, contamination_rate);
      strcat(detail, buffer);

      write_result(strcat(basename(query), ".tsv"), detail);
  }
}
