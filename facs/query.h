#ifndef _QUERY
#define _QUERY

#include "tool.h"
#include "bloom.h"
#include <zlib.h>
extern char* bq_main (int argc, char **argv);
extern char* query (char *query, char *bloom_filter, double tole_rate,
                  double sampling_rate, char *list, char *target_path, char* report_fmt);
#endif
