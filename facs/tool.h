#ifndef _TOOL
#define _TOOL

#include <stdlib.h>
#include "bloom.h"

int fq_read_length (char *data);
extern int dx_add (int k_mer);
extern int get_parainfo (char *full, Queue *head);
extern char *fastq_relocate (char *data, int offset, int length);
extern char *jump (char *target, int type, float sampling_rate);
int read_full_check (bloom * bl, char* begin, int length,  char model, float tole_rate, F_set *File_head);
extern int query_read(char *begin, int length, char model, bloom * bl, float tole_rate, F_set *File_head);
void isodate(char* buf);
char* report(F_set * File_head, char* query, char* fmt, char* prefix);
char* substr(const char* str, size_t begin, size_t len);

#endif
