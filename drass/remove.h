#ifndef _REMOVE_CONTAM_L
#define _REMOVE_CONTAM_L

#include "bloom.h"
extern void fasta_process_m (bloom * bl, Queue * info, Queue * tail, float tole_rate, F_set *File_head);
extern void fastq_process_m (bloom * bl, Queue * info, Queue * tail, float tole_rate, F_set *File_head);
extern int remove_main(float tole_rate, char *source, char *ref, char *list, char *prefix, int help);
extern void save_result (char *source, char *obj_file, int type, char *prefix, char *clean, char *clean2, char *contam, char *contam2);

#endif
