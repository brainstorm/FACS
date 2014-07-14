%module facs

%{
#include <Python.h>
#include "query.h"
%}

char *query (char *query, char *bloom_filter, double tole_rate, double sampling_rate, char *list, char *target_path, char *report_fmt, char mode);
