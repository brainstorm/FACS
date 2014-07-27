/* SWIG interface for FACS */

%module facs

%{
#include <Python.h>
#include "query.h"
%}

%init %{
    printf("Welcome to FACS\n");
%}
