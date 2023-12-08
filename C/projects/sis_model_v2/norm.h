#ifndef norm_h
#define norm_h

#include "dtype.h"

typedef struct conf {
    int* length;
    int* sigma;
    double* tau;
    int size;
    int nsite;
} conf;

conf* malloc_conf(int size, int nsite);

void free_conf(conf* c);

void realloc_conf(conf* c, int size_new);

void show_conf(conf* c);

void world_line_to_conf(conf* c, world_line* w, model* m);

#endif
