#include <stdlib.h>
#include <string.h>

#include "dtype.h"

model* malloc_model(int nsite, int nbond, int mhnspin) {
    model* m = (model*)malloc(sizeof(model));

    m->bond2type   = (int*)malloc(sizeof(int)*nbond);
    m->bond2hNspin = (int*)malloc(sizeof(int)*nbond);
    m->bond2weight = (double*)malloc(sizeof(double)*nbond);
    m->bond2index  = (int*)malloc(sizeof(int)*nbond*mhnspin);
    m->link        = (int*)malloc(sizeof(int)*mhnspin*4*20);
    m->insert      = (insert_rule*)malloc(sizeof(insert_rule)*20);

    m->nsite = nsite;
    m->nbond = nbond;
    m->mhnspin = mhnspin;
    m->sweight = 0;

    return m;
}

void free_model(model* m) {
    free(m->bond2type);
    free(m->bond2hNspin);
    free(m->bond2weight);
    free(m->bond2index);
    free(m->link);
    free(m->insert);
}

void copy_vertex(vertex* dist, vertex* src) {
    dist->tau = src->tau;
    dist->bond = src->bond;
    dist->hnspin = src->hnspin;
    for(int i=0;i<2*src->hnspin;i++) {
        dist->state[i] = src->state[i];
    }
}

world_line* malloc_world_line(int length, int mnspin, int nsite) {
    world_line* w = (world_line*)malloc(sizeof(world_line));

    w->sequenceA = (vertex*)malloc(sizeof(vertex)*length);
    w->sequenceB = (vertex*)malloc(sizeof(vertex)*length);
    w->cluster  = (int*)malloc(sizeof(int)*length*mnspin);
    w->weight   = (int*)malloc(sizeof(int)*length*mnspin);
    w->state    = (int*)malloc(sizeof(int)*nsite);
    w->first    = (int*)malloc(sizeof(int)*nsite);
    w->last     = (int*)malloc(sizeof(int)*nsite);

    w->length = length;
    w->nvertices = 0;
    w->mnspin = mnspin;
    w->flag = 0;
    w->nsite = nsite;

    return w;
}

void free_world_line(world_line* w) {
    free(w->sequenceA);
    free(w->sequenceB);
    free(w->cluster);
    free(w->weight);
    free(w->state);
    free(w->first);
    free(w->last);
}

void realloc_world_line(world_line* w, int length) {
    if(length > (w->length)) {
        vertex* sequence;
        sequence = (vertex*)malloc(sizeof(vertex)*length);
        for(int i=0;i<w->nvertices;i++) 
            copy_vertex(&(sequence[i]),&(w->sequenceA[i]));

        free(w->sequenceA);
        w->sequenceA = sequence;

        sequence = (vertex*)malloc(sizeof(vertex)*length);
        for(int i=0;i<w->nvertices;i++) 
            copy_vertex(&(sequence[i]),&(w->sequenceB[i]));

        free(w->sequenceB);
        w->sequenceB = sequence;

        free(w->cluster);
        free(w->weight);
        w->cluster = (int*)malloc(sizeof(int)*length*(w->mnspin));
        w->weight  = (int*)malloc(sizeof(int)*length*(w->mnspin));

        w->length = length;
    }
}

estimator* malloc_estimator(int length, char* name) {
    estimator* e = (estimator*)malloc(sizeof(estimator));

    e->sample = (double*)malloc(sizeof(double)*length);
    e->block  = (double*)malloc(sizeof(double)*1024);
    strcpy(e->name,name);
    e->length = length;
    e->n = 0;
    e->bsize  = 2;
    e->nblock = 0;

    return e;
}

void free_estimator(estimator* e) {
    free(e->sample);
    free(e->block);
}

void realloc_estimator(estimator* e, int length) {
    if(length>(e->length)) {
        double* sample = (double*)malloc(sizeof(double)*length);
        for(int i=0;i<(e->n);i++) 
            sample[i] = e->sample[i];

        free(e->sample);
        e->sample = sample;
    }
}
