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
    m->cmf         = (double*)malloc(sizeof(double)*nbond);

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
    free(m->cmf);
    free(m);
}

void copy_vertex(vertex* dist, vertex* src) {
    dist->tau = src->tau;
    dist->bond = src->bond;
    dist->hNspin = src->hNspin;
    for(int i=0;i<2*src->hNspin;i++) {
        dist->state[i] = src->state[i];
    }
}

world_line* malloc_world_line(int length, int mnspin, int nsite) {
    world_line* w = (world_line*)malloc(sizeof(world_line));

    w->sequenceA = (vertex*)malloc(sizeof(vertex)*length);
    w->sequenceB = (vertex*)malloc(sizeof(vertex)*length);
    w->cluster  = (int*)malloc(sizeof(int)*length*mnspin);
    w->weight   = (int*)malloc(sizeof(int)*length*mnspin);
    w->istate    = (int*)malloc(sizeof(int)*nsite);
    w->pstate    = (int*)malloc(sizeof(int)*nsite);
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
    free(w->istate);
    free(w->pstate);
    free(w->first);
    free(w->last);
    free(w);
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

estimator* malloc_estimator(int length, char name[128]) {
    estimator* e = (estimator*)malloc(sizeof(estimator));

    e->samples = (double*)malloc(sizeof(double)*length);
    e->blocks  = (double*)malloc(sizeof(double)*1024);
    strcpy(e->name,name);
    e->length = length;
    e->n = 0;
    e->bsize  = 2;
    e->nblock = 0;

    return e;
}

void free_estimator(estimator* e) {
    free(e->samples);
    free(e->blocks);
    free(e);
}

void realloc_estimator(estimator* e, int length) {
    if(length>(e->length)) {
        double* samples = (double*)malloc(sizeof(double)*length);
        for(int i=0;i<(e->n);i++) 
            samples[i] = e->samples[i];

        free(e->samples);
        e->samples = samples;
        e->length = length;
    }
}
