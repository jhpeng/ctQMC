#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "norm.h"
#include "dtype.h"


conf* malloc_conf(int size, int nsite) {
    int* length = (int*)malloc(sizeof(int)*nsite);
    int* sigma = (int*)malloc(sizeof(int)*nsite*size);
    double* tau = (double*)malloc(sizeof(int)*nsite*size);
    conf* c = (conf*)malloc(sizeof(conf));
    
    if(length==NULL) {
       printf("memory allocate error : malloc_conf\n");
       exit(-1); 
    }

    if(sigma==NULL) {
       printf("memory allocate error : malloc_conf\n");
       exit(-1); 
    }

    if(tau==NULL) {
        printf("mempry allocate error : malloc_conf\n");
        exit(-1);
    }

    if(c==NULL) {
        printf("mempry allocate error : malloc_conf\n");
        exit(-1);
    }
    
    for(int i=0; i<nsite; i++) {
        length[i]=0;
    }

    c->length = length;
    c->sigma  = sigma;
    c->tau    = tau;
    c->size   = size;
    c->nsite  = nsite;

    return c;
}

void free_conf(conf* c) {
    free(c->sigma);
    free(c->tau);
    free(c);
}

void realloc_conf(conf* c, int size_new) {
    if(c->size >= size_new) {
        return;
    }

    int nsite   = c->nsite;
    int size    = c->size;
    int* sigma  = c->sigma;
    double* tau = c->tau;

    int*  sigma_new = (int*)malloc(sizeof(int)*nsite*size_new);
    double* tau_new = (double*)malloc(sizeof(double)*nsite*size_new);

    if(sigma_new == NULL) {
       printf("memory allocate error : realloc_conf\n");
       exit(-1);
    }

    if(tau_new == NULL) {
        printf("mempry allocate error : realloc_conf\n");
        exit(-1);
    }

    for(int i_site = 0; i_site < nsite; i_site++) {
        for(int i = 0; i < size; i++) {
            sigma_new[i_site*size_new+i] = sigma[i_site*size+i];
            tau_new[i_site*size_new+i] = tau[i_site*size+i];
        }
    }

    free(sigma);
    free(tau);

    c->sigma = sigma_new;
    c->tau = tau_new;
    c->size = size_new;
}

void show_conf(conf* c) {

}

void world_line_to_conf(conf* c, world_line* w, model* m) {
    vertex* v;
    int bond, hNspin, type, index, size;
    int* indices;
    int* state;
    double tau;

    int nsite  = w->nsite;
    int mnspin = w->mnspin;

    vertex* sequence = w->sequenceB;
    if(w->flag)
        sequence = w->sequenceA;

    for(int i=0; i<nsite; i++) {
        c->length[i] = 0;
    }

    for(int i=0; i<(w->nvertices); i++) {
        v       = &(sequence[i]);
        bond    = v->bond;
        hNspin  = v->hNspin;
        tau     = v->tau;
        state   = v->state;
        indices = &(m->bond2index[bond*(m->mhnspin)]);

        for(int j=0; j<hNspin; j++) {
            if(state[j] != state[j+hNspin]) {
                size = c->size;
                int length = c->length[indices[j]];

                if((length+1)>size) {
                    realloc_conf(c, size+1000);
                    size = c->size;
                }

                c->sigma[size*indices[j]+length] = state[j+hNspin];
                c->length[indices[j]] = length+1;
            }
        }    
    }
}

