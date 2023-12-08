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
    int nsite = c->nsite;
    int size  = c->size;

    printf("--------------- show configuration --------------\n");
    printf("nsite = %d \n",nsite);

    printf("------------------------------\n");
    for(int i=0; i<nsite; i++) {
        int length = c->length[i];
        if(i<10) {
            printf("%d    | ",i);
        } else if(i<100) {
            printf("%d   | ",i);
        } else if(i<1000) {
            printf("%d  | ",i);
        } else {
            printf("%d | ",i);
        }

        for(int j=0; j<length; j++) {
            int sigma = c->sigma[size*i+j];
            double tau = c->tau[size*i+j];

            //printf("%d ",sigma);
            if(sigma==1) {
                printf("(+) ");
            } else {
                printf("(-) ");
            }
            printf("%.3f | ",tau);
        }
        printf("\n");
    }
}

void inner_product(conf* c1, conf* c2) {
    if(c1->nsite != c2->nsite) {
        printf("c1 and c2 should have the same nsite!\n");
        return;
    }

    int nsite = c1->nsite;
    double z=0;

    int length1, length2;
    int size1, size2;
    int sigma1, sigma2;
    double tau1, tau2;
    for(int i=0; i<nsite; i++) {
        length1 = c1->length[i];
        length2 = c2->length[i];

        size1 = c1->size;
        size2 = c2->size;

        sigma1 = c1->sigma[size1*i+0];
        sigma2 = c2->sigma[size2*i+0];

        tau1 = c1->tau[size1*i+0];
        tau2 = c2->tau[size2*i+0];
    }
}

void world_line_to_conf(conf* c, world_line* w, model* m) {
    vertex* v;
    int bond, hNspin, size;
    //int index, type;
    int* indices;
    int* state;
    double tau;

    int nsite  = w->nsite;
    //int mnspin = w->mnspin;

    vertex* sequence = w->sequenceB;
    if(w->flag)
        sequence = w->sequenceA;

    for(int i=0; i<nsite; i++) {
        size = c->size;

        c->sigma[size*i+0] = w->istate[i];
        c->tau[size*i+0] = 0;
        c->length[i] = 1;
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
                c->tau[size*indices[j]+length] = tau;
                c->length[indices[j]] = length+1;
            }
        }    
    }
}

