#include <math.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "union_find.h"

int weighted_sampling(double* cmf, int length, gsl_rng* rng) {
    double dis = gsl_rng_uniform_pos(rng)*cmf[length-1];

    int i = -1;
    int j = length;
    int d = j-1;

    while(d>1) {
        if(dis>cmf[i+d/2]){
            i = i+d/2;
        } else {
            j = i+d/2;
        }
        d = j-i;
    }

    return j;
}

int uniform_sequence_sampling(double** seq, int* size, double lam, double start, gsl_rng* rng) {
    double k=0;
    int n=0;

    double dis = gsl_rng_uniform_pos(rng);
    k -= log(dis)/lam;
    while((k<1.0) && (n<(*size))) {
        (*seq)[n] = k+start;
        dis = gsl_rng_uniform_pos(rng);
        k -= log(dis)/lam;
        n++;
    }
    while(k<1.0) {
        dis = gsl_rng_uniform_pos(rng);
        k -= log(dis)/lam;
        n++;
    }
    if(n>(*size)) {
        int ng = *size;
        double* temp = (double*)malloc(sizeof(double)*n*2);
        for(int i=0;i<(*size);i++) temp[i] = (*seq)[i];
        free(*seq);
        *seq = temp;
        *size = n*2;

        return ng;
    }

    return n;
}

void remove_vertices(world_line* w) {
    vertex* v;

    vertex* sequence1 = w->sequenceB;
    vertex* sequence2 = w->sequenceA;
    if(w->flag) {
        sequence1 = w->sequenceA;
        sequence2 = w->sequenceB;
    }

    int check_delete;
    int i,j,k;
    k=0;
    for(i=0;i<w->nvertices;i++) {
        v = &(sequence1[i]);
        check_delete = 1;

        for(j=0;j<(v->hNspin);j++) {
            if(v->state[j]!=v->state[j+v->hNspin])
                check_delete = 0;
        }

        if(!check_delete) {
            copy_vertex(&(sequence2[k]),&(sequence1[i]));
            k++;
        }
    }
    w->nvertices = k;
    w->flag = !(w->flag);
}

double* insert_seq;
int insert_length=0;
int insert_size=0;
void insert_vertices(world_line* w, model* m, gsl_rng* rng) {
    double lam = (m->sweight)*(w->beta);

    if(insert_size==0) {
        insert_size = (int)(lam+sqrt(lam)*10+1024);
        insert_seq = (double*)malloc(sizeof(double)*insert_size);
    }

    insert_length = uniform_sequence_sampling(&insert_seq,&insert_size,lam,0,rng);
    int length = insert_length+w->nvertices+1024;
    realloc_world_line(w,length);

    vertex* sequence1 = w->sequenceB;
    vertex* sequence2 = w->sequenceA;
    if(w->flag) {
        sequence1 = w->sequenceA;
        sequence2 = w->sequenceB;
    }

    int* pstate = w->pstate;
    int  nsite = w->nsite;

    for(int i=0;i<nsite;i++) pstate[i] = w->istate[i];

    vertex* v;
    int n,i,k,i_site,index;
    double tau1,tau2;

    int mhnspin = m->mhnspin;
    int lstate[mhnspin];

    for(i=0;i<insert_length;i++) {
        v = &(sequence1[k]);
        tau1 = v->tau;
        tau2 = insert_seq[i];

        while((tau1<tau2) && (k<(w->nvertices))) {
            for(i_site=0;i_site<(v->hNspin);i_site++) {
                index = m->bond2index[v->bond*mhnspin+i_site];
                pstate[index] = v->state[v->hNspin+i_site];
            }

            copy_vertex(&(sequence2[n]),v);
            n++;
            k++;
            v = &(sequence1[k]);
            tau1 = v->tau;
        }

        if(tau1!=tau2) {
            int bond     = weighted_sampling(m->cmf,m->nbond,rng);
            int t        = m->bond2type[bond];
            int* indices = &(m->bond2type[bond*mhnspin]);
            int hNspin   = m->bond2hNspin[bond];
            insert_rule rule = m->insert[t];

            for(i_site=0;i_site<hNspin;i_site++) 
                lstate[i_site] = pstate[indices[i_site]];

            if(rule(lstate)) {
                (sequence2[n]).tau    = tau2;
                (sequence2[n]).bond   = bond;
                (sequence2[n]).hNspin = hNspin;

                for(i_site=0;i_site<hNspin;i_site++) {
                    (sequence2[n]).state[i_site]        = lstate[i_site];
                    (sequence2[n]).state[i_site+hNspin] = lstate[i_site];
                }

                n++;
            }
        }
    }

    while(k<(w->nvertices)) {
        copy_vertex(&(sequence2[n]),&(sequence1[k]));
        n++;
        k++;
    }

    w->nvertices = n;
    w->flag = !(w->flag);
}
