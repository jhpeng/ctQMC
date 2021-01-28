#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "union_find.h"

//#define check_weighted_sampling

#ifdef check_weighted_sampling
int* count = NULL;
#endif

static int weighted_sampling(double* cmf, int length, gsl_rng* rng) {
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

#ifdef check_weighted_sampling
    if(count==NULL){
        count =  (int*)malloc(sizeof(int)*length);
        for(int k=0;k<length;k++) count[k]=0;
    }
    count[j]++;
#endif

    return j;
}

static int uniform_sequence_sampling(double** seq, int* size, double lam, double start, gsl_rng* rng) {
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

    k=0;
    n=0;

    tau1 = 0;
    if(w->nvertices!=0) tau1 = (sequence1[0]).tau;

    for(i=0;i<insert_length;i++) {
        tau2 = insert_seq[i];

        while((tau1<tau2) && (k<(w->nvertices))) {
            v = &(sequence1[k]);
            for(i_site=0;i_site<(v->hNspin);i_site++) {
                index = m->bond2index[v->bond*mhnspin+i_site];
                pstate[index] = v->state[v->hNspin+i_site];
            }

            copy_vertex(&(sequence2[n]),v);
            n++;
            k++;

            if(k<(w->nvertices)){
                tau1 = (sequence1[k]).tau;
            }
        }

        if(tau1!=tau2) {
            int bond     = weighted_sampling(m->cmf,m->nbond,rng);
            int t        = m->bond2type[bond];
            int hNspin   = m->bond2hNspin[bond];
            insert_rule rule = m->insert[t];

            for(i_site=0;i_site<hNspin;i_site++) { 
                index = m->bond2index[bond*mhnspin+i_site];
                lstate[i_site] = pstate[index];
            }

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

#ifdef check_weighted_sampling
    double tcount=0;
    for(i=0;i<(m->nbond);i++)
        tcount += count[i];

    for(i=0;i<(m->nbond);i++){
        double weight = count[i]/tcount*(m->sweight);
        printf("%.4f ",weight);
    }

    printf("\n");
#endif

}

void clustering(world_line* w, model* m) {
    vertex* v;
    int bond,hNspin,t,idn,idp;
    int i,j,index;
    int mnspin = w->mnspin;
    int nsite = w->nsite;
    int* rule;
    int* indices;

    int* first = w->first;
    int* last  = w->last;

    for(i=0;i<nsite;i++) {
        last[i]  = -1;
        first[i] = -1;
    }

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    for(i=0;i<(w->nvertices);i++) {
        v       = &(sequence[i]);
        bond    = v->bond;
        hNspin  = v->hNspin;
        t       = m->bond2type[bond];
        rule    = &(m->link[4*(m->mhnspin)*t]);
        indices = &(m->bond2index[bond*(m->mhnspin)]);

        for(j=0;j<2*hNspin;j++) {
            idn = i*mnspin+j;
            idp = i*mnspin+rule[j];
            w->cluster[idn] = idp;
            w->weight[idn]  = rule[2*hNspin+j];
        }

        for(j=0;j<hNspin;j++) {
            index = indices[j];
            idp = i*mnspin+j;
            idn = i*mnspin+j+hNspin;
            if(first[index]==-1) {
                first[index] = idp;
                last[index]  = idn;
            } else {
                merge(w->cluster,w->weight,last[index],idp);
                last[index] = idn;
            }
        }
    }

    for(i=0;i<nsite;i++) {
        if(first[i]!=-1) {
            merge(w->cluster,w->weight,first[i],last[i]);
        }
    }
}

void flip_cluster(world_line* w, gsl_rng* rng) {
    int* state;
    int hNspin,idv,idr,id,p,i,j;

    int mnspin = w->mnspin;
    int nsite  = w->nsite;

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    for(i=0;i<(w->nvertices);i++) {
        hNspin = (sequence[i]).hNspin;
        state  = (sequence[i]).state;

        for(j=0;j<2*hNspin;j++) {
            idv = i*mnspin+j;
            idr = root(w->cluster,idv);
            if(w->weight[idr]>0) {
                if(gsl_rng_uniform_pos(rng)<0.5) {
                    w->weight[idr] =  0;
                } else {
                    w->weight[idr] = -1;
                }
            }

            state[j] = state[j]*((w->weight[idr])*2+1);
        }
    }

    for(i=0;i<nsite;i++) {
        id = w->first[i];
        if(id!=-1) {
            p = id/mnspin;
            j  =id%mnspin;
            w->istate[i] = (sequence[p]).state[j];
        } else if(gsl_rng_uniform_pos(rng)<0.5) {
            w->istate[i] =  1;
        } else {
            w->istate[i] = -1;
        }
    }
}

int check_periodic(world_line* w, model* m) {
    vertex* sequence = w->sequenceB;
    if(w->flag)
        sequence = w->sequenceA;

    for(int i=0;i<(w->nsite);i++)
        w->pstate[i] = w->istate[i];

    vertex* v;
    for(int i=0;i<(w->nvertices);i++) {
        v = &(sequence[i]);
        int bond   = v->bond;
        int hNspin = v->hNspin;

        for(int j=0;j<hNspin;j++) {
            int i_site = m->bond2index[bond*(m->mhnspin)+j];
            if(w->pstate[i_site]!=(v->state[j])){
                printf("something wrong!\n");
            }
            w->pstate[i_site] = v->state[j+hNspin];
        }
    }

    int check=1;
    for(int i=0;i<(w->nsite);i++){
        if((w->pstate[i])!=(w->istate[i])) {
            check=0;
        }
    }

    return check;
}
