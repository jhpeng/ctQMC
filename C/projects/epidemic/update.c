#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "union_find.h"

//#define weighted_sampling

#ifdef weighted_sampling
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

    return j;
}
#endif

static double* insert_seq;
static int* insert_bond;
static int insert_len=0;
static int insert_cap=0;
static double d_ave=-1;
static double d_max=-1;
static void uniform_sequence_sampling(model* m, double lam, double start, gsl_rng* rng) {
    if(d_max<0){
        d_max = 0;
        d_ave = 0;
        for(int i=0;i<(m->nbond);i++) {
            d_ave += m->bond2weight[i];
            if((m->bond2weight[i])>d_max)
                d_max = m->bond2weight[i];
        }
        d_ave = d_ave/(m->nbond);
    }

    if(insert_cap==0) {
        insert_cap = (int)(lam+sqrt(lam)*10+1024);
        insert_seq  = (double*)malloc(sizeof(double)*insert_cap);
        insert_bond = (int*)malloc(sizeof(int)*insert_cap);
    }

    lam = lam*d_max/d_ave;

    double k=0;
    int n=0;
    int bond;

    double dis = gsl_rng_uniform_pos(rng);
    k -= log(dis)/lam;
    while((k<1.0) && (n<insert_cap)) {
        bond = gsl_rng_uniform_pos(rng)*(m->nbond);
        if(gsl_rng_uniform_pos(rng)*d_max<(m->bond2weight[bond])){
            insert_seq[n]  = k+start;
            insert_bond[n] = bond;
            n++;
        }

        dis = gsl_rng_uniform_pos(rng);
        k -= log(dis)/lam;
    }
    while(k<1.0) {
        dis = gsl_rng_uniform_pos(rng);
        k -= log(dis)/lam;
        n++;
    }
    insert_len = n;
    if(n>insert_cap) {
        double* seq = (double*)malloc(sizeof(double)*n*2);
        int* bond   = (int*)malloc(sizeof(int)*n*2);
        for(int i=0;i<insert_cap;i++){
            seq[i]  = insert_seq[i];
            bond[i] = insert_bond[i];
        }
        free(insert_seq);
        insert_seq  = seq;
        insert_bond = bond;
        insert_len = insert_cap;
        insert_cap = n*2;
    }
}

static int ninfection_count=0;
static int nrecover_count=0;
static double ninfection_ave=0;
static double nrecover_ave=0;

double ninfection_ave_value() {
    return ninfection_ave/ninfection_count;
}

double nrecover_ave_value() {
    return nrecover_ave/nrecover_count;
}

void print_ninfection() {
    ninfection_ave = ninfection_ave/ninfection_count;
    printf("# of infection = %.12e\n",ninfection_ave);
    ninfection_ave = 0;
    ninfection_count=0;
}

void print_nrecover() {
    nrecover_ave = nrecover_ave/nrecover_count;
    printf("# of recover = %.12e\n",nrecover_ave);
    nrecover_ave = 0;
    nrecover_count=0;
}

void remove_vertices(world_line* w) {
    vertex* v;

    vertex* sequence1 = w->sequenceB;
    vertex* sequence2 = w->sequenceA;
    if(w->flag) {
        sequence1 = w->sequenceA;
        sequence2 = w->sequenceB;
    }

    int ninfection=0;
    int nrecover=0;
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

            if(v->hNspin==1) {
                nrecover++;
            } else if(v->hNspin==2) {
                ninfection++;
            }
        }
    }
    w->nvertices = k;
    w->flag = !(w->flag);

    ninfection_ave += (double)ninfection;
    nrecover_ave += (double)nrecover;

    ninfection_count++;
    nrecover_count++;
}

void swapping_graphs(world_line* w, model* m, gsl_rng* rng) {
    int nnode = m->nsite;
    int nedge = (m->nbond-3*nnode)/7;
    vertex* v;

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    for(int i=0;i<(w->nvertices);i++) {
        v = &(sequence[i]);

        // swap between type (1,3,6) or (2,4,6)
        int* state = v->state;
        int type = m->bond2type[v->bond];
        if((type==1 || type==3) || (type==6 && state[0]==1)) {
            int i_edge = (v->bond)%nedge;
            double dis = gsl_rng_uniform_pos(rng)*3;
            if(dis<1) {
                v->bond = 1*nedge+i_edge;
            } else if(dis<2) {
                v->bond = 3*nedge+i_edge;
            } else {
                v->bond = 6*nedge+i_edge;
            }
        } else if((type==2 || type==4) || type==6) {
            int i_edge = (v->bond)%nedge;
            double dis = gsl_rng_uniform_pos(rng)*3;
            if(dis<1) {
                v->bond = 2*nedge+i_edge;
            } else if(dis<2) {
                v->bond = 4*nedge+i_edge;
            } else {
                v->bond = 6*nedge+i_edge;
            }
        }

        // swap between type 7&8
        else if(type==7 || type==8) {
            int i_node = (v->bond-7*nedge)%nnode;
            if(gsl_rng_uniform_pos(rng)<0.5) {
                v->bond = 7*nedge+i_node;
            } else {
                v->bond = 7*nedge+nnode+i_node;
            }
        }
    }
}

void insert_vertices(world_line* w, model* m, gsl_rng* rng) {
    double lam = (m->sweight)*(w->beta);

    uniform_sequence_sampling(m,lam,0,rng);

    int length = insert_cap+w->nvertices;
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

    for(i=0;i<insert_len;i++) {
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
            int bond     = insert_bond[i];
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

/*  disable for open boundary
**  for(i=0;i<nsite;i++) {
**      if(first[i]!=-1) {
**          merge(w->cluster,w->weight,first[i],last[i]);
**      }
**  }
*/
}


static int ncluster_flippable_count=0;
static double ncluster_flippable_ave=0;
double ncluster_flippable_ave_value() {
    return ncluster_flippable_ave/ncluster_flippable_count;
}

void print_ncluster_flippable() {
    ncluster_flippable_ave = ncluster_flippable_ave/ncluster_flippable_count;
    printf("# of flippable cluster = %.12e \n",ncluster_flippable_ave);
    ncluster_flippable_ave=0;
    ncluster_flippable_count=0;
}

void flip_cluster(world_line* w, gsl_rng* rng) {
    int* state;
    int hNspin,idv,idr,id,p,i,j;

    int mnspin = w->mnspin;
    int nsite  = w->nsite;
    int ncluster_flippable=0;

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
                ncluster_flippable++;
                if(gsl_rng_uniform_pos(rng)<0.5) {
                    w->weight[idr] =  0;
                } else {
                    w->weight[idr] = -1;
                }
            }
            if(w->weight[idr]==0) {
                state[j] = -state[j];
            }
        }
    }
    ncluster_flippable_ave += (double)ncluster_flippable;
    ncluster_flippable_count++;
    //printf("%d \n",ncluster_flippable);

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
